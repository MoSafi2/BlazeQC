"""CG content: composable Collector, Summarizer (with plotting), assembled CGModule."""

from collections.list import List
from math import sqrt, exp, pi
from python import Python, PythonObject
from blazeseq import FastqRecord, RefRecord
from blazeqc.stats.traits import (
    Collector,
    Summarizer,
    PlotOutput,
    FastqcHtmlOutput,
    FastqcDataOutput,
)
from blazeqc.stats.summary_utils import SummaryContext, GradeEntry, DefaultOutputter
from blazeqc.helpers import tensor_to_numpy_1d, list_float64_to_numpy
from blazeqc.html_maker import result_panel
from blazeqc.limits import GC_SEQUENCE_WARN, GC_SEQUENCE_ERROR


# ----- Collector: tally only -----

struct CGCollector(Collector, Copyable, Movable):
    """GC Collector: collect GC content counts per read."""
    var cg_content: List[Int64]

    fn __init__(out self) raises:
        self.cg_content = List[Int64](capacity=101)
        for _ in range(101):
            self.cg_content.append(0)

    # TODO: Optimize this to use a SIMD bitmask instead of a loop.
    @always_inline
    fn tally_read(mut self, record: FastqRecord):
        if len(record) == 0:
            return
        var cg_num = 0
        var seq_span = record.sequence()
        for index in range(0, len(record)):
            if (
                seq_span[index] & 0b111 == 3
                or seq_span[index] & 0b111 == 7
            ):
                cg_num += 1
        var read_cg_content = Int(
            round(cg_num * 100 / Int(len(record)))
        )
        self.cg_content[read_cg_content] += 1

    @always_inline
    fn tally_read(mut self, record: RefRecord):
        if len(record) == 0:
            return
        var cg_num = 0
        var seq_span = record.sequence().as_bytes()
        for index in range(0, len(record)):
            if (
                seq_span[index] & 0b111 == 3
                or seq_span[index] & 0b111 == 7
            ):
                cg_num += 1
        var read_cg_content = Int(
            round(cg_num * 100 / Int(len(record)))
        )
        self.cg_content[read_cg_content] += 1


# ----- Summarizer: prepare + grades + data for output + to_plot -----

struct CGSummarizer(Summarizer, PlotOutput, Copyable, Movable):
    """Summarization and plotting in one struct. Feed collector then prepare(ctx)."""
    var _cache_theoretical: List[Float64]
    var _cache_cg_content: List[Int64]
    var _cache_grade: String
    var _cache_ready: Bool

    fn __init__(out self) raises:
        self._cache_theoretical = List[Float64]()
        self._cache_cg_content = List[Int64]()
        self._cache_grade = ""
        self._cache_ready = False

    fn feed(mut self, collector: CGCollector):
        """Copy collector counts into cache for prepare()."""
        self._cache_cg_content = List[Int64](capacity=101)
        for i in range(len(collector.cg_content)):
            self._cache_cg_content.append(collector.cg_content[i])

    fn summerize(mut self, ctx: SummaryContext) raises:
        """Compute theoretical and grade from cached counts (call feed(collector) first)."""
        self._cache_theoretical = self._calculate_theoretical_distribution(self._cache_cg_content, ctx.num_reads)
        var max_dev = self._max_gc_deviation()
        if max_dev > GC_SEQUENCE_ERROR:
            self._cache_grade = "fail"
        elif max_dev > GC_SEQUENCE_WARN:
            self._cache_grade = "warn"
        else:
            self._cache_grade = "pass"
        self._cache_ready = True

    fn grade(self) -> GradeEntry:
        return GradeEntry("Per Sequence GC Content", self._cache_grade)

    fn data_block_body(self) -> String:
        """Module-specific lines (header + data); Outputter wraps with >>name\tgrade and >>END_MODULE."""
        if not self._cache_ready:
            return ""
        var body = "#GC Content\tCount\n"
        for i in range(len(self._cache_cg_content)):
            body += "{}\t{}\n".format(i, self._cache_cg_content[i])
        return body

    fn module_legend(self) -> String:
        return "Per Sequence GC Content"

    fn panel_id(self) -> String:
        return "cg_content"

    fn plot_result(self) raises -> PythonObject:
        """Build GC figure from cached data (summarizer is also the plotter)."""
        var plt = Python.import_module("matplotlib.pyplot")
        var arr = tensor_to_numpy_1d(self._cache_cg_content)
        var x = plt.subplots()
        var fig = x[0]
        var ax = x[1]
        ax.plot(arr, label="GC count per read")
        ax.plot(
            list_float64_to_numpy(self._cache_theoretical),
            label="Theoretical distribution",
        )
        ax.set_title("GC distribution over all sequences")
        ax.set_xlabel("Mean GC content (%)")
        return fig


    fn _calculate_theoretical_distribution(self, counts: List[Int64], total_counts: Int64) -> List[Float64]:
        """Compute a theoretical normal distribution fitted to the observed GC bin counts.

        Algorithm: 
        (1) Find the mode (GC bin with maximum count).
        (2) Compute weighted standard deviation: sum over bins of (bin - mode)^2 * count,
        divided by (total - 1).
        (3) Evaluate the normal PDF with that mean and stdev at
        each bin index, scaled by total count so the curve matches total reads. If
        total is 0 or stdev is 0, returns zeros or a spike at the mode respectively.
        """
        var n = len(counts)
        var result = List[Float64](length=n, fill=0.0)
        if total_counts == 0:
            return result^

        # Mode = bin index with highest count (location of the normal)
        var mode = _mode(counts, n)
        # Weighted standard deviation of GC bins around the mode
        var stdev = _weighted_stdev(counts, mode, total_counts)

        var total_f = Float64(total_counts)
        var mode_f = Float64(mode)
        # Degenerate case: all mass at mode
        if stdev == 0.0:
            result[mode] = total_f
            return result^

        # Use standalone normal PDF helper for non-degenerate case
        result = _normal_pdf(n, mode_f, stdev, total_f)
        return result^

    fn _max_gc_deviation(self) -> Float64:
        """Maximum absolute percentage-point deviation between observed and theoretical GC distribution.

        For each GC bin, observed and theoretical are expressed as percentages of their
        respective totals. The deviation is |obs_pct - theor_pct|. Returns the maximum
        over all bins; used to assign pass/warn/fail from GC_SEQUENCE_* limits.
        """
        var total_obs: Float64 = 0.0
        for i in range(len(self._cache_cg_content)):
            total_obs += Float64(self._cache_cg_content[i])

        var total_theor: Float64 = 0.0
        for i in range(len(self._cache_theoretical)):
            total_theor += self._cache_theoretical[i]

        if total_obs <= 0 or total_theor <= 0:
            return 0.0

        var max_dev: Float64 = 0.0
        for i in range(len(self._cache_cg_content)):
            var o_pct = (Float64(self._cache_cg_content[i]) / total_obs) * 100.0
            var t_val = self._cache_theoretical[i]
            var t_pct = (t_val / total_theor) * 100.0
            var dev = o_pct - t_pct
            if dev < 0:
                dev = -dev
            if dev > max_dev:
                max_dev = dev
        return max_dev



# ----- Assembled module: Collector + Summarizer + DefaultOutputter -----

struct CGModule(Collector, Summarizer, FastqcDataOutput, FastqcHtmlOutput, Copyable, Movable):
    """Module assembled from Collector + Summarizer; uses DefaultOutputter for text/HTML."""
    var collector: CGCollector
    var summarizer: CGSummarizer

    fn __init__(out self) raises:
        self.collector = CGCollector()
        self.summarizer = CGSummarizer()

    fn tally_read(mut self, record: FastqRecord):
        self.collector.tally_read(record)

    fn tally_read(mut self, record: RefRecord):
        self.collector.tally_read(record)

    fn prepare(mut self, ctx: SummaryContext) raises:
        self.summarizer.feed(self.collector)
        self.summarizer.summerize(ctx)

    fn summerize(mut self, ctx: SummaryContext) raises:
        """Trait-compatible alias for prepare(ctx)."""
        self.prepare(ctx)

    fn grade(self) raises -> GradeEntry:
        return self.summarizer.grade()

    fn grades(self) raises -> List[GradeEntry]:
        """Return list of grade entries (single-panel API)."""
        var out = List[GradeEntry]()
        out.append(self.summarizer.grade())
        return out^

    fn to_data_text(self, ctx: SummaryContext) raises -> String:
        """FastQC-style data block text for this module."""
        # ctx is currently unused but kept for trait compatibility.
        var body = self.summarizer.data_block_body()
        var g = self.summarizer.grade()
        var out = DefaultOutputter()
        return out.wrap_data_block(self.summarizer.module_legend(), g.grade, body)

    fn plot_result(self) raises -> PythonObject:
        return self.summarizer.plot_result()

    fn to_html(self) raises -> result_panel:
        var fig = self.summarizer.plot_result()
        var out = DefaultOutputter()
        return out.make_panel(
            self.summarizer.panel_id(),
            self.summarizer.grade().grade,
            self.summarizer.module_legend(),
            fig,
        )

    fn to_html_panels(self) raises -> List[result_panel]:
        """Return list of HTML panels (single GC content panel)."""
        var panels = List[result_panel]()
        panels.append(self.to_html())
        return panels^



@always_inline
fn _mode(counts: List[Int64], n: Int) -> Int:
    var max_count: Int64 = counts[0]
    var mode: Int = 0
    for i in range(1, n):
        if counts[i] > max_count:
            max_count = counts[i]
            mode = i
    return mode


@always_inline
fn _weighted_stdev(counts: List[Int64], mode: Int, total_counts: Int64) -> Float64:
    """Weighted standard deviation for histogram-like GC counts.

    Uses indices as bin centers and 'counts' as weights, computing:
        sqrt( sum_i (i - mode)^2 * counts[i] / (total_counts - 1) )
    Returns 0.0 when total_counts <= 1.
    """
    if total_counts <= 1:
        return 0.0
    var n = len(counts)
    var sum_sq: Float64 = 0.0
    for i in range(n):
        var diff = Float64(i) - Float64(mode)
        sum_sq += diff * diff * Float64(counts[i])
    return sqrt(sum_sq / Float64(total_counts - 1))


@always_inline
fn _normal_pdf(n: Int, mode: Float64, stdev: Float64, total: Float64) -> List[Float64]:
    """Evaluate a normal PDF at integer bin centers [0, n), centered at `mode`.

    Uses:
        (1 / (stdev * sqrt(2*pi))) * exp(-0.5 * z^2), scaled by `total`,
    where z = (x - mode) / stdev.
    """
    var result = List[Float64](length=n, fill=0.0)
    var scale = stdev * sqrt(2.0 * pi)
    for i in range(n):
        var x = Float64(i)
        var z = (x - mode) / stdev
        result[i] = (exp(-0.5 * z * z) / scale) * total
    return result^
