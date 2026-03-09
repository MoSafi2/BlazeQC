"""CG content: composable Collector, Summarizer (with plotting), assembled CGModule."""

from collections.list import List
from python import Python, PythonObject
from blazeseq import FastqRecord, RefRecord
from blazeqc.stats.analyser import Analyser
from blazeqc.stats.traits import (
    SummaryContext,
    GradeEntry,
    StatSummarizer,
    TextOutput,
    PlotOutput,
    HtmlOutput,
    DefaultOutputter,
)
from blazeqc.helpers import tensor_to_numpy_1d
from blazeqc.html_maker import result_panel
from blazeqc.limits import GC_SEQUENCE_WARN, GC_SEQUENCE_ERROR


# ----- Collector: tally only -----

struct CGCollector(Analyser, Copyable, Movable):
    """Collection only: raw cg_content counts. No prepare, grades, or output."""
    var cg_content: List[Int64]

    fn __init__(out self) raises:
        self.cg_content = List[Int64](capacity=101)
        for _ in range(101):
            self.cg_content.append(0)

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

struct CGSummarizer(StatSummarizer, Copyable, Movable):
    """Summarization and plotting in one struct. Feed collector then prepare(ctx)."""
    var _cache_theoretical: PythonObject
    var _cache_cg_content: List[Int64]
    var _cache_grade: String
    var _cache_ready: Bool

    fn __init__(out self) raises:
        self._cache_theoretical = Python.evaluate("None")
        self._cache_cg_content = List[Int64]()
        self._cache_grade = ""
        self._cache_ready = False

    fn feed(mut self, collector: CGCollector):
        """Copy collector counts into cache for prepare()."""
        self._cache_cg_content = List[Int64](capacity=101)
        for i in range(len(collector.cg_content)):
            self._cache_cg_content.append(collector.cg_content[i])

    fn _calculate_theoretical_distribution(self, counts: List[Int64]) raises -> PythonObject:
        var np = Python.import_module("numpy")
        var sc = Python.import_module("scipy")
        var arr = tensor_to_numpy_1d(counts)
        var total_counts = np.sum(arr)
        var x_categories = np.arange(len(arr))
        var mode = np.argmax(arr)
        var stdev = np.sqrt(
            np.sum((x_categories - mode) ** 2 * arr) / (total_counts - 1)
        )
        var nd = sc.stats.norm(loc=mode, scale=stdev)
        return nd.pdf(x_categories) * total_counts

    fn _max_gc_deviation(self) raises -> Float64:
        var np = Python.import_module("numpy")
        var obs_arr = tensor_to_numpy_1d(self._cache_cg_content)
        var total_obs = Float64(py=np.sum(obs_arr))
        var theor = self._cache_theoretical
        var total_theor = Float64(py=np.sum(theor))
        if total_obs <= 0 or total_theor <= 0:
            return 0.0
        var max_dev: Float64 = 0.0
        for i in range(len(self._cache_cg_content)):
            var o_pct = (Float64(self._cache_cg_content[i]) / total_obs) * 100.0
            var t_val = Float64(py=theor[Int(i)])
            var t_pct = (t_val / total_theor) * 100.0
            var dev = o_pct - t_pct
            if dev < 0:
                dev = -dev
            if dev > max_dev:
                max_dev = dev
        return max_dev

    fn prepare(mut self, ctx: SummaryContext) raises:
        """Compute theoretical and grade from cached counts (call feed(collector) first)."""
        self._cache_theoretical = self._calculate_theoretical_distribution(self._cache_cg_content)
        var max_dev = self._max_gc_deviation()
        if max_dev > GC_SEQUENCE_ERROR:
            self._cache_grade = "fail"
        elif max_dev > GC_SEQUENCE_WARN:
            self._cache_grade = "warn"
        else:
            self._cache_grade = "pass"
        self._cache_ready = True

    fn grades(self) raises -> List[GradeEntry]:
        var out = List[GradeEntry]()
        out.append(GradeEntry("Per Sequence GC Content", self._cache_grade))
        return out^

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

    fn to_plot(self) raises -> List[PythonObject]:
        """Build GC figure from cached data (summarizer is also the plotter)."""
        var plt = Python.import_module("matplotlib.pyplot")
        var arr = tensor_to_numpy_1d(self._cache_cg_content)
        var x = plt.subplots()
        var fig = x[0]
        var ax = x[1]
        ax.plot(arr, label="GC count per read")
        ax.plot(self._cache_theoretical, label="Theoritical Distribution")
        ax.set_title("GC distribution over all sequences")
        ax.set_xlabel("Mean GC content (%)")
        var figs = List[PythonObject]()
        figs.append(fig)
        return figs^


# ----- Assembled module: Collector + Summarizer + DefaultOutputter -----

struct CGModule(Analyser, StatSummarizer, TextOutput, PlotOutput, HtmlOutput, Copyable, Movable):
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
        self.summarizer.prepare(ctx)

    fn grades(self) raises -> List[GradeEntry]:
        return self.summarizer.grades()

    fn to_data_text(self, ctx: SummaryContext) raises -> String:
        var out = DefaultOutputter()
        return out.wrap_data_block(
            self.summarizer.module_legend(),
            self.summarizer.grades()[0].grade,
            self.summarizer.data_block_body(),
        )

    fn to_plot(self) raises -> List[PythonObject]:
        return self.summarizer.to_plot()

    fn to_html_panels(self) raises -> List[result_panel]:
        var figs = self.summarizer.to_plot()
        var out = DefaultOutputter()
        var panel = out.make_panel(
            self.summarizer.panel_id(),
            self.summarizer.grades()[0].grade,
            self.summarizer.module_legend(),
            figs[0],
        )
        var panels = List[result_panel]()
        panels.append(panel^)
        return panels^


comptime CGContent = CGModule
