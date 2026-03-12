"""Length distribution (split from stats_.mojo)."""

from python import Python, PythonObject
from collections.dict import Dict
from collections.list import List
from blazeseq import FastqRecord, RefRecord
from blazeqc.stats.traits import (
    Collector,
    Summarizer,
    PlotOutput,
    FastqcHtmlOutput,
    FastqcDataOutput,
    ModuleReport,
)
from blazeqc.stats.summary_utils import SummaryContext, GradeEntry, DefaultOutputter
from blazeqc.helpers import tensor_to_numpy_1d, bin_array
from blazeqc.html_maker import result_panel


# ----- Collector: tally only -----

struct LengthCollector(Collector, Copyable, Movable):
    """Collector for sequence length distribution (no summarization or output)."""
    var length_vector: List[Int64]
    var zero_length_count: Int

    fn __init__(out self) raises:
        self.length_vector = List[Int64]()
        self.zero_length_count = 0

    @always_inline
    fn tally_read(mut self, record: FastqRecord):
        if len(record) == 0:
            self.zero_length_count += 1
            return
        while len(self.length_vector) < len(record):
            self.length_vector.append(0)
        self.length_vector[len(record) - 1] += 1

    @always_inline
    fn tally_read(mut self, record: RefRecord):
        if len(record) == 0:
            self.zero_length_count += 1
            return
        while len(self.length_vector) < len(record):
            self.length_vector.append(0)
        self.length_vector[len(record) - 1] += 1

    @always_inline
    fn length_average(self, num_reads: Int) -> Float64:
        var cum: Int64 = 0
        for i in range(len(self.length_vector)):
            cum += self.length_vector[i] * (i + 1)
        return Int(cum) / num_reads



# ----- Summarizer: prepare + grades + data for output + to_plot -----

struct LengthSummarizer(Summarizer, PlotOutput, Copyable, Movable):
    """Summarization and plotting for sequence length distribution."""
    var _cache_length_vector: List[Int64]
    var _cache_zero_length_count: Int

    var _cache_binned_arr: PythonObject
    var _cache_ticks: PythonObject
    var _cache_labels: PythonObject
    var _cache_xlim_left: Int
    var _cache_xlim_right: Int
    var _cache_status: String
    var _cache_ready: Bool

    fn __init__(out self) raises:
        self._cache_length_vector = List[Int64]()
        self._cache_zero_length_count = 0
        self._cache_binned_arr = Python.evaluate("None")
        self._cache_ticks = Python.evaluate("None")
        self._cache_labels = Python.evaluate("None")
        self._cache_xlim_left = 0
        self._cache_xlim_right = 0
        self._cache_status = ""
        self._cache_ready = False

    fn feed_(mut self, collector: LengthCollector):
        """Copy collector counts into cache for summerize(ctx)."""
        self._cache_length_vector = List[Int64](capacity=len(collector.length_vector))
        for i in range(len(collector.length_vector)):
            self._cache_length_vector.append(collector.length_vector[i])
        self._cache_zero_length_count = collector.zero_length_count

    fn summerize(mut self, ctx: SummaryContext) raises:
        """Compute binned length distribution and status from cached counts."""
        var min_len, max_len = _min_max_lengths_from_vector(self._cache_length_vector)
        var starting, interval = self.get_size_distribution_(min_len, max_len)

        var bins = _build_bins(starting, interval, max_len)
        var x_categories = _build_x_categories(bins, interval, max_len)
        var binned_arr = _bin_length_counts(self._cache_length_vector, bins)
        var ticks, labels = _build_ticks_and_labels(x_categories)

        var num_bins = len(bins)
        var first_nonzero_index = _first_nonzero_bin_index(binned_arr, num_bins)
        var xlim_left, xlim_right = _xlim_from_bins(first_nonzero_index, num_bins)

        self._cache_binned_arr = binned_arr
        self._cache_ticks = ticks
        self._cache_labels = labels
        self._cache_xlim_left = xlim_left
        self._cache_xlim_right = xlim_right
        self._cache_status = self._get_status()
        self._cache_ready = True

    fn get_size_distribution_(
        self, min_val: Int, max_val: Int
    ) -> Tuple[Int, Int]:
        # We won't group if they've asked us not to
        var base = 1
        while base > (max_val - min_val):
            base //= 10
        var divisions: List[Int] = [1, 2, 5]

        while True:
            for d in divisions:
                var tester = base * d
                if (max_val - min_val) / tester <= 50:
                    interval = tester
                    break
            else:
                base *= 10
                continue
            break

        # Now we work out the first value to be plotted
        var basic_division = min_val // interval
        var test_start = basic_division * interval
        var starting = test_start

        return starting, interval

    fn grade(self) -> GradeEntry:
        return GradeEntry("Sequence Length Distribution", self._cache_status)


    fn data_block_body(self) -> String:
        """Module-specific lines (header + data)."""
        var out = "#Length\tCount\n"
        for i in range(len(self._cache_length_vector)):
            if self._cache_length_vector[i] > 0:
                out += "{}\t{}\n".format(i + 1, self._cache_length_vector[i])
        return out

    fn module_legend(self) -> String:
        return "Sequence Length Distribution"

    fn panel_id(self) -> String:
        return "seq_len_dis"

    fn plot_result(self) raises -> PythonObject:
        """Build figure from cached binned distribution (expects summerize(ctx) to be called first)."""
        var plt = Python.import_module("matplotlib.pyplot")
        var np = Python.import_module("numpy")
        var mtp = Python.import_module("matplotlib")

        var x = plt.subplots()
        var fig = x[0]
        var ax = x[1]
        ax.plot(self._cache_binned_arr)
        ax.set_xticks(self._cache_ticks)
        ax.set_xticklabels(self._cache_labels, rotation=45)
        ax.xaxis.set_major_locator(
            mtp.ticker.MaxNLocator(integer=True, nbins=15)
        )
        ax.set_xlim(self._cache_xlim_left, self._cache_xlim_right)
        ax.set_ylim(0)
        ax.set_title("Distribution of sequence lengths over all sequences")
        ax.set_xlabel("Sequence Length (bp)")
        ax.set_ylabel("Number of Reads")
        return fig

    fn _get_status(self) -> String:
        # Error if any zero-length sequences exist (matching FastQC raisesError)
        if self._cache_zero_length_count > 0:
            return "fail"
        # Warning if sequences have more than one distinct length (matching FastQC raisesWarning)
        var distinct_lengths: Int = 0
        for i in range(len(self._cache_length_vector)):
            if self._cache_length_vector[i] > 0:
                distinct_lengths += 1
        if distinct_lengths > 1:
            return "warn"
        return "pass"


# ----- Assembled module: Collector + Summarizer + DefaultOutputter -----

struct LengthModule(FastqcDataOutput, FastqcHtmlOutput, ModuleReport, Copyable, Movable):
    """Module assembling LengthCollector + LengthSummarizer; exposes only data/text/HTML helpers."""
    var collector: LengthCollector
    var summarizer: LengthSummarizer

    fn __init__(out self) raises:
        self.collector = LengthCollector()
        self.summarizer = LengthSummarizer()

    fn to_data_text(self, ctx: SummaryContext) raises -> String:
        """FastQC-style data block text for this module."""
        var body = self.summarizer.data_block_body()
        var g = self.summarizer.grade()
        var out = DefaultOutputter()
        return out.wrap_data_block(self.summarizer.module_legend(), g.grade, body)

    fn to_html(self) raises -> result_panel:
        var fig = self.summarizer.plot_result()
        var out = DefaultOutputter()
        return out.make_panel(
            self.summarizer.panel_id(),
            self.summarizer.grade().grade,
            self.summarizer.module_legend(),
            fig,
        )

    fn to_html_panels(self, figures: List[PythonObject]) raises -> Dict[String, result_panel]:
        var out = DefaultOutputter()
        var panel = out.make_panel(
            self.summarizer.panel_id(),
            self.summarizer.grade().grade,
            self.summarizer.module_legend(),
            figures[0],
        )
        var d = Dict[String, result_panel]()
        d[panel.legand] = panel^
        return d^

    fn plot_result(self) raises -> PythonObject:
        """Delegate to summarizer for plot; ModuleReport entry point."""
        return self.summarizer.plot_result()



# ----- Helper functions for length summarization -----

fn _min_max_lengths_from_vector(length_vector: List[Int64]) -> Tuple[Int, Int]:
    """Return (min_len, max_len) for the cached length vector (1-based min, exclusive max)."""
    var min_len: Int = 0
    var max_len: Int = len(length_vector)
    for i in range(len(length_vector)):
        if length_vector[i] > 0:
            min_len = i + 1
            break
    if min_len > 0:
        min_len -= 1
    max_len += 1
    return min_len, max_len

fn _build_bins(starting: Int, interval: Int, max_len: Int) -> List[Int]:
    """Build list of bin start values from starting up to max_len."""
    var bins = List[Int]()
    var start = starting
    while start <= max_len:
        bins.append(start)
        start += interval
    return bins^

fn _build_x_categories(
    bins: List[Int], interval: Int, max_len: Int
) -> List[String]:
    """Build display labels for each bin (e.g. '1' or '1-10')."""
    var x_categories = List[String]()
    for i in range(len(bins)):
        var min_value = bins[i]
        var max_value = bins[i] + interval - 1
        if max_value > max_len:
            max_value = max_len
        if interval == 1:
            x_categories.append(String(min_value))
        else:
            x_categories.append(String(min_value) + "-" + String(max_value))
    return x_categories^

fn _bin_length_counts(
    length_vector: List[Int64], bins: List[Int]
) raises -> PythonObject:
    """Pad length counts with zeros and bin them; returns binned numpy array."""
    var padded = List[Int64](capacity=len(length_vector) + 2)
    padded.append(0)
    for i in range(len(length_vector)):
        padded.append(length_vector[i])
    padded.append(0)
    var arr_np = tensor_to_numpy_1d(padded)
    var binned_arr, _ = bin_array(arr_np, bins, func="sum")
    return binned_arr

fn _build_ticks_and_labels(
    x_categories: List[String]
) raises -> Tuple[PythonObject, PythonObject]:
    """Build Python lists of tick indices and category labels for the plot."""
    var ticks = Python.list()
    var labels = Python.list()
    for i in range(len(x_categories)):
        ticks.append(i)
        labels.append(x_categories[i])
    return ticks, labels

fn _first_nonzero_bin_index(
    binned_arr: PythonObject, num_bins: Int
) raises -> Int:
    """Index of the first bin with count > 0, or -1 if none."""
    for i in range(num_bins):
        var bin_value = Int64(py=binned_arr.item(i))
        if bin_value > 0:
            return i
    return -1

fn _xlim_from_bins(
    first_nonzero_index: Int, num_bins: Int
) -> Tuple[Int, Int]:
    """Return (xlim_left, xlim_right) for the plot from first nonzero bin and num_bins."""
    var xlim_left = first_nonzero_index - 1
    var xlim_right = num_bins
    return xlim_left, xlim_right

