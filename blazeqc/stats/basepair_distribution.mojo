"""Base pair distribution: Collector + two Summarizers (N content, Sequence content) + BasepairModule."""

from collections.dict import Dict
from collections.list import List
from python import Python, PythonObject
from blazeseq import FastqRecord, RefRecord
from blazeqc.stats.traits import Collector, Summarizer, PlotOutput
from blazeqc.stats.reporting_traits import FastqcDataOutput, FastqcHtmlOutput, ModuleReport
from blazeqc.stats.summary_utils import SummaryContext, GradeEntry, DefaultOutputter, DataEntry, PanelEntry
from blazeqc.helpers import (
    Matrix2D,
    matrix_to_numpy,
    make_linear_base_groups,
    bin_array,
)
from blazeqc.html_maker import result_panel
from blazeqc.limits import N_CONTENT_WARN, N_CONTENT_ERROR, SEQUENCE_WARN, SEQUENCE_ERROR

comptime BP_WIDTH = 5

# ----- Shared prepared data (filled once, used by both summarizers) -----

struct BasepairPreparedData(Copyable, Movable):
    var bins: List[Int]
    var binned_counts: Matrix2D[DType.int64]
    var binned_pct: Matrix2D[DType.float64]
    var max_length: Int

    fn __init__(out self):
        self.bins = List[Int]()
        self.binned_counts = Matrix2D[DType.int64](0, 0)
        self.binned_pct = Matrix2D[DType.float64](0, 0)
        self.max_length = 0

    fn __copyinit__(out self, other: Self):
        self.bins = other.bins.copy()
        self.binned_counts = other.binned_counts
        self.binned_pct = other.binned_pct
        self.max_length = other.max_length

fn _compute_binned_basepair_data(
    bp_dist: Matrix2D[DType.int64], max_length: Int
) raises -> BasepairPreparedData:
    """Compute binned counts and percentages from raw bp_dist. Shared by both summarizers."""
    var bins = make_linear_base_groups(max_length)
    var n_bins = len(bins)
    var result = BasepairPreparedData()
    if n_bins == 0:
        result.max_length = max_length
        return result^
    var binned_counts = Matrix2D[DType.int64](n_bins, BP_WIDTH)
    for g in range(n_bins):
        var g_start = bins[g] - 1
        var g_end: Int
        if g + 1 < n_bins:
            g_end = bins[g + 1] - 1
        else:
            g_end = max_length
        for row in range(g_start, g_end):
            for col in range(BP_WIDTH):
                binned_counts.add(g, col, bp_dist.get(row, col))
    var binned_pct = Matrix2D[DType.float64](n_bins, BP_WIDTH)
    for g in range(n_bins):
        var row_sum = binned_counts.row_sum(g)
        if row_sum > 0:
            var t = Float64(row_sum)
            for col in range(BP_WIDTH):
                binned_pct.set(
                    g, col,
                    100.0 * Float64(binned_counts.get(g, col)) / t,
                )
    result.bins = bins^
    result.binned_counts = binned_counts^
    result.binned_pct = binned_pct^
    result.max_length = max_length
    return result^


# ----- Collector: tally only -----

struct BasepairCollector(Collector, Copyable, Movable):
    """Collector for per-base base composition (C,G,T,A,N). No summarization or output."""
    var bp_dist: Matrix2D[DType.int64]
    var max_length: Int
    var min_length: Int

    fn __init__(out self):
        self.bp_dist = Matrix2D[DType.int64](1, BP_WIDTH)
        self.max_length = 0
        self.min_length = Int.MAX

    @always_inline
    fn tally_read(mut self, record: FastqRecord):
        var rec_len = len(record)
        if rec_len > self.max_length:
            self.max_length = rec_len
            self.bp_dist.resize(self.max_length, BP_WIDTH)
        if rec_len < self.min_length:
            self.min_length = rec_len
        var seq_span = record.sequence()
        for i in range(rec_len):
            var base_val = Int((seq_span[i] & 0b11111) % BP_WIDTH)
            self.bp_dist.add(i, base_val, 1)

    @always_inline
    fn tally_read(mut self, record: RefRecord):
        var rec_len = len(record)
        if rec_len > self.max_length:
            self.max_length = rec_len
            self.bp_dist.resize(self.max_length, BP_WIDTH)
        if rec_len < self.min_length:
            self.min_length = rec_len
        var seq_span = record.sequence().as_bytes()
        for i in range(rec_len):
            var base_val = Int((seq_span[i] & 0b11111) % BP_WIDTH)
            self.bp_dist.add(i, base_val, 1)


# ----- Plot helpers (used by both summarizers) -----

fn _plot_n_content(
    arr2: PythonObject, py_bins: PythonObject, bins: List[Int]
) raises -> PythonObject:
    var plt = Python.import_module("matplotlib.pyplot")
    var mtp = Python.import_module("matplotlib")
    var y = plt.subplots()
    var fig2 = y[0]
    var ax2 = y[1]
    ax2.plot(arr2)
    var bins_range_2 = Python.list()
    for i in range(len(bins)):
        bins_range_2.append(i)
    ax2.set_xticks(bins_range_2)
    ax2.set_xticklabels(py_bins, rotation=45)
    ax2.xaxis.set_major_locator(
        mtp.ticker.MaxNLocator(integer=True, nbins=15)
    )
    ax2.set_ylim(0, 100)
    var legend_labels_n = Python.list()
    legend_labels_n.append("%N")
    ax2.legend(legend_labels_n)
    ax2.set_xlabel("Position in read (bp)")
    ax2.set_title("N content across all bases")
    return fig2


fn _plot_sequence_content(
    arr1: PythonObject, py_bins: PythonObject, bins: List[Int]
) raises -> PythonObject:
    var plt = Python.import_module("matplotlib.pyplot")
    var mtp = Python.import_module("matplotlib")
    var x = plt.subplots()
    var fig = x[0]
    var ax = x[1]
    var bins_range = Python.list()
    for i in range(len(bins)):
        bins_range.append(i)
    ax.set_xticks(bins_range)
    ax.set_xticklabels(py_bins, rotation=45)
    ax.plot(arr1)
    ax.set_ylim(0, 100)
    ax.xaxis.set_major_locator(
        mtp.ticker.MaxNLocator(integer=True, nbins=15)
    )
    var legend_labels = Python.list()
    legend_labels.append("%C")
    legend_labels.append("%G")
    legend_labels.append("%T")
    legend_labels.append("%A")
    ax.legend(legend_labels)
    ax.set_xlabel("Position in read (bp)")
    ax.set_title("Sequence content across all bases")
    return fig


fn _status_n_content_from_prepared(prepared: BasepairPreparedData) -> String:
    var max_n_pct: Float64 = 0.0
    for g in range(len(prepared.bins)):
        var row_sum = prepared.binned_counts.row_sum(g)
        if row_sum == 0:
            continue
        var n_pct = (Float64(prepared.binned_counts.get(g, 4)) / Float64(row_sum)) * 100.0
        if n_pct > max_n_pct:
            max_n_pct = n_pct
    if max_n_pct > N_CONTENT_ERROR:
        return "fail"
    if max_n_pct > N_CONTENT_WARN:
        return "warn"
    return "pass"


fn _status_sequence_content_from_prepared(prepared: BasepairPreparedData) -> String:
    var max_dev: Float64 = 0.0
    for g in range(len(prepared.bins)):
        var row_sum = prepared.binned_counts.row_sum(g)
        if row_sum == 0:
            continue
        var t = Float64(row_sum)
        var pct_a = (Float64(prepared.binned_counts.get(g, 3)) / t) * 100.0
        var pct_t = (Float64(prepared.binned_counts.get(g, 2)) / t) * 100.0
        var pct_c = (Float64(prepared.binned_counts.get(g, 0)) / t) * 100.0
        var pct_g = (Float64(prepared.binned_counts.get(g, 1)) / t) * 100.0
        var dev_at = pct_a - pct_t
        if dev_at < 0:
            dev_at = -dev_at
        var dev_cg = pct_c - pct_g
        if dev_cg < 0:
            dev_cg = -dev_cg
        var dev = dev_at if dev_at > dev_cg else dev_cg
        if dev > max_dev:
            max_dev = dev
    if max_dev > SEQUENCE_ERROR:
        return "fail"
    if max_dev > SEQUENCE_WARN:
        return "warn"
    return "pass"


# ----- Summarizer: Per Base N Content -----

struct PerBaseNContentSummarizer(Summarizer, PlotOutput, Copyable, Movable):
    var _prepared: BasepairPreparedData
    var _cache_grade: String
    var _cache_ready: Bool

    fn __init__(out self):
        self._prepared = BasepairPreparedData()
        self._cache_grade = ""
        self._cache_ready = False

    fn feed_prepared(mut self, prepared: BasepairPreparedData):
        self._prepared = prepared.copy()

    fn summerize(mut self, ctx: SummaryContext) raises:
        self._cache_grade = _status_n_content_from_prepared(self._prepared)
        self._cache_ready = True

    fn grade(self) raises -> GradeEntry:
        return GradeEntry("Per Base N Content", self._cache_grade)

    fn data_block_body(self) raises -> String:
        if not self._cache_ready:
            return ""
        var out = "#Base\tN-Count\n"
        var n_bins = len(self._prepared.bins)
        for g in range(n_bins):
            var start_bp = self._prepared.bins[g]
            var end_bp: Int
            if g + 1 < n_bins:
                end_bp = self._prepared.bins[g + 1] - 1
            else:
                end_bp = self._prepared.max_length
            var base_label: String
            if end_bp - start_bp + 1 == 1:
                base_label = String(start_bp)
            else:
                base_label = "{}-{}".format(start_bp, end_bp)
            var n_count = self._prepared.binned_counts.get(g, 4)
            out += "{}\t{}\n".format(base_label, n_count)
        return out

    fn module_legend(self) -> String:
        return "Per Base N Content"

    fn panel_id(self) -> String:
        return "n_percentage"

    fn plot_result(self) raises -> PythonObject:
        var np = Python.import_module("numpy")
        var arr = matrix_to_numpy(self._prepared.binned_pct)
        var slice_all = Python.evaluate("slice(None)")
        var arr2 = arr.__getitem__(Python.tuple(slice_all, Python.evaluate("slice(4, 5)")))
        var py_bins = Python.list()
        for i in range(len(self._prepared.bins)):
            py_bins.append(self._prepared.bins[i])
        return _plot_n_content(arr2, py_bins, self._prepared.bins)


# ----- Summarizer: Per Base Sequence Content -----

struct PerBaseSequenceContentSummarizer(Summarizer, PlotOutput, Copyable, Movable):
    var _prepared: BasepairPreparedData
    var _cache_grade: String
    var _cache_ready: Bool

    fn __init__(out self):
        self._prepared = BasepairPreparedData()
        self._cache_grade = ""
        self._cache_ready = False

    fn feed_prepared(mut self, prepared: BasepairPreparedData):
        self._prepared = prepared.copy()

    fn summerize(mut self, ctx: SummaryContext) raises:
        self._cache_grade = _status_sequence_content_from_prepared(self._prepared)
        self._cache_ready = True

    fn grade(self) raises -> GradeEntry:
        return GradeEntry("Per Base Sequence Content", self._cache_grade)

    fn data_block_body(self) raises -> String:
        if not self._cache_ready:
            return ""
        var out = "#Base\tG\tA\tT\tC\n"
        var n_bins = len(self._prepared.bins)
        for g in range(n_bins):
            var start_bp = self._prepared.bins[g]
            var end_bp: Int
            if g + 1 < n_bins:
                end_bp = self._prepared.bins[g + 1] - 1
            else:
                end_bp = self._prepared.max_length
            var base_label: String
            if end_bp - start_bp + 1 == 1:
                base_label = String(start_bp)
            else:
                base_label = "{}-{}".format(start_bp, end_bp)
            var gu = self._prepared.binned_pct.get(g, 1)
            var au = self._prepared.binned_pct.get(g, 3)
            var tu = self._prepared.binned_pct.get(g, 2)
            var cu = self._prepared.binned_pct.get(g, 0)
            out += "{}\t{}\t{}\t{}\t{}\n".format(base_label, gu, au, tu, cu)
        return out

    fn module_legend(self) -> String:
        return "Per Base Sequence Content"

    fn panel_id(self) -> String:
        return "base_pair_distribution"

    fn plot_result(self) raises -> PythonObject:
        var np = Python.import_module("numpy")
        var arr = matrix_to_numpy(self._prepared.binned_pct)
        var slice_all = Python.evaluate("slice(None)")
        var arr1 = arr.__getitem__(Python.tuple(slice_all, Python.evaluate("slice(0, 4)")))
        var py_bins = Python.list()
        for i in range(len(self._prepared.bins)):
            py_bins.append(self._prepared.bins[i])
        return _plot_sequence_content(arr1, py_bins, self._prepared.bins)


# ----- Assembled module -----

struct BasepairModule(FastqcDataOutput, FastqcHtmlOutput, Copyable, Movable):
    """Module assembling BasepairCollector + PerBaseNContentSummarizer + PerBaseSequenceContentSummarizer."""
    var collector: BasepairCollector
    var summarizer_n: PerBaseNContentSummarizer
    var summarizer_seq: PerBaseSequenceContentSummarizer

    fn __init__(out self):
        self.collector = BasepairCollector()
        self.summarizer_n = PerBaseNContentSummarizer()
        self.summarizer_seq = PerBaseSequenceContentSummarizer()

    fn prepare_summarizers(mut self, ctx: SummaryContext) raises:
        """Compute binned data and feed both summarizers; then summerize both."""
        var prepared_n = _compute_binned_basepair_data(
            self.collector.bp_dist, self.collector.max_length
        )
        self.summarizer_n.feed_prepared(prepared_n^)
        self.summarizer_n.summerize(ctx)
        var prepared_seq = _compute_binned_basepair_data(
            self.collector.bp_dist, self.collector.max_length
        )
        self.summarizer_seq.feed_prepared(prepared_seq^)
        self.summarizer_seq.summerize(ctx)

    fn to_html(self) raises -> result_panel:
        var out = DefaultOutputter()
        var fig = self.summarizer_seq.plot_result()
        return out.make_panel(
            self.summarizer_seq.panel_id(),
            self.summarizer_seq.grade().grade,
            self.summarizer_seq.module_legend(),
            fig,
        )

    fn data_entries(self, ctx: SummaryContext) raises -> List[DataEntry]:
        var body_seq = self.summarizer_seq.data_block_body()
        var g_seq = self.summarizer_seq.grade()
        var body_n = self.summarizer_n.data_block_body()
        var g_n = self.summarizer_n.grade()
        var entries = List[DataEntry]()
        entries.append(DataEntry(self.summarizer_seq.module_legend(), g_seq.grade, body_seq))
        entries.append(DataEntry(self.summarizer_n.module_legend(), g_n.grade, body_n))
        return entries^

    fn panel_entries(self, figures: List[PythonObject]) raises -> List[PanelEntry]:
        var entries = List[PanelEntry]()
        entries.append(PanelEntry(
            self.summarizer_seq.panel_id(),
            self.summarizer_seq.grade().grade,
            self.summarizer_seq.module_legend(),
            "image",
            figures[1],
            "",
        ))
        entries.append(PanelEntry(
            self.summarizer_n.panel_id(),
            self.summarizer_n.grade().grade,
            self.summarizer_n.module_legend(),
            "image",
            figures[0],
            "",
        ))
        return entries^

    fn plot_result(self) raises -> Tuple[PythonObject, PythonObject]:
        """Return (N figure, sequence content figure) for legacy plot() API."""
        return (
            self.summarizer_n.plot_result(),
            self.summarizer_seq.plot_result(),
        )
