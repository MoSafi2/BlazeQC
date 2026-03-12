"""Quality distribution: one collector, two summarizers, assembled QualityModule."""

from python import Python, PythonObject
from collections.list import List
from blazeseq import FastqRecord, RefRecord
from blazeqc.stats.traits import Collector, Summarizer, PlotOutput
from blazeqc.stats.summary_utils import SummaryContext, GradeEntry, DefaultOutputter
from blazeqc.helpers import (
    Matrix2D,
    matrix_to_numpy,
    tensor_to_numpy_1d,
    make_linear_base_groups,
    bin_array,
)
from blazeqc.CONSTS import QualitySchema
from blazeqc.html_maker import result_panel
from blazeqc.limits import (
    QUALITY_BASE_LOWER_WARN,
    QUALITY_BASE_LOWER_ERROR,
    QUALITY_BASE_MEDIAN_WARN,
    QUALITY_BASE_MEDIAN_ERROR,
    QUALITY_SEQUENCE_WARN,
    QUALITY_SEQUENCE_ERROR,
)


# ----- Collector: tally only -----

struct QualityCollector(Collector, Copyable, Movable):
    """Collector for per-base and per-sequence quality (no summarization or output)."""
    var qu_dist: Matrix2D[DType.int64]
    var qu_dist_seq: List[Int64]
    var max_length: Int
    var max_qu: UInt8
    var min_qu: UInt8

    fn __init__(out self) raises:
        self.qu_dist = Matrix2D[DType.int64](1, 128)
        self.qu_dist_seq = List[Int64](capacity=128)
        for _ in range(128):
            self.qu_dist_seq.append(0)
        self.max_length = 0
        self.max_qu = 0
        self.min_qu = 128

    fn tally_read(mut self, record: FastqRecord):
        if len(record) > self.max_length:
            self.max_length = len(record)
            self.qu_dist.resize(self.max_length, 128)

        var qu_span = record.quality()
        for i in range(len(record)):
            var base_qu = qu_span[i]
            self.qu_dist.add(i, Int(base_qu), 1)
            if base_qu > self.max_qu:
                self.max_qu = base_qu
            if base_qu < self.min_qu:
                self.min_qu = base_qu

        var qu_sum: Int = 0
        for i in range(len(record)):
            qu_sum += Int(qu_span[i])
        var average = Int(qu_sum / len(record))
        while len(self.qu_dist_seq) <= average:
            self.qu_dist_seq.append(0)
        self.qu_dist_seq[average] += 1

    fn tally_read(mut self, record: RefRecord):
        if len(record) > self.max_length:
            self.max_length = len(record)
            self.qu_dist.resize(self.max_length, 128)

        var qu_span = record.quality().as_bytes()
        for i in range(len(record)):
            var base_qu = qu_span[i]
            self.qu_dist.add(i, Int(base_qu), 1)
            if base_qu > self.max_qu:
                self.max_qu = base_qu
            if base_qu < self.min_qu:
                self.min_qu = base_qu

        var qu_sum: Int = 0
        for i in range(len(record)):
            qu_sum += Int(qu_span[i])
        var average = Int(qu_sum / len(record))
        while len(self.qu_dist_seq) <= average:
            self.qu_dist_seq.append(0)
        self.qu_dist_seq[average] += 1

    fn _guess_schema(self) -> QualitySchema:
        comptime ILLUMINA_1_3_ENCODING_OFFSET = 64
        if self.min_qu < 64:
            return QualitySchema("Illumina v1.8", 33, 126, 33)
        elif self.min_qu == ILLUMINA_1_3_ENCODING_OFFSET + 1:
            return QualitySchema("Illumina v1.3", 64, 126, 64)
        elif self.min_qu <= 126:
            return QualitySchema("Illumina v1.5", 66, 126, 64)
        else:
            print("Unable to parse Quality Schema, returning generic schema")
            return QualitySchema("Generic", 33, 126, 33)


# ----- Helpers (used by summarizers) -----

fn _slice_array(arr: PythonObject, min_index: Int, max_index: Int) raises -> PythonObject:
    var np = Python.import_module("numpy")
    var indices = np.arange(min_index, max_index)
    return np.take(arr, indices, axis=1)


fn _plot_per_base_quality(
    mean_line: PythonObject,
    boxplot_stats: PythonObject,
    bins: List[Int],
    py_bins: PythonObject,
) raises -> PythonObject:
    """Plot quality scores across all bases (boxplot + mean line)."""
    var plt = Python.import_module("matplotlib.pyplot")
    var mtp = Python.import_module("matplotlib")
    var x = plt.subplots()
    var fig = x[0]
    var ax = x[1]
    ax.bxp(boxplot_stats, showfliers=False)
    ax.plot(mean_line)
    ax.set_ylim(0, 60)
    ax.set_title("Quality scores across all bases")
    ax.set_xlabel("Position in read (bp)")
    var bins_range = Python.list()
    for i in range(len(bins)):
        bins_range.append(i)
    ax.set_xticks(bins_range)
    ax.set_xticklabels(py_bins, rotation=45)
    ax.xaxis.set_major_locator(
        mtp.ticker.MaxNLocator(integer=True, nbins=15)
    )
    return fig


fn _plot_per_sequence_quality(arr2: PythonObject) raises -> PythonObject:
    """Plot quality score distribution over all sequences."""
    var plt = Python.import_module("matplotlib.pyplot")
    var z = plt.subplots()
    var fig2 = z[0]
    var ax2 = z[1]
    ax2.plot(arr2)
    ax2.set_xlabel("Mean Sequence Quality (Phred Score)")
    ax2.set_title("Quality score distribution over all sequences")
    return fig2


# ----- Summarizer: Per Base Sequence Quality -----

struct PerBaseQualitySummarizer(Summarizer, PlotOutput, Copyable, Movable):
    """Summarization and plotting for per-base sequence quality. Caches only heavy numpy/plot results."""
    var _cache_qu_dist: Matrix2D[DType.int64]
    var _cache_max_length: Int
    var _cache_mean_line: PythonObject
    var _cache_boxplot_stats: PythonObject
    var _cache_py_bins: PythonObject
    var _cache_bins: List[Int]
    var _cache_p10: PythonObject
    var _cache_p90: PythonObject
    var _cache_nrows: Int
    var _cache_ready: Bool

    fn __init__(out self) raises:
        self._cache_qu_dist = Matrix2D[DType.int64](1, 128)
        self._cache_max_length = 0
        self._cache_mean_line = Python.evaluate("None")
        self._cache_boxplot_stats = Python.evaluate("None")
        self._cache_py_bins = Python.evaluate("None")
        self._cache_bins = List[Int]()
        self._cache_p10 = Python.evaluate("None")
        self._cache_p90 = Python.evaluate("None")
        self._cache_nrows = 0
        self._cache_ready = False

    fn feed(mut self, collector: QualityCollector):
        """Copy collector state needed for per-base summarization."""
        self._cache_qu_dist = collector.qu_dist.copy()
        self._cache_max_length = collector.max_length

    fn summerize(mut self, ctx: SummaryContext) raises:
        """Compute binned per-base stats (mean line, boxplot, percentiles)."""
        var np = Python.import_module("numpy")
        var schema = _guess_schema_from_collector(self._cache_qu_dist, self._cache_max_length)
        var arr = matrix_to_numpy(self._cache_qu_dist)
        var min_index = schema.OFFSET
        var max_qu_int = 0
        for row in range(self._cache_qu_dist.rows):
            for col in range(self._cache_qu_dist.cols):
                if self._cache_qu_dist.get(row, col) != 0 and col > max_qu_int:
                    max_qu_int = col
        var max_index = 40 if 40 > max_qu_int else max_qu_int
        arr = _slice_array(arr, Int(min_index), Int(max_index))
        var bins = make_linear_base_groups(self._cache_max_length)
        var py_bins: PythonObject
        arr, py_bins = bin_array(arr, bins, func="mean")

        var nrows = Int(py=arr.shape[0])
        var ncols = Int(py=arr.shape[1])
        var mean_line = np.sum(
            arr * np.arange(1, ncols + 1), axis=1
        ) / np.sum(arr, axis=1)
        var cum_sum = np.cumsum(arr, axis=1)
        var total_counts = np.reshape(np.sum(arr, axis=1), Python.tuple(nrows, 1))
        var median = np.argmax(cum_sum > total_counts / 2, axis=1)
        var Q75 = np.argmax(cum_sum > total_counts * 0.75, axis=1)
        var Q25 = np.argmax(cum_sum > total_counts * 0.25, axis=1)
        var p10 = np.argmax(cum_sum > total_counts * 0.1, axis=1)
        var p90 = np.argmax(cum_sum > total_counts * 0.9, axis=1)
        var IQR = Q75 - Q25
        var py_none = Python.evaluate("None")
        var whislo = np.full(len(IQR), py_none)
        var whishi = np.full(len(IQR), py_none)
        var boxplot_stats = Python.list()
        for i in range(len(IQR)):
            var stat: PythonObject = Python.dict()
            stat["med"] = median[i]
            stat["q1"] = Q25[i]
            stat["q3"] = Q75[i]
            stat["whislo"] = whislo[i]
            stat["whishi"] = whishi[i]
            boxplot_stats.append(stat)

        self._cache_mean_line = mean_line
        self._cache_boxplot_stats = boxplot_stats
        self._cache_py_bins = py_bins
        self._cache_bins = bins^
        self._cache_p10 = p10
        self._cache_p90 = p90
        self._cache_nrows = nrows
        self._cache_ready = True

    fn _build_data_block_body(self) raises -> String:
        """Build data block body from cached numpy arrays (called from summerize)."""
        var out = String()
        out += "#Base\tMean\tMedian\tLower Quartile\tUpper Quartile\t10th Percentile\t90th Percentile\n"
        for i in range(self._cache_nrows):
            var start_bp = self._cache_bins[i]
            var end_bp: Int = self._cache_max_length
            if i + 1 < len(self._cache_bins):
                end_bp = self._cache_bins[i + 1] - 1
            var base_label: String
            if start_bp == end_bp:
                base_label = String(start_bp)
            else:
                base_label = "{}-{}".format(start_bp, end_bp)
            var mean_val = Float64(py=self._cache_mean_line[i]) - 1.0
            var med_val = Float64(py=self._cache_boxplot_stats[i]["med"])
            var q1_val = Float64(py=self._cache_boxplot_stats[i]["q1"])
            var q3_val = Float64(py=self._cache_boxplot_stats[i]["q3"])
            var p10_val = Float64(py=self._cache_p10[i])
            var p90_val = Float64(py=self._cache_p90[i])
            out += "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(base_label, mean_val, med_val, q1_val, q3_val, p10_val, p90_val)
        return out

    fn _get_status_per_base(self) -> String:
        """Status from per-base quality: lower quartile and median per group. Schema computed on the fly."""
        var schema = _guess_schema_from_collector(self._cache_qu_dist, self._cache_max_length)
        var offset = Int(schema.OFFSET)
        var bins = make_linear_base_groups(self._cache_max_length)
        var min_q25_phred: Float64 = 1e9
        var min_median_phred: Float64 = 1e9
        for b in range(len(bins)):
            var start_row = bins[b] - 1
            var end_row: Int = self._cache_qu_dist.rows
            if b + 1 < len(bins):
                end_row = bins[b + 1] - 1
            if start_row >= end_row:
                continue
            var total: Int64 = 0
            for col in range(self._cache_qu_dist.cols):
                var s: Int64 = 0
                for row in range(start_row, end_row):
                    s += self._cache_qu_dist.get(row, col)
                total += s
            if total == 0:
                continue
            var cum: Int64 = 0
            var q25_phred: Float64 = -1.0
            var median_phred: Float64 = -1.0
            for col in range(self._cache_qu_dist.cols):
                for row in range(start_row, end_row):
                    cum += self._cache_qu_dist.get(row, col)
                if q25_phred < 0 and Float64(cum) >= 0.25 * Float64(total):
                    q25_phred = Float64(col) - Float64(offset)
                if median_phred < 0 and Float64(cum) >= 0.5 * Float64(total):
                    median_phred = Float64(col) - Float64(offset)
                    break
            if q25_phred >= 0 and q25_phred < min_q25_phred:
                min_q25_phred = q25_phred
            if median_phred >= 0 and median_phred < min_median_phred:
                min_median_phred = median_phred
        if min_q25_phred > 1e8:
            min_q25_phred = 0.0
        if min_median_phred > 1e8:
            min_median_phred = 0.0
        if min_q25_phred < QUALITY_BASE_LOWER_ERROR or min_median_phred < QUALITY_BASE_MEDIAN_ERROR:
            return "fail"
        if min_q25_phred < QUALITY_BASE_LOWER_WARN or min_median_phred < QUALITY_BASE_MEDIAN_WARN:
            return "warn"
        return "pass"

    fn grade(self) -> GradeEntry:
        return GradeEntry("Per Base Sequence Quality", self._get_status_per_base())

    fn data_block_body(self) raises -> String:
        """Build data block on the fly from cached arrays (cheap)."""
        if not self._cache_ready:
            return ""
        return self._build_data_block_body()

    fn module_legend(self) -> String:
        return "Per Base Sequence Quality"

    fn panel_id(self) -> String:
        return "qu_score_dis_base"

    fn plot_result(self) raises -> PythonObject:
        return _plot_per_base_quality(
            self._cache_mean_line,
            self._cache_boxplot_stats,
            self._cache_bins,
            self._cache_py_bins,
        )


# ----- Summarizer: Per Sequence Quality Scores -----

struct PerSequenceQualitySummarizer(Summarizer, PlotOutput, Copyable, Movable):
    """Summarization and plotting for per-sequence quality scores. Caches only heavy arr2; data block built on the fly."""
    var _cache_qu_dist_seq: List[Int64]
    var _cache_arr2: PythonObject
    var _cache_schema_offset: Int
    var _cache_ready: Bool

    fn __init__(out self) raises:
        self._cache_qu_dist_seq = List[Int64]()
        self._cache_arr2 = Python.evaluate("None")
        self._cache_schema_offset = 0
        self._cache_ready = False

    fn feed(mut self, collector: QualityCollector):
        """Copy collector per-sequence counts and schema offset."""
        self._cache_qu_dist_seq = List[Int64](capacity=len(collector.qu_dist_seq))
        for i in range(len(collector.qu_dist_seq)):
            self._cache_qu_dist_seq.append(collector.qu_dist_seq[i])
        var schema = collector._guess_schema()
        self._cache_schema_offset = Int(schema.OFFSET)

    fn summerize(mut self, ctx: SummaryContext) raises:
        """Build trimmed 1d array for plot and data block."""
        var schema_offset = self._cache_schema_offset
        var index = 0
        for i in range(len(self._cache_qu_dist_seq) - 1, -1, -1):
            if self._cache_qu_dist_seq[i] != 0:
                index = i
                break
        var arr2 = tensor_to_numpy_1d(self._cache_qu_dist_seq)
        arr2 = arr2[schema_offset : index + 2]
        self._cache_arr2 = arr2
        self._cache_ready = True

    fn _build_data_block_body(self) raises -> String:
        """Build data block body from cached arr2 (called from summerize)."""
        var out = String()
        out += "#Quality\tCount\n"
        var arr2 = self._cache_arr2
        var n_seq = Int(py=arr2.shape[0])
        for i in range(n_seq):
            var q_val = Float64(i) + Float64(self._cache_schema_offset)
            var c_val = Float64(py=arr2[i])
            out += "{}\t{}\n".format(q_val, c_val)
        return out

    fn _get_status_per_sequence(self) -> String:
        """Status from per-sequence quality: most frequent mean Phred."""
        var offset = self._cache_schema_offset
        var max_count: Int64 = 0
        var mode_phred: Float64 = 0.0
        for i in range(len(self._cache_qu_dist_seq)):
            var c = self._cache_qu_dist_seq[i]
            if c > max_count:
                max_count = c
                mode_phred = Float64(i) - Float64(offset)
        if max_count == 0:
            return "pass"
        if mode_phred < QUALITY_SEQUENCE_ERROR:
            return "fail"
        if mode_phred < QUALITY_SEQUENCE_WARN:
            return "warn"
        return "pass"

    fn grade(self) -> GradeEntry:
        return GradeEntry("Per Sequence Quality Scores", self._get_status_per_sequence())

    fn data_block_body(self) raises -> String:
        """Build data block on the fly from cached arr2 (cheap)."""
        if not self._cache_ready:
            return ""
        return self._build_data_block_body()

    fn module_legend(self) -> String:
        return "Per Sequence Quality Scores"

    fn panel_id(self) -> String:
        return "qu_score_dis_seq"

    fn plot_result(self) raises -> PythonObject:
        return _plot_per_sequence_quality(self._cache_arr2)


# ----- Helper: guess schema from cached qu_dist (no collector reference) -----
# PerBaseQualitySummarizer has no reference to collector after feed; we infer from data.
fn _guess_schema_from_collector(qu_dist: Matrix2D[DType.int64], max_length: Int) -> QualitySchema:
    comptime ILLUMINA_1_3_ENCODING_OFFSET = 64
    var min_qu: UInt8 = 128
    for row in range(qu_dist.rows):
        for col in range(qu_dist.cols):
            if qu_dist.get(row, col) != 0 and UInt8(col) < min_qu:
                min_qu = UInt8(col)
    if min_qu < 64:
        return QualitySchema("Illumina v1.8", 33, 126, 33)
    elif min_qu == ILLUMINA_1_3_ENCODING_OFFSET + 1:
        return QualitySchema("Illumina v1.3", 64, 126, 64)
    elif min_qu <= 126:
        return QualitySchema("Illumina v1.5", 66, 126, 64)
    else:
        return QualitySchema("Generic", 33, 126, 33)


# ----- Assembled module: Collector + two Summarizers -----

struct QualityModule(Copyable, Movable):
    """Module assembling QualityCollector + PerBaseQualitySummarizer + PerSequenceQualitySummarizer."""
    var collector: QualityCollector
    var summarizer_base: PerBaseQualitySummarizer
    var summarizer_seq: PerSequenceQualitySummarizer

    fn __init__(out self) raises:
        self.collector = QualityCollector()
        self.summarizer_base = PerBaseQualitySummarizer()
        self.summarizer_seq = PerSequenceQualitySummarizer()

    fn to_data_text(self, ctx: SummaryContext) raises -> String:
        var out = DefaultOutputter()
        var body_base = self.summarizer_base.data_block_body()
        var g_base = self.summarizer_base.grade()
        var body_seq = self.summarizer_seq.data_block_body()
        var g_seq = self.summarizer_seq.grade()
        return out.wrap_data_block(
            self.summarizer_base.module_legend(), g_base.grade, body_base
        ) + out.wrap_data_block(
            self.summarizer_seq.module_legend(), g_seq.grade, body_seq
        )

    fn to_html_panels(self) raises -> List[result_panel]:
        var panels = List[result_panel]()
        var out = DefaultOutputter()
        var fig_base = self.summarizer_base.plot_result()
        var panel_base = out.make_panel(
            self.summarizer_base.panel_id(),
            self.summarizer_base.grade().grade,
            self.summarizer_base.module_legend(),
            fig_base,
        )
        panels.append(panel_base^)
        var fig_seq = self.summarizer_seq.plot_result()
        var panel_seq = out.make_panel(
            self.summarizer_seq.panel_id(),
            self.summarizer_seq.grade().grade,
            self.summarizer_seq.module_legend(),
            fig_seq,
        )
        panels.append(panel_seq^)
        return panels^

    fn plot(self) raises -> Tuple[PythonObject, PythonObject]:
        return Tuple(
            self.summarizer_base.plot_result(),
            self.summarizer_seq.plot_result(),
        )

    fn _guess_schema(self) -> QualitySchema:
        return self.collector._guess_schema()
