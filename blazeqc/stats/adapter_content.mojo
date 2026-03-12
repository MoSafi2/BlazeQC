"""Adapter content (split from stats_.mojo). Collector / Summarizer separation."""

from utils import Index
from python import Python, PythonObject
from collections.list import List
from blazeqc.stats.traits import Collector, Summarizer, PlotOutput
from blazeqc.stats.summary_utils import SummaryContext, GradeEntry, DefaultOutputter
from blazeqc.helpers import (
    Matrix2D,
    matrix_to_numpy,
    encode_img_b64,
)
from blazeqc.html_maker import result_panel
from blazeqc.limits import ADAPTER_WARN, ADAPTER_ERROR
from blazeseq import FastqRecord, RefRecord


# TODO: Also plot the Over-represented sequences.
# TODO: Add binning
@fieldwise_init
struct AdapterContentCollector[bits: Int = 3](Collector, Copyable, Movable):
    """Collector for adapter content (kmer hashes along positions). No summarization or output."""
    var kmer_len: Int
    var hash_counts: Matrix2D[DType.int64]
    var hash_list: List[UInt64]
    var max_length: Int

    fn __init__(out self, var hashes: List[UInt64], kmer_len: Int = 0):
        self.kmer_len = min(kmer_len, 64 // Self.bits)
        self.hash_list = hashes^
        self.hash_counts = Matrix2D[DType.int64](len(self.hash_list), 1)
        self.max_length = 0

    @always_inline
    fn tally_read(mut self, record: FastqRecord):
        self.tally_read(record, 0)

    @always_inline
    fn tally_read(mut self, record: RefRecord):
        self.tally_read(record, 0)

    fn tally_read(mut self, record: RefRecord, read_no: Int64):
        var hash: UInt64 = 0
        var end = 0
        var mask: UInt64 = (0b1 << self.kmer_len * Self.bits) - 1
        var neg_mask = mask >> Self.bits
        var bit_shift = (0b1 << Self.bits) - 1

        var rec_len = len(record)
        if rec_len > self.max_length:
            self.max_length = rec_len
            self.hash_counts.resize(len(self.hash_list), self.max_length)

        if len(self.hash_list) > 0:
            self._check_hashes(hash, 1)

        var seq_span = record.sequence().as_bytes()
        for i in range(end, rec_len):
            hash = hash & neg_mask
            var rem = seq_span[i] & bit_shift
            hash = (hash << Self.bits) + Int(rem)
            if len(self.hash_list) > 0:
                self._check_hashes(hash, i + 1)

    # TODO: Check if it will be easier to use the bool_tuple and hashes as a list instead
    @always_inline
    fn tally_read(mut self, record: FastqRecord, read_no: Int64):
        var hash: UInt64 = 0
        var end = 0
        # Make a custom bit mask of 1s by certain length
        var mask: UInt64 = (0b1 << self.kmer_len * Self.bits) - 1
        var neg_mask = mask >> Self.bits
        var bit_shift = (0b1 << Self.bits) - 1

        var rec_len = len(record)
        if rec_len > self.max_length:
            self.max_length = rec_len
            self.hash_counts.resize(len(self.hash_list), self.max_length)

        # Check initial Kmer
        if len(self.hash_list) > 0:
            self._check_hashes(hash, 1)

        var seq_span = record.sequence()
        for i in range(end, rec_len):
            # Remove the most signifcant xx bits
            hash = hash & neg_mask

            # Mask for the least sig. three bits, add to hash
            var rem = seq_span[i] & bit_shift
            hash = (hash << Self.bits) + Int(rem)
            if len(self.hash_list) > 0:
                self._check_hashes(hash, i + 1)

    @always_inline
    fn _check_hashes(mut self, hash: UInt64, pos: Int):
        for i in range(len(self.hash_list)):
            if hash == self.hash_list[i]:
                self.hash_counts[Index(i, pos)] += 1


# ----- Adapter names (shared by summarizer output) -----
fn _adapter_names() -> List[String]:
    var names = List[String]()
    names.append("Illumina Universal Adapter")
    names.append("Illumina Small RNA 3' Adapter")
    names.append("Illumina Small RNA 5' Adapter")
    names.append("Nextera Transposase Sequence")
    names.append("PolyA")
    names.append("PolyG")
    return names^


# ----- Summarizer: prepare + grade + data block + plot -----

struct AdapterContentSummarizer(Summarizer, PlotOutput, Copyable, Movable):
    """Summarization and plotting for adapter content. Uses cached counts from collector (feed) and ctx.num_reads."""
    var _cache_hash_counts: Matrix2D[DType.int64]
    var _cache_kmer_len: Int
    var _cache_pct: Matrix2D[DType.float64]
    var _cache_status: String
    var _cache_ready: Bool

    fn __init__(out self):
        self._cache_hash_counts = Matrix2D[DType.int64](0, 0)
        self._cache_kmer_len = 0
        self._cache_pct = Matrix2D[DType.float64](0, 0)
        self._cache_status = "pass"
        self._cache_ready = False

    fn feed(mut self, collector: AdapterContentCollector[3]):
        """Copy collector state needed for summarization."""
        self._cache_hash_counts = collector.hash_counts.copy()
        self._cache_kmer_len = collector.kmer_len

    fn summerize(mut self, ctx: SummaryContext) raises:
        """Compute position-wise adapter percentages and status from cached counts."""
        var total_reads = ctx.num_reads
        if total_reads <= 0:
            self._cache_ready = True
            self._cache_status = "pass"
            return
        var rows = self._cache_hash_counts.rows
        var cols = self._cache_hash_counts.cols
        self._cache_pct = Matrix2D[DType.float64](rows, cols)
        var t = Float64(total_reads)
        for i in range(rows):
            for j in range(cols):
                self._cache_pct.set(
                    i, j, 100.0 * Float64(self._cache_hash_counts.get(i, j)) / t
                )
        self._cache_status = self._get_status(total_reads)
        self._cache_ready = True

    fn _get_status(self, total_reads: Int64) -> String:
        if total_reads == 0:
            return "pass"
        var max_length = self._cache_hash_counts.cols
        if max_length < self._cache_kmer_len:
            return "warn"
        var max_count: Int64 = 0
        for i in range(self._cache_hash_counts.rows):
            for j in range(self._cache_hash_counts.cols):
                var c = self._cache_hash_counts.get(i, j)
                if c > max_count:
                    max_count = c
        var pct = (Float64(max_count) / Float64(total_reads)) * 100.0
        if pct > ADAPTER_ERROR:
            return "fail"
        if pct > ADAPTER_WARN:
            return "warn"
        return "pass"

    fn grade(self) raises -> GradeEntry:
        return GradeEntry("Adapter Content", self._cache_status)

    fn data_block_body(self) raises -> String:
        if not self._cache_ready:
            return ""
        var adapter_names = _adapter_names()
        var out = String()
        out += "#Position"
        for i in range(min(6, self._cache_pct.rows)):
            out += "\t" + adapter_names[i]
        out += "\n"
        for j in range(self._cache_pct.cols):
            out += String(j + 1)
            for i in range(self._cache_pct.rows):
                out += "\t" + String(self._cache_pct.get(i, j))
            out += "\n"
        return out

    fn module_legend(self) -> String:
        return "Adapter Content"

    fn panel_id(self) -> String:
        return "adapter_content"

    fn plot_result(self) raises -> PythonObject:
        var plt = Python.import_module("matplotlib.pyplot")
        var np = Python.import_module("numpy")
        var arr = matrix_to_numpy(self._cache_pct)
        var arr_t = np.transpose(arr)
        var z = plt.subplots()
        var fig = z[0]
        var ax = z[1]
        ax.plot(arr_t)
        ax.set_ylim(0, 100)
        var legend_labels = Python.list()
        var names = _adapter_names()
        for i in range(len(names)):
            legend_labels.append(names[i])
        plt.legend(legend_labels)
        plt.xlabel("Position")
        plt.ylabel("Percentage of Reads")
        plt.title("Adapter content")
        return fig


# ----- Assembled module: Collector + Summarizer -----

@fieldwise_init
struct AdapterContentModule[bits: Int = 3](Copyable, Movable):
    """Module assembling AdapterContentCollector + AdapterContentSummarizer."""
    var collector: AdapterContentCollector[Self.bits]
    var summarizer: AdapterContentSummarizer

    fn __init__(out self, var hashes: List[UInt64], kmer_len: Int = 0):
        self.collector = AdapterContentCollector[Self.bits](hashes^, kmer_len)
        self.summarizer = AdapterContentSummarizer()

    @always_inline
    fn tally_read(mut self, record: FastqRecord, read_no: Int64):
        self.collector.tally_read(record, read_no)

    @always_inline
    fn tally_read(mut self, record: RefRecord, read_no: Int64):
        self.collector.tally_read(record, read_no)

    fn to_data_text(self, ctx: SummaryContext) raises -> String:
        var out = DefaultOutputter()
        var body = self.summarizer.data_block_body()
        var g = self.summarizer.grade()
        return out.wrap_data_block(self.summarizer.module_legend(), g.grade, body)

    fn to_html_panels(self) raises -> List[result_panel]:
        var panels = List[result_panel]()
        var out = DefaultOutputter()
        var fig = self.summarizer.plot_result()
        var panel = out.make_panel(
            self.summarizer.panel_id(),
            self.summarizer.grade().grade,
            self.summarizer.module_legend(),
            fig,
        )
        panels.append(panel^)
        return panels^

    fn plot(self, total_reads: Int64) raises -> PythonObject:
        """Legacy: plot from summarizer cache. Call feed + summerize(ctx) first."""
        return self.summarizer.plot_result()
