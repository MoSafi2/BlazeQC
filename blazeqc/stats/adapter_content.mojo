"""Adapter content (split from stats_.mojo)."""

from utils import Index
from python import Python, PythonObject
from blazeqc.stats.analyser import Analyser
from blazeqc.helpers import (
    Matrix2D,
    matrix_to_numpy,
    encode_img_b64,
)
from blazeqc.html_maker import result_panel
from blazeqc.limits import ADAPTER_WARN, ADAPTER_ERROR
from blazeseq import FastqRecord, RefRecord


# TODO: Check how to add the analyzer Trait again
# TODO: Also plot the Over-represented sequences.
# TODO: Add binning
@fieldwise_init
struct AdapterContent[bits: Int = 3](Analyser):
    var kmer_len: Int
    var hash_counts: Matrix2D[DType.int64]
    var hash_list: List[UInt64]
    var max_length: Int
    var _cache_pct: Matrix2D[DType.float64]
    var _cache_ready: Bool

    fn __init__(out self, var hashes: List[UInt64], kmer_len: Int = 0):
        self.kmer_len = min(kmer_len, 64 // Self.bits)
        self.hash_list = hashes^
        self.hash_counts = Matrix2D[DType.int64](len(self.hash_list), 1)
        self.max_length = 0
        self._cache_pct = Matrix2D[DType.float64](0, 0)
        self._cache_ready = False

    fn __copyinit__(out self, existing: Self):
        self.kmer_len = existing.kmer_len
        self.hash_counts = Matrix2D[DType.int64](
            existing.hash_counts.rows, existing.hash_counts.cols
        )
        for i in range(existing.hash_counts.rows):
            for j in range(existing.hash_counts.cols):
                self.hash_counts.set(i, j, existing.hash_counts.get(i, j))
        self.hash_list = List[UInt64](capacity=len(existing.hash_list))
        for i in range(len(existing.hash_list)):
            self.hash_list.append(existing.hash_list[i])
        self.max_length = existing.max_length
        self._cache_pct = Matrix2D[DType.float64](
            existing._cache_pct.rows, existing._cache_pct.cols
        )
        for i in range(existing._cache_pct.rows):
            for j in range(existing._cache_pct.cols):
                self._cache_pct.set(i, j, existing._cache_pct.get(i, j))
        self._cache_ready = existing._cache_ready

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

        for i in range(end, rec_len):
            hash = hash & neg_mask
            var rem = record.sequence[i] & bit_shift
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

        for i in range(end, rec_len):
            # Remove the most signifcant xx bits
            hash = hash & neg_mask

            # Mask for the least sig. three bits, add to hash
            var rem = record.sequence.as_span()[i] & bit_shift
            hash = (hash << Self.bits) + Int(rem)
            if len(self.hash_list) > 0:
                self._check_hashes(hash, i + 1)

    @always_inline
    fn plot(self, total_reads: Int64) raises -> PythonObject:
        var plt = Python.import_module("matplotlib.pyplot")
        var np = Python.import_module("numpy")
        var arr: PythonObject = matrix_to_numpy(self.hash_counts)
        arr = (arr / total_reads) * 100
        var arr_t: PythonObject = np.transpose(arr)

        var z = plt.subplots()
        var fig: PythonObject = z[0]
        var ax: PythonObject = z[1]
        ax.plot(arr_t)
        ax.set_ylim(0, 100)
        var legend_labels = Python.list()
        legend_labels.append("Illumina Universal Adapter")
        legend_labels.append("Illumina Small RNA 3' Adapter")
        legend_labels.append("Illumina Small RNA 5' Adapter")
        legend_labels.append("Nextera Transposase Sequence")
        legend_labels.append("PolyA")
        legend_labels.append("PolyG")
        plt.legend(legend_labels)
        plt.xlabel("Position")
        plt.ylabel("Percentage of Reads")
        plt.title("Adapter content")

        return fig

    fn prepare_data(mut self, total_reads: Int64):
        """Compute position-wise adapter percentages and cache for get_module_data/data_plot."""
        if total_reads <= 0:
            return
        var rows = self.hash_counts.rows
        var cols = self.hash_counts.cols
        self._cache_pct = Matrix2D[DType.float64](rows, cols)
        var t = Float64(total_reads)
        for i in range(rows):
            for j in range(cols):
                self._cache_pct.set(
                    i, j, 100.0 * Float64(self.hash_counts.get(i, j)) / t
                )
        self._cache_ready = True

    fn get_module_data(self, total_reads: Int64) -> String:
        """Return FastQC-style block text for Adapter Content from cache (#Position + adapter columns)."""
        if not self._cache_ready:
            return ""
        var adapter_names = List[String]()
        adapter_names.append("Illumina Universal Adapter")
        adapter_names.append("Illumina Small RNA 3' Adapter")
        adapter_names.append("Illumina Small RNA 5' Adapter")
        adapter_names.append("Nextera Transposase Sequence")
        adapter_names.append("PolyA")
        adapter_names.append("PolyG")
        var out = ">>Adapter Content\t{}\n".format(self._get_status(total_reads))
        out += "#Position"
        for i in range(min(6, self._cache_pct.rows)):
            out += "\t" + adapter_names[i]
        out += "\n"
        for j in range(self._cache_pct.cols):
            out += String(j + 1)
            for i in range(self._cache_pct.rows):
                out += "\t" + String(self._cache_pct.get(i, j))
            out += "\n"
        out += ">>END_MODULE\n"
        return out

    fn data_plot(self, total_reads: Int64) raises -> PythonObject:
        """Plot from cache when ready; otherwise delegate to plot()."""
        if self._cache_ready:
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
            legend_labels.append("Illumina Universal Adapter")
            legend_labels.append("Illumina Small RNA 3' Adapter")
            legend_labels.append("Illumina Small RNA 5' Adapter")
            legend_labels.append("Nextera Transposase Sequence")
            legend_labels.append("PolyA")
            legend_labels.append("PolyG")
            plt.legend(legend_labels)
            plt.xlabel("Position")
            plt.ylabel("Percentage of Reads")
            plt.title("Adapter content")
            return fig
        return self.plot(total_reads)

    fn _get_status(self, total_reads: Int64) -> String:
        # Warn if read length too short to analyse adapter (kmer_len)
        if total_reads == 0:
            return "pass"
        if self.max_length < self.kmer_len:
            return "warn"
        var max_count: Int64 = 0
        for i in range(self.hash_counts.rows):
            for j in range(self.hash_counts.cols):
                var c = self.hash_counts.get(i, j)
                if c > max_count:
                    max_count = c
        var pct = (Float64(max_count) / Float64(total_reads)) * 100.0
        if pct > ADAPTER_ERROR:
            return "fail"
        if pct > ADAPTER_WARN:
            return "warn"
        return "pass"

    fn make_html(self, total_reads: Int64) raises -> result_panel:
        fig1 = self.data_plot(total_reads)
        var encoded_fig1 = encode_img_b64(fig1)
        var result_1 = result_panel(
            "adapter_content",
            self._get_status(total_reads),
            "Adapter Content",
            encoded_fig1,
        )
        return result_1^

    @always_inline
    fn _check_hashes(mut self, hash: UInt64, pos: Int):
        for i in range(len(self.hash_list)):
            if hash == self.hash_list[i]:
                self.hash_counts[Index(i, pos)] += 1
