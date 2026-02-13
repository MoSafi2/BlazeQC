"""Adapter content (split from stats_.mojo)."""

from utils import Index
from python import Python, PythonObject
from blazeqc.stats.analyser import Analyser
from blazeqc.helpers import (
    Matrix2D,
    grow_matrix,
    matrix_to_numpy,
    encode_img_b64,
)
from blazeqc.html_maker import result_panel
#from blazeseq.record import FastqRecord


# TODO: Check how to add the analyzer Trait again
# TODO: Also plot the Over-represented sequences.
# TODO: Add binning
@fieldwise_init
struct AdapterContent[bits: Int = 3](Analyser):
    var kmer_len: Int
    var hash_counts: Matrix2D
    var hash_list: List[UInt64]
    var max_length: Int

    fn __init__(out self, var hashes: List[UInt64], kmer_len: Int = 0):
        self.kmer_len = min(kmer_len, 64 // Self.bits)
        self.hash_list = hashes^
        self.hash_counts = Matrix2D(len(self.hash_list), 1)
        self.max_length = 0

    fn __copyinit__(out self, existing: Self):
        self.kmer_len = existing.kmer_len
        self.hash_counts = Matrix2D(
            existing.hash_counts.rows, existing.hash_counts.cols
        )
        for i in range(existing.hash_counts.rows):
            for j in range(existing.hash_counts.cols):
                self.hash_counts.set(i, j, existing.hash_counts.get(i, j))
        self.hash_list = List[UInt64](capacity=len(existing.hash_list))
        for i in range(len(existing.hash_list)):
            self.hash_list.append(existing.hash_list[i])
        self.max_length = existing.max_length

    @always_inline
    fn tally_read(mut self, record: FastqRecord):
        self.tally_read(record, 0)

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
            self.hash_counts = grow_matrix(
                self.hash_counts, len(self.hash_list), self.max_length
            )

        # Check initial Kmer
        if len(self.hash_list) > 0:
            self._check_hashes(hash, 1)

        for i in range(end, rec_len):
            # Remove the most signifcant xx bits
            hash = hash & neg_mask

            # Mask for the least sig. three bits, add to hash
            var rem = record.SeqStr[i] & bit_shift
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
        plt.title("Adapter Content")

        return fig

    fn make_html(self, total_reads: Int64) raises -> result_panel:
        fig1 = self.plot(total_reads)
        var encoded_fig1 = encode_img_b64(fig1)
        var result_1 = result_panel(
            "adapter_content",
            "pass",
            "Adapter Content",
            encoded_fig1,
        )
        return result_1^

    @always_inline
    fn _check_hashes(mut self, hash: UInt64, pos: Int):
        for i in range(len(self.hash_list)):
            if hash == self.hash_list[i]:
                self.hash_counts[Index(i, pos)] += 1
