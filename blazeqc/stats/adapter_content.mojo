"""Adapter content (split from stats_.mojo)."""

from utils import Index
from python import Python, PythonObject
from blazeseq import FastqRecord
from blazeqc.stats.analyser import py_lib
from blazeqc.helpers import grow_matrix, matrix_to_numpy, encode_img_b64
from blazeqc.html_maker import result_panel


# TODO: Check how to add the analyzer Trait again
# TODO: Also plot the Over-represented sequences.
# TODO: Add binning
@value
struct AdapterContent[bits: Int = 3]():
    var kmer_len: Int
    var hash_counts: Tensor[DType.int64]
    var hash_list: List[UInt64]
    var max_length: Int

    fn __init__(out self, hashes: List[UInt64], kmer_len: Int = 0):
        self.kmer_len = min(kmer_len, 64 // bits)
        self.hash_list = hashes
        shape = TensorShape(len(self.hash_list), 1)
        self.hash_counts = Tensor[DType.int64](shape)
        self.max_length = 0

    # TODO: Check if it will be easier to use the bool_tuple and hashes as a list instead
    @always_inline
    fn tally_read(mut self, record: FastqRecord, read_no: Int64):
        var hash: UInt64 = 0
        var end = 0
        # Make a custom bit mask of 1s by certain length
        var mask: UInt64 = (0b1 << self.kmer_len * bits) - 1
        var neg_mask = mask >> bits
        var bit_shift = (0b1 << bits) - 1

        if record.len_record() > self.max_length:
            self.max_length = record.len_record()
            new_shape = TensorShape(len(self.hash_list), self.max_length)
            self.hash_counts = grow_matrix(self.hash_counts, new_shape)

        # Check initial Kmer
        if len(self.hash_list) > 0:
            self._check_hashes(hash, 1)

        for i in range(end, record.len_record()):
            # Remove the most signifcant xx bits
            hash = hash & neg_mask

            # Mask for the least sig. three bits, add to hash
            var rem = record.SeqStr[i] & bit_shift
            hash = (hash << bits) + Int(rem)
            if len(self.hash_list) > 0:
                self._check_hashes(hash, i + 1)

    @always_inline
    fn plot(self, total_reads: Int64) raises -> PythonObject:
        Python.add_to_path(py_lib.as_string_slice())
        plt = Python.import_module("matplotlib.pyplot")
        arr = matrix_to_numpy(self.hash_counts)
        arr = (arr / total_reads) * 100

        z = plt.subplots()
        fig = z[0]
        ax = z[1]

        ax.plot(arr.T)
        ax.set_ylim(0, 100)
        plt.legend(
            [
                "Illumina Universal Adapter",
                "Illumina Small RNA 3' Adapter",
                "Illumina Small RNA 5' Adapter",
                "Nextera Transposase Sequence",
                "PolyA",
                "PolyG",
            ]
        )
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
        return result_1

    @always_inline
    fn _check_hashes(mut self, hash: UInt64, pos: Int):
        for i in range(len(self.hash_list)):
            if hash == self.hash_list[i]:
                self.hash_counts[Index(i, pos)] += 1
