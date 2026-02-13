"""Base pair distribution (split from stats_.mojo)."""
from collections.dict import Dict
from python import Python, PythonObject
from blazeseq import FastqRecord, RecordCoord
from blazeqc.stats.analyser import Analyser
from blazeqc.helpers import (
    Matrix2D,
    matrix_to_numpy,
    grow_matrix,
    make_linear_base_groups,
    bin_array,
    encode_img_b64,
)
from blazeqc.html_maker import result_panel


@fieldwise_init
struct BasepairDistribution(Analyser):
    var bp_dist: Matrix2D
    var max_length: Int
    var min_length: Int
    comptime WIDTH = 5

    fn __init__(out self):
        self.bp_dist = Matrix2D(1, self.WIDTH)
        self.max_length = 0
        self.min_length = Int.MAX

    fn tally_read(mut self, record: FastqRecord):
        var rec_len = len(record)
        if rec_len > self.max_length:
            self.max_length = rec_len
            var new_mat = grow_matrix(self.bp_dist, self.max_length, self.WIDTH)
            swap(self.bp_dist, new_mat)

        if rec_len < self.min_length:
            self.min_length = rec_len

        for i in range(rec_len):
            var base_val = Int((record.SeqStr[i] & 0b11111) % self.WIDTH)
            self.bp_dist.add(i, base_val, 1)

    fn tally_read(mut self, record: RecordCoord):
        if record.seq_len() > self.max_length:
            self.max_length = Int(record.seq_len())
            var new_mat = grow_matrix(self.bp_dist, self.max_length, self.WIDTH)
            swap(self.bp_dist, new_mat)

        for i in range(Int(record.seq_len())):
            var base_val = Int((record.SeqStr[i] & 0b11111) % self.WIDTH)
            self.bp_dist.add(i, base_val, 1)

    fn plot(
        self, total_reads: Int64
    ) raises -> Tuple[PythonObject, PythonObject]:
        var plt = Python.import_module("matplotlib.pyplot")
        var mtp = Python.import_module("matplotlib")
        var np = Python.import_module("numpy")
        var arr = matrix_to_numpy(self.bp_dist)
        var bins = make_linear_base_groups(self.max_length)
        arr, py_bins = bin_array(arr, bins, func="sum")
        arr = (np.divide(arr.T, arr.sum(axis=1)).T) * 100

        var arr1 = arr[:, 0:4]  # C,G,T,A pairs
        var arr2 = arr[:, 4:5]  # N content
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
        ax.set_title("Base Distribution")

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
        ax2.set_title("N percentage")

        return fig, fig2

    @always_inline
    fn make_grade(self, grades: Dict[String, Int]):
        pass

    fn make_html(
        self, total_reads: Int64
    ) raises -> Tuple[result_panel, result_panel]:
        fig1, fig2 = self.plot(total_reads)
        var encoded_fig1 = encode_img_b64(fig1)
        var encoded_fig2 = encode_img_b64(fig2)
        var result_1 = result_panel(
            "base_pair_distribution",
            "pass",
            "Base Pair Distribtion",
            encoded_fig1,
        )
        var result_2 = result_panel(
            "n_percentage",
            "pass",
            "N Percentage (%)",
            encoded_fig2,
        )

        return (result_1^, result_2^)
