"""Base pair distribution (split from stats_.mojo)."""
from collections.dict import Dict
from python import Python, PythonObject
from blazeseq import FastqRecord, RefRecord
from blazeqc.stats.analyser import Analyser
from blazeqc.helpers import (
    Matrix2D,
    matrix_to_numpy,
    make_linear_base_groups,
    bin_array,
    encode_img_b64,
)
from blazeqc.html_maker import result_panel
from blazeqc.limits import N_CONTENT_WARN, N_CONTENT_ERROR, SEQUENCE_WARN, SEQUENCE_ERROR


@fieldwise_init
struct BasepairDistribution(Analyser):
    var bp_dist: Matrix2D[DType.int64]
    var max_length: Int
    var min_length: Int
    comptime WIDTH = 5

    fn __init__(out self):
        self.bp_dist = Matrix2D[DType.int64](1, self.WIDTH)
        self.max_length = 0
        self.min_length = Int.MAX

    fn tally_read(mut self, record: FastqRecord):
        var rec_len = len(record)
        if rec_len > self.max_length:
            self.max_length = rec_len
            self.bp_dist.resize(self.max_length, self.WIDTH)

        if rec_len < self.min_length:
            self.min_length = rec_len

        for i in range(rec_len):
            var base_val = Int((record.sequence[i] & 0b11111) % self.WIDTH)
            self.bp_dist.add(i, base_val, 1)

    fn tally_read(mut self, record: RefRecord):
        var rec_len = len(record)
        if rec_len > self.max_length:
            self.max_length = rec_len
            self.bp_dist.resize(self.max_length, self.WIDTH)

        if rec_len < self.min_length:
            self.min_length = rec_len

        for i in range(rec_len):
            var base_val = Int((record.sequence[i] & 0b11111) % self.WIDTH)
            self.bp_dist.add(i, base_val, 1)

    fn _plot_sequence_content(
        self, arr1: PythonObject, py_bins: PythonObject, bins: List[Int]
    ) raises -> PythonObject:
        """Plot per-base sequence content (C, G, T, A) across all bases."""
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

    fn _plot_n_content(
        self, arr2: PythonObject, py_bins: PythonObject, bins: List[Int]
    ) raises -> PythonObject:
        """Plot per-base N content across all bases."""
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

    fn plot(
        self, total_reads: Int64
    ) raises -> Tuple[PythonObject, PythonObject]:
        var np = Python.import_module("numpy")
        var arr = matrix_to_numpy(self.bp_dist)
        var bins = make_linear_base_groups(self.max_length)
        arr, py_bins = bin_array(arr, bins, func="sum")
        arr = (np.divide(arr.T, arr.sum(axis=1)).T) * 100
        var arr1 = arr[:, 0:4]  # C,G,T,A pairs
        var arr2 = arr[:, 4:5]  # N content
        # Return (N figure, base figure) so make_html receives fig1=N, fig2=base
        return (
            self._plot_n_content(arr2, py_bins, bins),
            self._plot_sequence_content(arr1, py_bins, bins),
        )

    fn _get_status_n_content(self) -> String:
        """Max N% at any position (FastQC n_content)."""
        var max_n_pct: Float64 = 0.0
        for row in range(self.bp_dist.rows):
            var total: Int64 = 0
            for col in range(self.bp_dist.cols):
                total += self.bp_dist.get(row, col)
            if total == 0:
                continue
            var n_pct = (Float64(self.bp_dist.get(row, 4)) / Float64(total)) * 100.0
            if n_pct > max_n_pct:
                max_n_pct = n_pct
        if max_n_pct > N_CONTENT_ERROR:
            return "fail"
        if max_n_pct > N_CONTENT_WARN:
            return "warn"
        return "pass"

    fn _get_status_sequence_content(self) -> String:
        """Max deviation A-T or C-G at any position (FastQC per-base sequence content)."""
        var max_dev: Float64 = 0.0
        for row in range(self.bp_dist.rows):
            var total: Int64 = 0
            for col in range(self.bp_dist.cols):
                total += self.bp_dist.get(row, col)
            if total == 0:
                continue
            var t = Float64(total)
            var pct_a = (Float64(self.bp_dist.get(row, 3)) / t) * 100.0
            var pct_t = (Float64(self.bp_dist.get(row, 2)) / t) * 100.0
            var pct_c = (Float64(self.bp_dist.get(row, 0)) / t) * 100.0
            var pct_g = (Float64(self.bp_dist.get(row, 1)) / t) * 100.0
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

    @always_inline
    fn make_grade(self, grades: Dict[String, Int]):
        pass

    fn make_html(
        self, total_reads: Int64
    ) raises -> Tuple[result_panel, result_panel]:
        fig1, fig2 = self.plot(total_reads)
        var encoded_fig1 = encode_img_b64(fig1)
        var encoded_fig2 = encode_img_b64(fig2)
        # fig1 = N figure, fig2 = base figure (see plot() return order)
        var result_1 = result_panel(
            "n_percentage",
            self._get_status_n_content(),
            "Per Base N Content",
            encoded_fig1,
        )
        var result_2 = result_panel(
            "base_pair_distribution",
            self._get_status_sequence_content(),
            "Per Base Sequence Content",
            encoded_fig2,
        )

        return (result_1^, result_2^)
