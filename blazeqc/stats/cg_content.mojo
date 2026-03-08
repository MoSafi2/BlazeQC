"""CG content (split from stats_.mojo)."""

from python import Python, PythonObject
from blazeseq import FastqRecord, RefRecord
from blazeqc.stats.analyser import Analyser
from blazeqc.helpers import tensor_to_numpy_1d, encode_img_b64
from blazeqc.html_maker import result_panel
from blazeqc.limits import GC_SEQUENCE_WARN, GC_SEQUENCE_ERROR


# Done!
struct CGContent(Analyser, Copyable, Movable):
    var cg_content: List[Int64]
    var theoritical_distribution: List[Int64]
    var _cache_theoretical: PythonObject
    var _cache_ready: Bool

    fn __init__(out self) raises:
        self.cg_content = List[Int64](capacity=101)
        for _ in range(101):
            self.cg_content.append(0)
        self.theoritical_distribution = List[Int64](capacity=101)
        for _ in range(101):
            self.theoritical_distribution.append(0)
        self._cache_theoretical = Python.evaluate("None")
        self._cache_ready = False

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

    # TODO: Convert as much as possible away from numpy
    fn calculate_theoritical_distribution(self) raises -> PythonObject:
        var np = Python.import_module("numpy")
        var sc = Python.import_module("scipy")
        var arr = tensor_to_numpy_1d(self.cg_content)
        var total_counts = np.sum(arr)
        var x_categories = np.arange(len(arr))
        var mode = np.argmax(arr)

        var stdev = np.sqrt(
            np.sum((x_categories - mode) ** 2 * arr) / (total_counts - 1)
        )
        var nd = sc.stats.norm(loc=mode, scale=stdev)
        var theoritical_distribution = nd.pdf(x_categories) * total_counts
        return theoritical_distribution

    fn _max_gc_deviation(self) raises -> Float64:
        """Max absolute deviation (%) between observed and theoretical GC distribution."""
        var theor = self.calculate_theoritical_distribution()
        var np = Python.import_module("numpy")
        var obs_arr = tensor_to_numpy_1d(self.cg_content)
        var total_obs = Float64(py=np.sum(obs_arr))
        var total_theor = Float64(py=np.sum(theor))
        if total_obs <= 0 or total_theor <= 0:
            return 0.0
        var max_dev: Float64 = 0.0
        for i in range(len(self.cg_content)):
            var o_pct = (Float64(self.cg_content[i]) / total_obs) * 100.0
            var t_val = Float64(py=theor[Int(i)])
            var t_pct = (t_val / total_theor) * 100.0
            var dev = o_pct - t_pct
            if dev < 0:
                dev = -dev
            if dev > max_dev:
                max_dev = dev
        return max_dev

    fn prepare_data(mut self) raises:
        """Compute theoretical distribution and cache for get_module_data/data_plot."""
        self._cache_theoretical = self.calculate_theoritical_distribution()
        self._cache_ready = True

    fn get_module_data(self) raises -> String:
        """Return FastQC-style block text for Per Sequence GC Content from cache (cg_content counts)."""
        if not self._cache_ready:
            return ""
        var out = ">>Per Sequence GC Content\t{}\n".format(self._get_status())
        out += "#GC Content\tCount\n"
        for i in range(len(self.cg_content)):
            out += "{}\t{}\n".format(i, self.cg_content[i])
        out += ">>END_MODULE\n"
        return out

    fn data_plot(self) raises -> PythonObject:
        """Plot from cache when ready; otherwise delegate to plot()."""
        if self._cache_ready:
            var plt = Python.import_module("matplotlib.pyplot")
            var arr = tensor_to_numpy_1d(self.cg_content)
            var x = plt.subplots()
            var fig = x[0]
            var ax = x[1]
            ax.plot(arr, label="GC count per read")
            ax.plot(self._cache_theoretical, label="Theoritical Distribution")
            ax.set_title("GC distribution over all sequences")
            ax.set_xlabel("Mean GC content (%)")
            return fig
        return self.plot()

    fn _get_status(self) raises -> String:
        var max_dev = self._max_gc_deviation()
        if max_dev > GC_SEQUENCE_ERROR:
            return "fail"
        if max_dev > GC_SEQUENCE_WARN:
            return "warn"
        return "pass"

    fn plot(self) raises -> PythonObject:
        var plt = Python.import_module("matplotlib.pyplot")
        var arr = tensor_to_numpy_1d(self.cg_content)
        var theoritical_distribution = self.calculate_theoritical_distribution()
        var x = plt.subplots()
        var fig = x[0]
        var ax = x[1]
        ax.plot(arr, label="GC count per read")
        ax.plot(theoritical_distribution, label="Theoritical Distribution")
        ax.set_title("GC distribution over all sequences")
        ax.set_xlabel("Mean GC content (%)")

        return fig

    fn make_html(
        self,
    ) raises -> result_panel:
        fig = self.data_plot()
        var encoded_fig1 = encode_img_b64(fig)
        var result_1 = result_panel(
            "cg_content",
            self._get_status(),
            "Per Sequence GC Content",
            encoded_fig1,
        )

        return result_1^
