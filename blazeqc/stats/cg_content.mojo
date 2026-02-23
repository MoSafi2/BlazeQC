"""CG content (split from stats_.mojo)."""

from python import Python, PythonObject
from blazeseq import FastqRecord
from blazeqc.stats.analyser import Analyser
from blazeqc.helpers import tensor_to_numpy_1d, encode_img_b64
from blazeqc.html_maker import result_panel


# Done!
struct CGContent(Analyser, Copyable, Movable):
    var cg_content: List[Int64]
    var theoritical_distribution: List[Int64]

    fn __init__(out self):
        self.cg_content = List[Int64](capacity=101)
        for _ in range(101):
            self.cg_content.append(0)
        self.theoritical_distribution = List[Int64](capacity=101)
        for _ in range(101):
            self.theoritical_distribution.append(0)

    fn tally_read(mut self, record: FastqRecord):
        if len(record) == 0:
            return

        var cg_num = 0

        for index in range(0, len(record)):
            if (
                record.sequence[index] & 0b111 == 3
                or record.sequence[index] & 0b111 == 7
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

    fn plot(self) raises -> PythonObject:
        var plt = Python.import_module("matplotlib.pyplot")
        var arr = tensor_to_numpy_1d(self.cg_content)
        var theoritical_distribution = self.calculate_theoritical_distribution()
        var x = plt.subplots()
        var fig = x[0]
        var ax = x[1]
        ax.plot(arr, label="GC count per read")
        ax.plot(theoritical_distribution, label="Theoritical Distribution")
        ax.set_title("Per sequence GC content")
        ax.set_xlabel("Mean GC content (%)")

        return fig

    fn make_html(
        self,
    ) raises -> result_panel:
        fig = self.plot()
        var encoded_fig1 = encode_img_b64(fig)
        var result_1 = result_panel(
            "cg_content",
            "pass",
            "Per sequence GC content",
            encoded_fig1,
        )

        return result_1^
