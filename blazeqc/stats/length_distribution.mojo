"""Length distribution (split from stats_.mojo)."""

from python import Python, PythonObject
from blazeseq import FastqRecord, RecordCoord
from blazeqc.stats.analyser import Analyser, py_lib
from blazeqc.helpers import tensor_to_numpy_1d, bin_array, encode_img_b64
from blazeqc.html_maker import result_panel


struct LengthDistribution(Analyser, Copyable, Movable):
    var length_vector: List[Int64]

    fn __init__(out self):
        self.length_vector = List[Int64]()

    @always_inline
    fn tally_read(mut self, record: FastqRecord):
        while len(self.length_vector) < len(record):
            self.length_vector.append(0)
        self.length_vector[len(record) - 1] += 1

    @always_inline
    fn tally_read(mut self, record: RecordCoord):
        while len(self.length_vector) < Int(record.seq_len()):
            self.length_vector.append(0)
        self.length_vector[Int(record.seq_len() - 1)] += 1

    @always_inline
    fn length_average(self, num_reads: Int) -> Float64:
        var cum: Int64 = 0
        for i in range(len(self.length_vector)):
            cum += self.length_vector[i] * (i + 1)
        return Int(cum) / num_reads

    fn get_size_distribution(
        self, min_val: Int, max_val: Int
    ) -> Tuple[Int, Int]:
        # We won't group if they've asked us not to

        base = 1

        while base > (max_val - min_val):
            base //= 10

        var divisions: List[Int] = [1, 2, 5]

        while True:
            for d in divisions:
                tester = base * d
                if (max_val - min_val) / tester <= 50:
                    interval = tester
                    break
            else:
                base *= 10
                continue
            break

        # Now we work out the first value to be plotted
        basic_division = min_val // interval
        test_start = basic_division * interval
        starting = test_start

        return starting, interval

    fn plot(self) raises -> PythonObject:
        Python.add_to_path(py_lib.as_string_slice())
        var plt = Python.import_module("matplotlib.pyplot")
        var np = Python.import_module("numpy")
        var mtp = Python.import_module("matplotlib")

        var min_val: Int = 0
        var max_val: Int = len(self.length_vector)

        for i in range(len(self.length_vector)):
            if self.length_vector[i] > 0:
                min_val = i
                break

        _, interval = self.get_size_distribution(min_val, max_val)
        bins = List[Int]()
        start = 0
        while start <= max_val:
            bins.append(start)
            start += interval

        var arr = tensor_to_numpy_1d(self.length_vector)

        var x = plt.subplots()  # Create a figure
        var fig = x[0]
        var ax = x[1]

        var arr2 = np.insert(arr, 0, 0)
        var arr3 = np.append(arr2, 0)
        # bins = make_linear_base_groups(self.length_vector.num_elements())
        arr3, py_bins = bin_array(arr3, bins, func="mean")

        ticks = Python.list()
        for i in range(len(bins)):
            ticks.append(i)

        ax.plot(arr3)

        ax.set_xticks(ticks)
        ax.set_xticklabels(py_bins, rotation=45)
        ax.xaxis.set_major_locator(
            mtp.ticker.MaxNLocator(integer=True, nbins=15)
        )
        # ax.set_xlim(np.argmax(arr3 > 0) - 1, len(arr3) - 1)
        ax.set_ylim(0)
        ax.set_title("Distribution of sequence lengths over all sequences")
        ax.set_xlabel("Sequence Length (bp)")

        return fig

    fn make_html(self) raises -> result_panel:
        fig = self.plot()
        var encoded_fig1 = encode_img_b64(fig)
        var result_1 = result_panel(
            "seq_len_dis",
            "pass",
            "Sequence Duplication Levels",
            encoded_fig1,
        )

        return result_1^
