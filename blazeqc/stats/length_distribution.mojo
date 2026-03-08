"""Length distribution (split from stats_.mojo)."""

from python import Python, PythonObject
from blazeseq import FastqRecord, RefRecord
from blazeqc.stats.analyser import Analyser
from blazeqc.helpers import tensor_to_numpy_1d, bin_array, encode_img_b64
from blazeqc.html_maker import result_panel


struct LengthDistribution(Analyser, Copyable, Movable):
    var length_vector: List[Int64]
    var zero_length_count: Int
    var _cache_binned_arr: PythonObject
    var _cache_ticks: PythonObject
    var _cache_labels: PythonObject
    var _cache_xlim_left: Int
    var _cache_xlim_right: Int
    var _cache_ready: Bool

    fn __init__(out self) raises:
        self.length_vector = List[Int64]()
        self.zero_length_count = 0
        self._cache_binned_arr = Python.evaluate("None")
        self._cache_ticks = Python.evaluate("None")
        self._cache_labels = Python.evaluate("None")
        self._cache_xlim_left = 0
        self._cache_xlim_right = 0
        self._cache_ready = False

    @always_inline
    fn tally_read(mut self, record: FastqRecord):
        if len(record) == 0:
            self.zero_length_count += 1
            return
        while len(self.length_vector) < len(record):
            self.length_vector.append(0)
        self.length_vector[len(record) - 1] += 1

    @always_inline
    fn tally_read(mut self, record: RefRecord):
        if len(record) == 0:
            self.zero_length_count += 1
            return
        while len(self.length_vector) < len(record):
            self.length_vector.append(0)
        self.length_vector[len(record) - 1] += 1

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
        # if max_val <= min_val:
        #     return (min_val, 1)

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
        var plt = Python.import_module("matplotlib.pyplot")
        var np = Python.import_module("numpy")
        var mtp = Python.import_module("matplotlib")

        # Find actual min/max lengths (length_vector[i] = count of seqs with length i+1)
        var min_len: Int = 0
        var max_len: Int = len(self.length_vector)

        for i in range(len(self.length_vector)):
            if self.length_vector[i] > 0:
                min_len = i + 1  # actual length
                break

        # Add padding either side (matching Java: minLen--, maxLen++)
        if min_len > 0:
            min_len -= 1
        max_len += 1

        starting, interval = self.get_size_distribution(min_len, max_len)

        # Build bins starting at computed starting point (matching Java)
        bins = List[Int]()
        var start = starting
        while start <= max_len:
            bins.append(start)
            start += interval

        # Build x-axis category labels matching Java: "n" if interval==1, "n-m" otherwise
        var x_categories = List[String]()
        for i in range(len(bins)):
            var min_value = bins[i]
            var max_value = bins[i] + interval - 1
            if max_value > max_len:
                max_value = max_len
            if interval == 1:
                x_categories.append(String(min_value))
            else:
                x_categories.append(String(min_value) + "-" + String(max_value))

        # arr3[j] = count of seqs with length j (matches Java's lengthCounts[seqLen])
        var arr = tensor_to_numpy_1d(self.length_vector)
        var arr2 = np.insert(arr, 0, 0)   # prepend so arr2[1] = count for length 1
        var arr3 = np.append(arr2, 0)      # append trailing zero for padding

        arr3, _ = bin_array(arr3, bins, func="sum")

        var x = plt.subplots()
        var fig = x[0]
        var ax = x[1]

        var ticks = Python.list()
        var labels = Python.list()
        for i in range(len(x_categories)):
            ticks.append(i)
            labels.append(x_categories[i])

        ax.plot(arr3)

        ax.set_xticks(ticks)
        ax.set_xticklabels(labels, rotation=45)
        ax.xaxis.set_major_locator(
            mtp.ticker.MaxNLocator(integer=True, nbins=15)
        )
        ax.set_xlim(np.argmax(arr3 > 0) - 1, len(arr3))
        ax.set_ylim(0)
        ax.set_title("Distribution of sequence lengths over all sequences")
        ax.set_xlabel("Sequence Length (bp)")
        ax.set_ylabel("Number of Reads")

        return fig

    fn prepare_data(mut self) raises:
        """Compute binned length distribution and cache for data_plot; get_module_data uses length_vector directly."""
        var plt = Python.import_module("matplotlib.pyplot")
        var np = Python.import_module("numpy")

        var min_len: Int = 0
        var max_len: Int = len(self.length_vector)
        for i in range(len(self.length_vector)):
            if self.length_vector[i] > 0:
                min_len = i + 1
                break
        if min_len > 0:
            min_len -= 1
        max_len += 1

        var starting = 0
        var interval = 1
        starting, interval = self.get_size_distribution(min_len, max_len)

        var bins = List[Int]()
        var start = starting
        while start <= max_len:
            bins.append(start)
            start += interval

        var x_categories = List[String]()
        for i in range(len(bins)):
            var min_value = bins[i]
            var max_value = bins[i] + interval - 1
            if max_value > max_len:
                max_value = max_len
            if interval == 1:
                x_categories.append(String(min_value))
            else:
                x_categories.append(String(min_value) + "-" + String(max_value))

        var arr = tensor_to_numpy_1d(self.length_vector)
        var arr2 = np.insert(arr, 0, 0)
        var arr3 = np.append(arr2, 0)
        arr3, _ = bin_array(arr3, bins, func="sum")

        var ticks = Python.list()
        var labels = Python.list()
        for i in range(len(x_categories)):
            ticks.append(i)
            labels.append(x_categories[i])

        var xlim_left = Int(py=np.argmax(arr3 > 0)) - 1
        var xlim_right = len(arr3)

        self._cache_binned_arr = arr3
        self._cache_ticks = ticks
        self._cache_labels = labels
        self._cache_xlim_left = xlim_left
        self._cache_xlim_right = xlim_right
        self._cache_ready = True

    fn get_module_data(self) -> String:
        """Return FastQC-style block text for Sequence Length Distribution (#Length, Count)."""
        var out = ">>Sequence Length Distribution\t{}\n".format(self._get_status())
        out += "#Length\tCount\n"
        for i in range(len(self.length_vector)):
            if self.length_vector[i] > 0:
                out += "{}\t{}\n".format(i + 1, self.length_vector[i])
        out += ">>END_MODULE\n"
        return out

    fn data_plot(self) raises -> PythonObject:
        """Plot from cache when ready; otherwise delegate to plot()."""
        if self._cache_ready:
            var plt = Python.import_module("matplotlib.pyplot")
            var np = Python.import_module("numpy")
            var mtp = Python.import_module("matplotlib")
            var x = plt.subplots()
            var fig = x[0]
            var ax = x[1]
            ax.plot(self._cache_binned_arr)
            ax.set_xticks(self._cache_ticks)
            ax.set_xticklabels(self._cache_labels, rotation=45)
            ax.xaxis.set_major_locator(
                mtp.ticker.MaxNLocator(integer=True, nbins=15)
            )
            ax.set_xlim(self._cache_xlim_left, self._cache_xlim_right)
            ax.set_ylim(0)
            ax.set_title("Distribution of sequence lengths over all sequences")
            ax.set_xlabel("Sequence Length (bp)")
            ax.set_ylabel("Number of Reads")
            return fig
        return self.plot()

    fn _get_status(self) -> String:
        # Error if any zero-length sequences exist (matching FastQC raisesError)
        if self.zero_length_count > 0:
            return "fail"
        # Warning if sequences have more than one distinct length (matching FastQC raisesWarning)
        var distinct_lengths: Int = 0
        for i in range(len(self.length_vector)):
            if self.length_vector[i] > 0:
                distinct_lengths += 1
        if distinct_lengths > 1:
            return "warn"
        return "pass"

    fn make_html(self) raises -> result_panel:
        fig = self.data_plot()
        var encoded_fig1 = encode_img_b64(fig)
        var result_1 = result_panel(
            "seq_len_dis",
            self._get_status(),
            "Sequence Length Distribution",
            encoded_fig1,
        )

        return result_1^
