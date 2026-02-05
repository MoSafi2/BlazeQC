"""Quality distribution (split from stats_.mojo)."""

from python import Python, PythonObject
from blazeseq import FastqRecord, RecordCoord
from blazeqc.stats.analyser import Analyser, py_lib
from blazeqc.helpers import (
    Matrix2D,
    QualitySchema,
    matrix_to_numpy,
    grow_matrix,
    tensor_to_numpy_1d,
    make_linear_base_groups,
    bin_array,
    encode_img_b64,
)
from blazeqc.CONSTS import (
    illumina_1_5_schema,
    illumina_1_3_schema,
    illumina_1_8_schema,
    generic_schema,
)
from blazeqc.html_maker import result_panel


struct QualityDistribution(Analyser, Copyable, Movable):
    var qu_dist: Matrix2D
    var qu_dist_seq: List[Int64]
    var max_length: Int
    var max_qu: UInt8
    var min_qu: UInt8

    fn __init__(out self):
        self.qu_dist = Matrix2D(1, 128)
        self.qu_dist_seq = List[Int64](capacity=128)
        for _ in range(128):
            self.qu_dist_seq.append(0)
        self.max_length = 0
        self.max_qu = 0
        self.min_qu = 128

    fn tally_read(mut self, record: FastqRecord):
        if len(record) > self.max_length:
            self.max_length = len(record)
            var new_qu_dist = grow_matrix(self.qu_dist, self.max_length, 128)
            swap(self.qu_dist, new_qu_dist)

        for i in range(len(record)):
            var base_qu = record.QuStr[i]
            self.qu_dist.add(i, Int(base_qu), 1)
            if base_qu > self.max_qu:
                self.max_qu = base_qu
            if base_qu < self.min_qu:
                self.min_qu = base_qu

        var qu_sum: Int = 0
        for i in range(len(record)):
            qu_sum += Int(record.QuStr[i])
        var average = Int(qu_sum / len(record))
        while len(self.qu_dist_seq) <= average:
            self.qu_dist_seq.append(0)
        self.qu_dist_seq[average] += 1

    fn tally_read(mut self, record: RecordCoord):
        if record.qu_len() > self.max_length:
            self.max_length = Int(record.seq_len())
            var new_qu_dist = grow_matrix(self.qu_dist, self.max_length, 128)
            swap(self.qu_dist, new_qu_dist)

        for i in range(Int(record.qu_len())):
            var base_qu = record.QuStr[i]
            self.qu_dist.add(i, Int(base_qu), 1)
            if base_qu > self.max_qu:
                self.max_qu = base_qu
            if base_qu < self.min_qu:
                self.min_qu = base_qu
        var qu_sum: Int = 0
        for i in range(Int(record.qu_len())):
            qu_sum += Int(record.QuStr[i])
        var average = Int(qu_sum / record.qu_len())
        while len(self.qu_dist_seq) <= average:
            self.qu_dist_seq.append(0)
        self.qu_dist_seq[average] += 1

    # Use this answer for plotting: https://stackoverflow.com/questions/58053594/how-to-create-a-boxplot-from-data-with-weights
    fn slice_array(
        self, arr: PythonObject, min_index: Int, max_index: Int
    ) raises -> PythonObject:
        var np = Python.import_module("numpy")
        indices = np.arange(min_index, max_index)
        return np.take(arr, indices, axis=1)

    fn plot(self) raises -> Tuple[PythonObject, PythonObject]:
        Python.add_to_path(py_lib.as_string_slice())
        var np = Python.import_module("numpy")
        var plt = Python.import_module("matplotlib.pyplot")
        var mtp = Python.import_module("matplotlib")

        var schema = self._guess_schema()
        var arr = matrix_to_numpy(self.qu_dist)
        var min_index = schema.OFFSET
        var max_index = max(40, self.max_qu)
        arr = self.slice_array(arr, Int(min_index), Int(max_index))
        # Convert the raw array to binned array to account for very long seqs.
        var bins = make_linear_base_groups(arr.shape[0])
        arr, var py_bins = bin_array(arr, bins, func="mean")

        ################ Quality Boxplot ##################

        # TODO: Convert as much as possible away from numpy
        mean_line = np.sum(
            arr * np.arange(1, arr.shape[1] + 1), axis=1
        ) / np.sum(arr, axis=1)
        cum_sum = np.cumsum(arr, axis=1)
        total_counts = np.reshape(np.sum(arr, axis=1), (len(arr), 1))
        median = np.argmax(cum_sum > total_counts / 2, axis=1)
        Q75 = np.argmax(cum_sum > total_counts * 0.75, axis=1)
        Q25 = np.argmax(cum_sum > total_counts * 0.25, axis=1)
        IQR = Q75 - Q25

        whislo = np.full(len(IQR), None)
        whishi = np.full(len(IQR), None)

        x = plt.subplots()
        fig = x[0]
        ax = x[1]
        l = Python.list()
        for i in range(len(IQR)):
            var stat: PythonObject = Python.dict()
            stat["med"] = median[i]
            stat["q1"] = Q25[i]
            stat["q3"] = Q75[i]
            stat["whislo"] = whislo[i]
            stat["whishi"] = whishi[i]
            l.append(stat)

        ax.bxp(l, showfliers=False)
        ax.plot(mean_line)
        ax.set_ylim(0, 60)
        ax.set_title("Quality Scores across all bases")
        ax.set_xlabel("Position in read (bp)")

        bins_range = Python.list()
        for i in range(len(bins)):
            bins_range.append(i)

        ax.set_xticks(bins_range)
        ax.set_xticklabels(py_bins, rotation=45)
        ax.xaxis.set_major_locator(
            mtp.ticker.MaxNLocator(integer=True, nbins=15)
        )

        ###############################################################
        ####                Average quality /seq                   ####
        ###############################################################

        # Finding the last non-zero index
        index = 0
        for i in range(len(self.qu_dist_seq) - 1, -1, -1):
            if self.qu_dist_seq[i] != 0:
                index = i
                break

        arr2 = tensor_to_numpy_1d(self.qu_dist_seq)
        arr2 = arr2[Int(schema.OFFSET) : index + 2]
        z = plt.subplots()
        fig2 = z[0]
        ax2 = z[1]
        ax2.plot(arr2)
        ax2.set_xlabel("Mean Sequence Quality (Phred Score)")
        ax2.set_title("Quality score distribution over all sequences")

        return Tuple(fig, fig2)

    fn make_html(self) raises -> Tuple[result_panel, result_panel]:
        fig1, fig2 = self.plot()
        var encoded_fig1 = encode_img_b64(fig1)
        var encoded_fig2 = encode_img_b64(fig2)
        var result_1 = result_panel(
            "qu_score_dis_base",
            "pass",
            "Quality Scores Distribtion",
            encoded_fig1,
        )

        var result_2 = result_panel(
            "qu_score_dis_seq",
            "pass",
            "Mean Quality distribution",
            encoded_fig2,
        )

        return result_1, result_2

    fn _guess_schema(self) -> QualitySchema:
        comptime SANGER_ENCODING_OFFSET = 33
        comptime ILLUMINA_1_3_ENCODING_OFFSET = 64

        if self.min_qu < 64:
            return materialize(illumina_1_8_schema)
        elif self.min_qu == ILLUMINA_1_3_ENCODING_OFFSET + 1:
            return materialize(illumina_1_3_schema)
        elif self.min_qu <= 126:
            return materialize(illumina_1_5_schema)
        else:
            print("Unable to parse Quality Schema, returning generic schema")
            return materialize(generic_schema)
