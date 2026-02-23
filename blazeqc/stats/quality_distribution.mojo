"""Quality distribution (split from stats_.mojo)."""

from python import Python, PythonObject
from blazeseq import FastqRecord, RefRecord
from blazeqc.stats.analyser import Analyser
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
            var base_qu = record.quality[i]
            self.qu_dist.add(i, Int(base_qu), 1)
            if base_qu > self.max_qu:
                self.max_qu = base_qu
            if base_qu < self.min_qu:
                self.min_qu = base_qu

        var qu_sum: Int = 0
        for i in range(len(record)):
            qu_sum += Int(record.quality[i])
        var average = Int(qu_sum / len(record))
        while len(self.qu_dist_seq) <= average:
            self.qu_dist_seq.append(0)
        self.qu_dist_seq[average] += 1

    fn tally_read(mut self, record: RefRecord):
        if len(record) > self.max_length:
            self.max_length = len(record)
            var new_qu_dist = grow_matrix(self.qu_dist, self.max_length, 128)
            swap(self.qu_dist, new_qu_dist)

        for i in range(len(record)):
            var base_qu = record.quality[i]
            self.qu_dist.add(i, Int(base_qu), 1)
            if base_qu > self.max_qu:
                self.max_qu = base_qu
            if base_qu < self.min_qu:
                self.min_qu = base_qu

        var qu_sum: Int = 0
        for i in range(len(record)):
            qu_sum += Int(record.quality[i])
        var average = Int(qu_sum / len(record))
        while len(self.qu_dist_seq) <= average:
            self.qu_dist_seq.append(0)
        self.qu_dist_seq[average] += 1

    # Use this answer for plotting: https://stackoverflow.com/questions/58053594/how-to-create-a-boxplot-from-data-with-weights
    fn slice_array(
        self, arr: PythonObject, min_index: Int, max_index: Int
    ) raises -> PythonObject:
        var np = Python.import_module("numpy")
        var indices = np.arange(min_index, max_index)
        return np.take(arr, indices, axis=1)

    fn plot(self) raises -> Tuple[PythonObject, PythonObject]:
        var np = Python.import_module("numpy")
        var plt = Python.import_module("matplotlib.pyplot")
        var mtp = Python.import_module("matplotlib")

        var schema = self._guess_schema()
        var arr = matrix_to_numpy(self.qu_dist)
        var min_index = schema.OFFSET
        var max_qu_int = Int(self.max_qu)
        var max_index = 40 if 40 > max_qu_int else max_qu_int
        arr = self.slice_array(arr, Int(min_index), Int(max_index))
        # Convert the raw array to binned array to account for very long seqs.
        var bins = make_linear_base_groups(self.max_length)
        var py_bins: PythonObject
        arr, py_bins = bin_array(arr, bins, func="mean")

        ################ Quality Boxplot ##################

        # TODO: Convert as much as possible away from numpy
        var nrows = Int(py=arr.shape[0])
        var ncols = Int(py=arr.shape[1])
        var mean_line = np.sum(
            arr * np.arange(1, ncols + 1), axis=1
        ) / np.sum(arr, axis=1)
        var cum_sum = np.cumsum(arr, axis=1)
        var total_counts = np.reshape(np.sum(arr, axis=1), Python.tuple(nrows, 1))
        var median = np.argmax(cum_sum > total_counts / 2, axis=1)
        var Q75 = np.argmax(cum_sum > total_counts * 0.75, axis=1)
        var Q25 = np.argmax(cum_sum > total_counts * 0.25, axis=1)
        var IQR = Q75 - Q25

        # Get Python None for boxplot default whiskers (matplotlib bxp uses None to auto whiskers)
        var py_none = Python.evaluate("None")
        var whislo = np.full(len(IQR), py_none)
        var whishi = np.full(len(IQR), py_none)

        var x = plt.subplots()
        var fig = x[0]
        var ax = x[1]
        var l = Python.list()
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

        var bins_range = Python.list()
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
        var index = 0
        for i in range(len(self.qu_dist_seq) - 1, -1, -1):
            if self.qu_dist_seq[i] != 0:
                index = i
                break

        var arr2 = tensor_to_numpy_1d(self.qu_dist_seq)
        arr2 = arr2[Int(schema.OFFSET) : index + 2]
        var z = plt.subplots()
        var fig2 = z[0]
        var ax2 = z[1]
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

        return (result_1^, result_2^)

    fn _guess_schema(self) -> QualitySchema:
        comptime SANGER_ENCODING_OFFSET = 33
        comptime ILLUMINA_1_3_ENCODING_OFFSET = 64

        if self.min_qu < 64:
            return QualitySchema("Illumina v1.8", 33, 126, 33)
        elif self.min_qu == ILLUMINA_1_3_ENCODING_OFFSET + 1:
            return QualitySchema("Illumina v1.3", 64, 126, 64)
        elif self.min_qu <= 126:
            return QualitySchema("Illumina v1.5", 66, 126, 64)
        else:
            print("Unable to parse Quality Schema, returning generic schema")
            return QualitySchema("Generic", 33, 126, 33)
