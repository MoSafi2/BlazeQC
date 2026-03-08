"""BlazeQC helpers: matrix/list utilities, base grouping, quality schema, and Python/NumPy interop."""

from memory import UnsafePointer, memset_zero, memcpy, alloc
from python import PythonObject, Python
from utils import Index, IndexList

comptime py_lib: String = "./.pixi/envs/default/lib/python3.12/site-packages/"

# 2D matrix backed by a flat heap buffer (row-major).
# Generic over DType — use DType.int64 for count accumulators and
# DType.float64 for floating-point matrices (e.g. tile quality deviations).
struct Matrix2D[dtype: DType](Copyable, Movable):
    var data: UnsafePointer[Scalar[Self.dtype], MutExternalOrigin]
    var rows: Int
    var cols: Int

    fn __init__(out self, rows: Int, cols: Int):
        self.rows = rows
        self.cols = cols
        self.data = alloc[Scalar[Self.dtype]](rows * cols)
        memset_zero(self.data, rows * cols)

    fn __copyinit__(out self, other: Self):
        self.rows = other.rows
        self.cols = other.cols
        var n = other.rows * other.cols
        self.data = alloc[Scalar[Self.dtype]](n)
        memcpy(dest=self.data, src=other.data, count=n)

    fn __moveinit__(out self, deinit other: Self):
        self.rows = other.rows
        self.cols = other.cols
        self.data = other.data

    fn __del__(deinit self):
        if self.data:
            self.data.free()

    # --- shape / indexing ---

    fn shape(self) -> Tuple[Int, Int]:
        return (self.rows, self.cols)

    fn num_elements(self) -> Int:
        return self.rows * self.cols

    @always_inline
    fn _idx(self, i: Int, j: Int) -> Int:
        return i * self.cols + j

    @always_inline
    fn get(self, i: Int, j: Int) -> Scalar[Self.dtype]:
        return self.data[self._idx(i, j)]

    @always_inline
    fn set(mut self, i: Int, j: Int, v: Scalar[Self.dtype]):
        self.data[self._idx(i, j)] = v

    @always_inline
    fn add(mut self, i: Int, j: Int, delta: Scalar[Self.dtype]):
        self.data[self._idx(i, j)] += delta

    fn __getitem__(self, index: IndexList[2]) -> Scalar[Self.dtype]:
        return self.get(index[0], index[1])

    fn __setitem__(mut self, index: IndexList[2], value: Scalar[Self.dtype]):
        self.set(index[0], index[1], value)

    fn __eq__(self, other: Self) -> Bool:
        if self.rows != other.rows or self.cols != other.cols:
            return False
        for k in range(self.rows * self.cols):
            if self.data[k] != other.data[k]:
                return False
        return True

    # --- bulk operations ---

    fn fill(mut self, v: Scalar[Self.dtype]):
        """Set every element to v."""
        for k in range(self.rows * self.cols):
            self.data[k] = v

    fn zero(mut self):
        """Zero all elements."""
        memset_zero(self.data, self.rows * self.cols)

    fn row_sum(self, i: Int) -> Scalar[Self.dtype]:
        """Return the sum of row i."""
        var acc: Scalar[Self.dtype] = 0
        for j in range(self.cols):
            acc += self.data[self._idx(i, j)]
        return acc

    fn col_sum(self, j: Int) -> Scalar[Self.dtype]:
        """Return the sum of column j."""
        var acc: Scalar[Self.dtype] = 0
        for i in range(self.rows):
            acc += self.data[self._idx(i, j)]
        return acc

    fn resize(mut self, new_rows: Int, new_cols: Int):
        """Grow or shrink in-place, preserving the top-left overlap; new cells are zero."""
        var new_data = alloc[Scalar[Self.dtype]](new_rows * new_cols)
        memset_zero(new_data, new_rows * new_cols)
        var copy_rows = min(self.rows, new_rows)
        var copy_cols = min(self.cols, new_cols)
        for i in range(copy_rows):
            for j in range(copy_cols):
                new_data[i * new_cols + j] = self.data[i * self.cols + j]
        self.data.free()
        self.data = new_data
        self.rows = new_rows
        self.cols = new_cols

fn encode_img_b64(fig: PythonObject) raises -> String:
    """Encode a matplotlib figure as a base64 PNG string."""
    var py_io = Python.import_module("io")
    var py_base64 = Python.import_module("base64")
    var plt = Python.import_module("matplotlib.pyplot")

    var buf = py_io.BytesIO()
    plt.savefig(buf, format="png", dpi=150)
    buf.seek(0)
    plt.close()
    var base64_image = py_base64.b64encode(buf.read()).decode("utf-8")
    buf.close()

    return String(base64_image)


# Same logic as exponential BaseGroups in FASTQC
# Ported from FASTQC: https://github.com/s-andrews/FastQC/blob/1faeea0412093224d7f6a07f777fad60a5650795/uk/ac/babraham/FastQC/Graphs/BaseGroup.java
fn get_exp_interval(max_len: Int) -> List[Int]:
    var pos: Int = 1
    var interval: Int = 1
    var bins = List[Int]()
    while pos <= max_len:
        if pos > max_len:
            pos = max_len
        bins.append(pos)
        pos += interval
        if pos == 10 and max_len > 75:
            interval = 5
        if pos == 50 and max_len > 200:
            interval = 10
        if pos == 100 and max_len > 300:
            interval = 50
        if pos == 500 and max_len > 1000:
            interval = 100
        if pos == 1000 and max_len > 2000:
            interval = 500

    if bins[-1] < max_len:
        bins.append(max_len)

    return bins^


# Ported from FastQC: https://github.com/s-andrews/FastQC/blob/1faeea0412093224d7f6a07f777fad60a5650795/uk/ac/babraham/FastQC/Graphs/BaseGroup.java
fn get_linear_interval(length: Int) -> Int:
    var base_values: List[Int] = [2, 5, 10]
    var multiplier = 1

    while True:
        for base in base_values:
            var interval = base * multiplier
            var group_count = 9 + ((length - 9) // interval)
            if (length - 9) % interval != 0:
                group_count += 1

            if group_count < 75:
                return interval

        multiplier *= 10

        # if multiplier == 10_000_000:
        #     raise Error("Couldn't find a sensible interval grouping for length")


fn make_ungrouped_groups(max_length: Int) -> List[Int]:
    var ungrouped_list = List[Int]()
    for i in range(1, max_length + 1):
        ungrouped_list.append(i)
    return ungrouped_list^

fn make_linear_base_groups(max_length: Int) -> List[Int]:
    if max_length <= 75:
        return make_ungrouped_groups(max_length)

    var interval = get_linear_interval(max_length)
    var starting_base = 1
    var groups = List[Int]()

    while starting_base <= max_length:
        groups.append(starting_base)

        if starting_base < 10:
            starting_base += 1
        elif starting_base == 10 and interval > 10:
            starting_base = interval
        else:
            starting_base += interval

    return groups^


@always_inline
fn bin_array(
    arr: PythonObject, bins: List[Int], func: String = "sum"
) raises -> Tuple[PythonObject, PythonObject]:
    var np = Python.import_module("numpy")

    var py_bins = Python.list()
    for i in bins:
        py_bins.append(i)
    var x = np.linspace(0, arr.shape[0], arr.shape[0])
    var binned_array = np.digitize(x, py_bins)
    var py_binned_slices = Python.list()
    for i in range(1, len(bins) + 1):
        var mask = np.equal(binned_array, i)
        var t1 = np.compress(mask, arr, axis=0)
        var t2: PythonObject
        if func == "mean":
            t2 = np.mean(t1, axis=0)
        elif func == "sum":
            t2 = np.sum(t1, axis=0)
        else:
            t2 = t1
        py_binned_slices.append(t2)
    var new_arr = np.vstack(py_binned_slices)
    return new_arr, py_bins


@always_inline
fn base2int(byte: Byte) -> UInt8:
    comptime A_b = 65
    comptime a_b = 95
    comptime C_b = 67
    comptime c_b = 97
    comptime G_b = 71
    comptime g_b = 103
    comptime T_b = 84
    comptime t_b = 116
    if byte == A_b or byte == a_b:
        return 0
    if byte == C_b or byte == c_b:
        return 1
    if byte == G_b or byte == g_b:
        return 2
    if byte == T_b or byte == t_b:
        return 3
    return 4


# TODO: Make this also parametrized on the number of bits per bp
fn _seq_to_hash(seq: String) -> UInt64:
    var hash = 0
    var bytes = seq.as_bytes()
    for i in range(len(bytes)):
        # Remove the most significant 3 bits
        hash = hash & 0x1FFFFFFFFFFFFFFF
        # Mask for the least significant three bits, add to hash
        var rem = Int(bytes[i]) & 0b111
        hash = (hash << 3) + rem
    return hash


fn tensor_to_numpy_1d(ref data: List[Int64]) raises -> PythonObject:
    """Convert a 1D List[Int64] to a NumPy array."""
    var np = Python.import_module("numpy")
    var py_list = Python.list()
    for val in data:
        py_list.append(val)
    return np.array(py_list)


fn list_float64_to_numpy(ref data: List[Float64]) raises -> PythonObject:
    """Convert a List[Float64] to a NumPy float64 array."""
    var np = Python.import_module("numpy")
    var py_list = Python.list()
    for val in data:
        py_list.append(val)
    return np.array(py_list)


fn matrix_to_numpy[dtype: DType](m: Matrix2D[dtype]) raises -> PythonObject:
    """Convert a Matrix2D to a NumPy 2D array."""
    var np = Python.import_module("numpy")
    var py_rows = Python.list()
    for i in range(m.rows):
        var py_row = Python.list()
        for j in range(m.cols):
            py_row.append(m.get(i, j))
        py_rows.append(py_row)
    return np.array(py_rows)


fn from_numpy[dtype: DType = DType.int64](arr: PythonObject) raises -> Matrix2D[dtype]:
    """Create a Matrix2D[dtype] from a NumPy 2D array."""
    if arr.ndim != 2:
        raise Error("from_numpy expects a 2D array")
    var rows = Int(py=arr.shape[0])
    var cols = Int(py=arr.shape[1])
    var mat = Matrix2D[dtype](rows, cols)
    for i in range(rows):
        for j in range(cols):
            mat.set(i, j, Scalar[dtype](py=arr.item(i, j)))
    return mat^


fn grow_tensor(old_list: List[Int64], new_size: Int) -> List[Int64]:
    """Return a new list of length new_size with old data copied and remainder zero-filled."""
    var new_list = List[Int64](capacity=new_size)
    var old_len = len(old_list)
    for i in range(new_size):
        if i < old_len:
            new_list.append(old_list[i])
        else:
            new_list.append(0)
    return new_list^


fn sum_tensor(list: List[Int64]) -> Int64:
    """Sum all elements of a List[Int64]."""
    var acc: Int64 = 0
    for i in range(len(list)):
        acc += list[i]
    return acc
