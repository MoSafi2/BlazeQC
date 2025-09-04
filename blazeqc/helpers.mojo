from python import PythonObject, Python
from memory import UnsafePointer, memcpy

alias py_lib: String = "./.pixi/envs/default/lib/python3.12/site-packages/"


fn encode_img_b64(fig: PythonObject) raises -> String:
    Python.add_to_path(py_lib.as_string_slice())
    py_io = Python.import_module("io")
    py_base64 = Python.import_module("base64")
    plt = Python.import_module("matplotlib.pyplot")

    buf = py_io.BytesIO()
    plt.savefig(buf, format="png", dpi=150)
    buf.seek(0)
    plt.close()
    base64_image = py_base64.b64encode(buf.read()).decode("utf-8")
    buf.close()

    return String(base64_image)


# Same logic as exponential BaseGroups in FASTQC
# Ported from FASTQC: https://github.com/s-andrews/FastQC/blob/1faeea0412093224d7f6a07f777fad60a5650795/uk/ac/babraham/FastQC/Graphs/BaseGroup.java
fn get_exp_interval(max_len: Int) -> List[Int]:
    var pos: Int = 1
    var interval: Int = 1
    bins = List[Int]()
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

    return bins


# Ported from FastQC: https://github.com/s-andrews/FastQC/blob/1faeea0412093224d7f6a07f777fad60a5650795/uk/ac/babraham/FastQC/Graphs/BaseGroup.java
fn get_linear_interval(length: Int) -> Int:
    base_values = List[Int](2, 5, 10)
    multiplier = 1

    while True:
        for base in base_values:
            interval = base * multiplier
            group_count = 9 + ((length - 9) // interval)
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
    return ungrouped_list


fn make_linear_base_groups(max_length: Int) -> List[Int]:
    if max_length <= 75:
        return make_ungrouped_groups(max_length)

    interval = get_linear_interval(max_length)
    starting_base = 1
    groups = List[Int]()

    while starting_base <= max_length:
        end_base = starting_base + (interval - 1)

        if starting_base < 10:
            end_base = starting_base
        elif starting_base == 10 and interval > 10:
            end_base = interval - 1

        if end_base > max_length:
            end_base = max_length

        groups.append(starting_base)

        if starting_base < 10:
            starting_base += 1
        elif starting_base == 10 and interval > 10:
            starting_base = interval
        else:
            starting_base += interval

    return groups


@always_inline
fn bin_array(
    arr: PythonObject, bins: List[Int], func: String = "sum"
) raises -> Tuple[PythonObject, PythonObject]:
    Python.add_to_path(py_lib.as_string_slice())
    np = Python.import_module("numpy")

    py_bins = Python.list()
    for i in bins:
        py_bins.append(i[])
    x = np.linspace(0, arr.shape[0], arr.shape[0])
    binned_array = np.digitize(x, py_bins)
    py_binned_slices = Python.list()
    for i in range(len(bins)):
        mask = np.equal(binned_array, i)
        t1 = np.compress(mask, arr, axis=0)
        if func == "mean":
            t2 = np.mean(t1, axis=0)
        elif func == "sum":
            t2 = np.sum(t1, axis=0)
        else:
            t2 = t1
        py_binned_slices.append(t2)
    new_arr = np.vstack(py_binned_slices)
    return new_arr, py_bins


@fieldwise_init
struct QualitySchema(Copyable & Movable, Stringable, Writable):
    var SCHEMA: StringSlice[StaticConstantOrigin]
    var LOWER: UInt8
    var UPPER: UInt8
    var OFFSET: UInt8

    fn __init__(
        out self,
        schema: StringSlice[StaticConstantOrigin],
        lower: Int,
        upper: Int,
        offset: Int,
    ):
        self.SCHEMA = schema
        self.UPPER = upper
        self.LOWER = lower
        self.OFFSET = offset

    fn write_to[w: Writer](self, mut writer: w) -> None:
        writer.write(self.__str__())

    fn __str__(self) -> String:
        return (
            String("Quality schema: ")
            + self.SCHEMA
            + "\nLower: "
            + String(self.LOWER)
            + "\nUpper: "
            + String(self.UPPER)
            + "\nOffset: "
            + String(self.OFFSET)
        )


@always_inline
fn base2int(byte: Byte) -> UInt8:
    alias A_b = 65
    alias a_b = 95
    alias C_b = 67
    alias c_b = 97
    alias G_b = 71
    alias g_b = 103
    alias T_b = 84
    alias t_b = 116
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
    for i in range(0, len(seq)):
        # Remove the most signifcant 3 bits
        hash = hash & 0x1FFFFFFFFFFFFFFF
        # Mask for the least sig. three bits, add to hash
        var rem = ord(seq[i]) & 0b111
        hash = (hash << 3) + Int(rem)
    return hash


def tensor_to_numpy_1d[T: DType](tensor: Tensor[T]) -> PythonObject:
    Python.add_to_path(py_lib.as_string_slice())
    np = Python.import_module("numpy")
    ar = np.zeros(tensor.num_elements())
    for i in range(tensor.num_elements()):
        ar.itemset(i, tensor[i])
    return ar


def matrix_to_numpy[T: DType](tensor: Tensor[T]) -> PythonObject:
    np = Python.import_module("numpy")
    ar = np.zeros([tensor.shape()[0], tensor.shape()[1]])
    for i in range(tensor.shape()[0]):
        for j in range(tensor.shape()[1]):
            ar.itemset((i, j), tensor[i, j])
    return ar


fn grow_tensor[
    T: DType,
](old_tensor: Tensor[T], num_ele: Int) -> Tensor[T]:
    var new_tensor = Tensor[T](num_ele)
    cpy_tensor(new_tensor, old_tensor, old_tensor.num_elements(), 0, 0)
    return new_tensor


fn grow_matrix[
    T: DType
](old_tensor: Tensor[T], new_shape: TensorShape) -> Tensor[T]:
    var new_tensor = Tensor[T](new_shape)
    for i in range(old_tensor.shape()[0]):
        for j in range(old_tensor.shape()[1]):
            new_tensor[VariadicList(i, j)] = old_tensor[VariadicList(i, j)]
    return new_tensor


fn sum_tensor[T: DType](tensor: Tensor[T]) -> Int:
    acc = 0
    for i in range(tensor.num_elements()):
        acc += Int(tensor[i])
    return acc


# The Function does not provide bounds checks on purpose, the bounds checks is the callers responsibility
@always_inline
fn cpy_tensor[
    T: DType
    # simd_width: Int
](
    mut dest: Tensor[T],
    src: Tensor[T],
    num_elements: Int,
    dest_strt: Int = 0,
    src_strt: Int = 0,
):
    var dest_ptr: UnsafePointer[Scalar[T]] = dest._ptr + dest_strt
    var src_ptr: UnsafePointer[Scalar[T]] = src._ptr + src_strt
    memcpy(dest_ptr, src_ptr, num_elements)
