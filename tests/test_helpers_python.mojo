"""Unit tests for blazeqc.helpers Python interop (NumPy, matplotlib, etc.)."""

from blazeqc.helpers import (
    Matrix2D,
    bin_array,
    encode_img_b64,
    from_numpy,
    list_float64_to_numpy,
    matrix_to_numpy,
    tensor_to_numpy_1d,
)
from python import Python, PythonObject
from testing import (
    assert_equal,
    assert_raises,
    assert_true,
    TestSuite,
)


# ----- tensor_to_numpy_1d -----


def test_tensor_to_numpy_1d():
    var data = List[Int64]()
    data.append(1)
    data.append(2)
    data.append(3)
    var arr = tensor_to_numpy_1d(data)
    var np = Python.import_module("numpy")
    assert_true(np.issubdtype(arr.dtype, np.integer))
    assert_equal(Int(py=arr.shape[0]), 3)
    assert_equal(Int64(py=arr.item(0)), 1)
    assert_equal(Int64(py=arr.item(1)), 2)
    assert_equal(Int64(py=arr.item(2)), 3)


# ----- list_float64_to_numpy -----


def test_list_float64_to_numpy():
    var data = List[Float64]()
    data.append(1.0)
    data.append(2.5)
    data.append(3.0)
    var arr = list_float64_to_numpy(data)
    var np = Python.import_module("numpy")
    assert_equal(Int(py=arr.shape[0]), 3)
    assert_equal(Float64(py=arr.item(0)), 1.0)
    assert_equal(Float64(py=arr.item(1)), 2.5)
    assert_equal(Float64(py=arr.item(2)), 3.0)


# ----- matrix_to_numpy -----


def test_matrix_to_numpy():
    var m = Matrix2D(2, 2)
    m.set(0, 0, 1)
    m.set(0, 1, 2)
    m.set(1, 0, 3)
    m.set(1, 1, 4)
    var arr = matrix_to_numpy(m)
    assert_equal(Int(py=arr.shape[0]), 2)
    assert_equal(Int(py=arr.shape[1]), 2)
    assert_equal(Int64(py=arr.item(0, 0)), 1)
    assert_equal(Int64(py=arr.item(0, 1)), 2)
    assert_equal(Int64(py=arr.item(1, 0)), 3)
    assert_equal(Int64(py=arr.item(1, 1)), 4)


# ----- from_numpy -----


def test_from_numpy_valid():
    var np = Python.import_module("numpy")
    var py_rows = Python.list()
    var row0 = Python.list()
    row0.append(1)
    row0.append(2)
    var row1 = Python.list()
    row1.append(3)
    row1.append(4)
    py_rows.append(row0)
    py_rows.append(row1)
    var arr = np.array(py_rows)
    var mat = from_numpy(arr)
    assert_equal(mat.shape()[0], 2)
    assert_equal(mat.shape()[1], 2)
    assert_equal(mat.get(0, 0), 1)
    assert_equal(mat.get(0, 1), 2)
    assert_equal(mat.get(1, 0), 3)
    assert_equal(mat.get(1, 1), 4)


def test_from_numpy_1d_raises():
    var np = Python.import_module("numpy")
    var py_list = Python.list()
    py_list.append(1)
    py_list.append(2)
    py_list.append(3)
    var arr_1d = np.array(py_list)
    with assert_raises(contains="2D"):
        _ = from_numpy(arr_1d)


# ----- bin_array -----


def test_bin_array_sum():
    var np = Python.import_module("numpy")
    var py_list = Python.list()
    py_list.append(1.0)
    py_list.append(2.0)
    py_list.append(3.0)
    py_list.append(4.0)
    var arr = np.array(py_list)
    var bins = List[Int]()
    bins.append(0)
    bins.append(2)
    bins.append(4)
    var binned, _ = bin_array(arr, bins, "sum")
    assert_equal(Int(py=binned.shape[0]), 3)
    assert_equal(Float64(py=binned.item(0)), 0.0)
    assert_equal(Float64(py=binned.item(1)), 3.0)
    assert_equal(Float64(py=binned.item(2)), 3.0)


def test_bin_array_mean():
    var np = Python.import_module("numpy")
    var py_list = Python.list()
    py_list.append(2.0)
    py_list.append(4.0)
    var arr = np.array(py_list)
    var bins = List[Int]()
    bins.append(0)
    bins.append(2)
    var binned, _ = bin_array(arr, bins, "mean")
    assert_equal(Int(py=binned.shape[0]), 2)


# ----- encode_img_b64 -----


def test_encode_img_b64():
    var matplotlib = Python.import_module("matplotlib")
    matplotlib.use("Agg")
    var plt = Python.import_module("matplotlib.pyplot")
    var fig = plt.figure()
    var result = encode_img_b64(fig)
    assert_true(len(result) > 0)
    assert_true(result.startswith("iVBOR") or len(result) > 100)


def main():
    TestSuite.discover_tests[__functions_in_module()]().run()
