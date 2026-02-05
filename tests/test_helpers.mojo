"""Unit tests for blazeqc.helpers (pure Mojo APIs only)."""

from blazeqc.helpers import (
    Matrix2D,
    base2int,
    get_exp_interval,
    get_linear_interval,
    grow_matrix,
    grow_tensor,
    make_linear_base_groups,
    make_ungrouped_groups,
    sum_tensor,
    _seq_to_hash,
)
from testing import assert_equal, assert_true, TestSuite


# ----- 1. Matrix2D -----


def test_matrix2d_shape_and_num_elements():
    var m = Matrix2D(2, 3)
    var sh = m.shape()
    assert_equal(sh[0], 2)
    assert_equal(sh[1], 3)
    assert_equal(m.num_elements(), 6)


def test_matrix2d_init_zeroed():
    var m = Matrix2D(2, 3)
    assert_equal(m.get(0, 0), 0)
    assert_equal(m.get(1, 2), 0)
    assert_equal(m.get(0, 1), 0)


def test_matrix2d_get_set():
    var m = Matrix2D(2, 3)
    m.set(0, 0, 10)
    m.set(1, 2, -5)
    m.set(0, 1, 7)
    assert_equal(m.get(0, 0), 10)
    assert_equal(m.get(1, 2), -5)
    assert_equal(m.get(0, 1), 7)


def test_matrix2d_add():
    var m = Matrix2D(2, 2)
    m.set(0, 0, 1)
    m.add(0, 0, 2)
    m.add(1, 1, 5)
    assert_equal(m.get(0, 0), 3)
    assert_equal(m.get(1, 1), 5)


def test_matrix2d_1x1():
    var m = Matrix2D(1, 1)
    assert_equal(m.get(0, 0), 0)
    m.set(0, 0, 99)
    assert_equal(m.get(0, 0), 99)
    m.add(0, 0, 1)
    assert_equal(m.get(0, 0), 100)


# ----- 2. get_exp_interval -----


def test_get_exp_interval_small():
    var bins = get_exp_interval(10)
    assert_true(len(bins) > 0)
    assert_equal(bins[0], 1)
    assert_equal(bins[-1], 10)
    for i in range(len(bins) - 1):
        assert_true(bins[i] < bins[i + 1])


def test_get_exp_interval_last_is_max_len():
    for max_len in [50, 100]:
        var bins = get_exp_interval(max_len)
        assert_equal(bins[-1], max_len)
        for i in range(len(bins) - 1):
            assert_true(bins[i] < bins[i + 1])


def test_get_exp_interval_large():
    var bins = get_exp_interval(500)
    assert_equal(bins[-1], 500)
    for i in range(len(bins) - 1):
        assert_true(bins[i] < bins[i + 1])


# ----- 3. get_linear_interval and grouping -----


def test_get_linear_interval():
    var i76 = get_linear_interval(76)
    var i100 = get_linear_interval(100)
    var i300 = get_linear_interval(300)
    assert_true(i76 > 0)
    assert_true(i100 > 0)
    assert_true(i300 > 0)


def test_make_ungrouped_groups():
    var g = make_ungrouped_groups(5)
    assert_equal(len(g), 5)
    assert_equal(g[0], 1)
    assert_equal(g[1], 2)
    assert_equal(g[4], 5)


def test_make_linear_base_groups_short():
    var max_length = 10
    var got = make_linear_base_groups(max_length)
    var expected = make_ungrouped_groups(max_length)
    assert_equal(len(got), len(expected))
    for i in range(len(got)):
        assert_equal(got[i], expected[i])


def test_make_linear_base_groups_long():
    var max_length = 200
    var groups = make_linear_base_groups(max_length)
    assert_true(len(groups) > 0)
    assert_equal(groups[0], 1)
    assert_equal(groups[-1], max_length)
    for i in range(len(groups) - 1):
        assert_true(groups[i] <= groups[i + 1])


# ----- 4. base2int -----


def test_base2int_acgt():
    # Uppercase A,C,G,T -> 0,1,2,3 (ASCII 65,67,71,84)
    assert_equal(base2int("A".as_bytes()[0]), 0)
    assert_equal(base2int("C".as_bytes()[0]), 1)
    assert_equal(base2int("G".as_bytes()[0]), 2)
    assert_equal(base2int("T".as_bytes()[0]), 3)


def test_base2int_invalid():
    assert_equal(base2int("N".as_bytes()[0]), 4)


# ----- 5. _seq_to_hash -----


def test_seq_to_hash_deterministic():
    var h1 = _seq_to_hash("ACGT")
    var h2 = _seq_to_hash("ACGT")
    assert_equal(h1, h2)


def test_seq_to_hash_different():
    var h1 = _seq_to_hash("A")
    var h2 = _seq_to_hash("C")
    assert_true(h1 != h2)


# ----- 6. grow_tensor -----


def test_grow_tensor_larger():
    var old_list = List[Int64]()
    old_list.append(1)
    old_list.append(2)
    old_list.append(3)
    var new_list = grow_tensor(old_list, 5)
    assert_equal(len(new_list), 5)
    assert_equal(new_list[0], 1)
    assert_equal(new_list[1], 2)
    assert_equal(new_list[2], 3)
    assert_equal(new_list[3], 0)
    assert_equal(new_list[4], 0)


def test_grow_tensor_same_size():
    var old_list = List[Int64]()
    old_list.append(1)
    old_list.append(2)
    var new_list = grow_tensor(old_list, 2)
    assert_equal(len(new_list), 2)
    assert_equal(new_list[0], 1)
    assert_equal(new_list[1], 2)


def test_grow_tensor_empty():
    var old_list = List[Int64]()
    var new_list = grow_tensor(old_list, 3)
    assert_equal(len(new_list), 3)
    assert_equal(new_list[0], 0)
    assert_equal(new_list[1], 0)
    assert_equal(new_list[2], 0)


# ----- 7. grow_matrix -----


def test_grow_matrix_larger():
    var old_mat = Matrix2D(2, 2)
    old_mat.set(0, 0, 1)
    old_mat.set(0, 1, 1)
    old_mat.set(1, 0, 1)
    old_mat.set(1, 1, 1)
    var new_mat = grow_matrix(old_mat, 3, 3)
    assert_equal(new_mat.shape()[0], 3)
    assert_equal(new_mat.shape()[1], 3)
    assert_equal(new_mat.get(0, 0), 1)
    assert_equal(new_mat.get(0, 1), 1)
    assert_equal(new_mat.get(1, 0), 1)
    assert_equal(new_mat.get(1, 1), 1)
    assert_equal(new_mat.get(0, 2), 0)
    assert_equal(new_mat.get(2, 2), 0)


def test_grow_matrix_same_size():
    var old_mat = Matrix2D(2, 2)
    old_mat.set(0, 0, 5)
    old_mat.set(1, 1, 6)
    var new_mat = grow_matrix(old_mat, 2, 2)
    assert_equal(new_mat.get(0, 0), 5)
    assert_equal(new_mat.get(1, 1), 6)


# ----- 8. sum_tensor -----


def test_sum_tensor_empty():
    var empty = List[Int64]()
    assert_equal(sum_tensor(empty), 0)


def test_sum_tensor_single():
    var one = List[Int64]()
    one.append(42)
    assert_equal(sum_tensor(one), 42)


def test_sum_tensor_multiple():
    var many = List[Int64]()
    many.append(1)
    many.append(2)
    many.append(3)
    assert_equal(sum_tensor(many), 6)


def main():
    TestSuite.discover_tests[__functions_in_module()]().run()
