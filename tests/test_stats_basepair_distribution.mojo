"""Unit tests for blazeqc.stats.basepair_distribution (pure Mojo)."""

from blazeseq import FastqRecord
from blazeqc.stats.basepair_distribution import BasepairDistribution
from testing import assert_equal, assert_true, TestSuite


# ----- 1. Initialisation -----


def test_bp_dist_init_max_length_zero():
    var bd = BasepairDistribution()
    assert_equal(bd.max_length, 0)


def test_bp_dist_init_min_length_max_int():
    var bd = BasepairDistribution()
    assert_equal(bd.min_length, Int.MAX)


def test_bp_dist_init_matrix_shape():
    # bp_dist is initialized as Matrix2D(1, WIDTH=5)
    var bd = BasepairDistribution()
    var sh = bd.bp_dist.shape()
    assert_equal(sh[0], 1)
    assert_equal(sh[1], 5)


def test_bp_dist_init_matrix_zeros():
    var bd = BasepairDistribution()
    assert_equal(bd.bp_dist.get(0, 0), 0)
    assert_equal(bd.bp_dist.get(0, 4), 0)


# ----- 2. tally_read — length tracking -----


def test_bp_dist_tally_updates_max_length():
    var bd = BasepairDistribution()
    var rec = FastqRecord("r1", "ACGT", "IIII")
    bd.tally_read(rec)
    assert_equal(bd.max_length, 4)


def test_bp_dist_tally_updates_min_length():
    var bd = BasepairDistribution()
    var rec = FastqRecord("r1", "ACGT", "IIII")
    bd.tally_read(rec)
    assert_equal(bd.min_length, 4)


def test_bp_dist_tally_max_length_grows():
    var bd = BasepairDistribution()
    var rec2 = FastqRecord("r1", "AC", "II")
    var rec6 = FastqRecord("r2", "ACGTAC", "IIIIII")
    bd.tally_read(rec2)
    assert_equal(bd.max_length, 2)
    bd.tally_read(rec6)
    assert_equal(bd.max_length, 6)


def test_bp_dist_tally_min_length_tracks_shorter():
    var bd = BasepairDistribution()
    var rec6 = FastqRecord("r1", "ACGTAC", "IIIIII")
    var rec2 = FastqRecord("r2", "AC", "II")
    bd.tally_read(rec6)
    assert_equal(bd.min_length, 6)
    bd.tally_read(rec2)
    assert_equal(bd.min_length, 2)


# ----- 3. tally_read — base-to-column mapping -----
# Column = (ASCII & 0b11111) % 5:
#   'T' (84) → & 31 = 20 → % 5 = 0
#   'A' (65) → & 31 =  1 → % 5 = 1
#   'G' (71) → & 31 =  7 → % 5 = 2
#   'C' (67) → & 31 =  3 → % 5 = 3
#   'N' (78) → & 31 = 14 → % 5 = 4


def test_bp_dist_tally_base_T_maps_to_col0():
    var bd = BasepairDistribution()
    var rec = FastqRecord("r1", "T", "I")
    bd.tally_read(rec)
    assert_equal(bd.bp_dist.get(0, 0), 1)   # col 0 (T)
    assert_equal(bd.bp_dist.get(0, 1), 0)
    assert_equal(bd.bp_dist.get(0, 2), 0)
    assert_equal(bd.bp_dist.get(0, 3), 0)
    assert_equal(bd.bp_dist.get(0, 4), 0)


def test_bp_dist_tally_base_A_maps_to_col1():
    var bd = BasepairDistribution()
    var rec = FastqRecord("r1", "A", "I")
    bd.tally_read(rec)
    assert_equal(bd.bp_dist.get(0, 1), 1)   # col 1 (A)


def test_bp_dist_tally_base_G_maps_to_col2():
    var bd = BasepairDistribution()
    var rec = FastqRecord("r1", "G", "I")
    bd.tally_read(rec)
    assert_equal(bd.bp_dist.get(0, 2), 1)   # col 2 (G)


def test_bp_dist_tally_base_C_maps_to_col3():
    var bd = BasepairDistribution()
    var rec = FastqRecord("r1", "C", "I")
    bd.tally_read(rec)
    assert_equal(bd.bp_dist.get(0, 3), 1)   # col 3 (C)


def test_bp_dist_tally_base_N_maps_to_col4():
    var bd = BasepairDistribution()
    var rec = FastqRecord("r1", "N", "I")
    bd.tally_read(rec)
    assert_equal(bd.bp_dist.get(0, 4), 1)   # col 4 (N)


def test_bp_dist_tally_multi_position_mapping():
    # "TAGCN": T→col0, A→col1, G→col2, C→col3, N→col4 — each at a different position
    var bd = BasepairDistribution()
    var rec = FastqRecord("r1", "TAGCN", "IIIII")
    bd.tally_read(rec)
    assert_equal(bd.bp_dist.get(0, 0), 1)   # pos 0: T → col 0
    assert_equal(bd.bp_dist.get(1, 1), 1)   # pos 1: A → col 1
    assert_equal(bd.bp_dist.get(2, 2), 1)   # pos 2: G → col 2
    assert_equal(bd.bp_dist.get(3, 3), 1)   # pos 3: C → col 3
    assert_equal(bd.bp_dist.get(4, 4), 1)   # pos 4: N → col 4


def test_bp_dist_tally_accumulates_same_position():
    # Two reads, same position, same base → count doubles
    var bd = BasepairDistribution()
    var rec = FastqRecord("r1", "A", "I")
    bd.tally_read(rec)
    bd.tally_read(rec)
    assert_equal(bd.bp_dist.get(0, 1), 2)   # col 1 (A), two counts


def main():
    TestSuite.discover_tests[__functions_in_module()]().run()
