"""Unit tests for blazeqc.stats.tile_quality (pure Mojo)."""

from blazeseq import FastqRecord
from blazeqc.helpers import Matrix2D
from blazeqc.stats.tile_quality import PerTileQuality, TileQualityEntry
from testing import assert_equal, assert_true, TestSuite


# ----- 1. Initialisation -----


def test_per_tile_init_n_zero():
    var ptq = PerTileQuality()
    assert_equal(ptq.n, 0)


def test_per_tile_init_max_length_zero():
    var ptq = PerTileQuality()
    assert_equal(ptq.max_length, 0)


# ----- 2. TileQualityEntry -----


def test_tile_entry_init():
    var entry = TileQualityEntry(42, 1, 4)
    assert_equal(entry.tile, 42)
    assert_equal(entry.count, 1)
    assert_equal(len(entry.quality), 4)


def test_tile_entry_quality_zeros():
    var entry = TileQualityEntry(1, 0, 3)
    for i in range(3):
        assert_equal(entry.quality[i], 0)


def test_tile_entry_iadd():
    var entry = TileQualityEntry(1, 5, 0)
    entry += 3
    assert_equal(entry.count, 8)


def test_tile_entry_add():
    var entry = TileQualityEntry(1, 5, 0)
    var result = entry + 3
    assert_equal(result, 8)


# ----- 3. _find_tile_info — colon counting -----
# < 4 colons → -1
# >= 4 colons → 2
# >= 6 colons → 4


def test_find_tile_info_no_colons():
    # "readname" has 0 colons → -1
    var ptq = PerTileQuality()
    var rec = FastqRecord("readname", "ACGT", "IIII")
    var result = ptq._find_tile_info(rec)
    assert_equal(result, -1)


def test_find_tile_info_three_colons():
    # "SIM:1:3:15" has 3 colons → -1
    var ptq = PerTileQuality()
    var rec = FastqRecord("SIM:1:3:15", "ACGT", "IIII")
    var result = ptq._find_tile_info(rec)
    assert_equal(result, -1)


def test_find_tile_info_four_colons():
    # "SIM:1:3:15:42" has 4 colons → 2
    var ptq = PerTileQuality()
    var rec = FastqRecord("SIM:1:3:15:42", "ACGT", "IIII")
    var result = ptq._find_tile_info(rec)
    assert_equal(result, 2)


def test_find_tile_info_five_colons():
    # "SIM:1:3:15:42:99" has 5 colons → 2 (>= 4 but < 6)
    var ptq = PerTileQuality()
    var rec = FastqRecord("SIM:1:3:15:42:99", "ACGT", "IIII")
    var result = ptq._find_tile_info(rec)
    assert_equal(result, 2)


def test_find_tile_info_six_colons():
    # "SIM:1:FCX:1:15:6329:1045" has 6 colons → 4
    var ptq = PerTileQuality()
    var rec = FastqRecord("SIM:1:FCX:1:15:6329:1045", "ACGT", "IIII")
    var result = ptq._find_tile_info(rec)
    assert_equal(result, 4)


# ----- 4. _find_tile_value — numeric field extraction -----


def test_find_tile_value_four_colon_id():
    # ID = "SIM:1:3:15:42", pos=2 → extract field between colon #2 and colon #3
    # "SIM:1:3:15:42"
    #      ^ ^    → colons at idx 3, 5, 7, 10 (count 1,2,3,4)
    # pos=2 → index_1 = colon2_pos+1 = 6, index_2 = colon3_pos = 7 → s[6:7] = "3"
    var ptq = PerTileQuality()
    var rec = FastqRecord("SIM:1:3:15:42", "ACGT", "IIII")
    var result = ptq._find_tile_value(rec, 2)
    assert_equal(result, 3)


def test_find_tile_value_six_colon_id():
    # ID = "SIM:1:FCX:1:15:6329:1045", pos=4
    # Colon positions: 3,5,9,11,14,19 (counts 1,2,3,4,5,6)
    # pos=4 → index_1 = 12, index_2 = 14 → s[12:14] = "15"
    var ptq = PerTileQuality()
    var rec = FastqRecord("SIM:1:FCX:1:15:6329:1045", "ACGT", "IIII")
    var result = ptq._find_tile_value(rec, 4)
    assert_equal(result, 15)


def test_find_tile_value_returns_zero_on_non_numeric():
    # If the field between colons is not numeric, atol raises and returns 0
    var ptq = PerTileQuality()
    # "A:B:C:D:E" → pos=2 → field between colon #2 and #3 = "C" → atol raises → 0
    var rec = FastqRecord("A:B:C:D:E", "ACGT", "IIII")
    var result = ptq._find_tile_value(rec, 2)
    assert_equal(result, 0)


# ----- 5. tally_read — basic state tracking -----


def test_per_tile_tally_increments_n():
    var ptq = PerTileQuality()
    var rec = FastqRecord("SIM:1:FCX:1:15:6329:1045", "ACGT", "IIII")
    ptq.tally_read(rec)
    assert_equal(ptq.n, 1)


def test_per_tile_tally_updates_max_length():
    var ptq = PerTileQuality()
    var rec = FastqRecord("SIM:1:FCX:1:15:6329:1045", "ACGT", "IIII")
    ptq.tally_read(rec)
    assert_equal(ptq.max_length, 4)


def test_per_tile_tally_n_increments_on_no_tile():
    # A record with no tile info (-1) still increments n
    var ptq = PerTileQuality()
    var rec = FastqRecord("readname", "ACGT", "IIII")
    ptq.tally_read(rec)
    assert_equal(ptq.n, 1)


# ----- 6. FastQC alignment: normalization and plot -----


def test_tile_quality_normalization_column_sums_zero():
    """After _subtract_group_averages, each column (group) sums to zero across tiles (FastQC normalization)."""
    var ptq = PerTileQuality()
    var means = Matrix2D[DType.float64](2, 3)
    means.set(0, 0, 10.0)
    means.set(0, 1, 20.0)
    means.set(0, 2, 30.0)
    means.set(1, 0, 14.0)
    means.set(1, 1, 22.0)
    means.set(1, 2, 26.0)
    _ = ptq._subtract_group_averages(means, 2, 3)
    assert_equal(means.col_sum(0), Float64(0.0))
    assert_equal(means.col_sum(1), Float64(0.0))
    assert_equal(means.col_sum(2), Float64(0.0))


def test_tile_quality_plot_runs_and_sets_max_deviation():
    """Plot runs on minimal tile data and sets max_deviation (sanity check for FastQC-aligned pipeline)."""
    var ptq = PerTileQuality()
    var rec1 = FastqRecord("SIM:1:FCX:1:15:6329:1045", "ACGT", "IIII")
    var rec2 = FastqRecord("SIM:1:FCX:2:15:6329:1046", "ACGT", "IIII")
    ptq.tally_read(rec1)
    ptq.tally_read(rec2)
    _ = ptq.plot()
    assert_true(ptq.max_deviation >= 0.0)
    assert_true(ptq.max_length > 0)


# ----- 6. _get_status (pass/warn/fail) -----
# Limits: TILE_WARN=5, TILE_ERROR=10


def test_tile_quality_status_pass():
    var ptq = PerTileQuality()
    ptq.max_deviation = 2.0
    assert_equal(ptq._get_status(), "pass")


def test_tile_quality_status_warn():
    var ptq = PerTileQuality()
    ptq.max_deviation = 6.0
    assert_equal(ptq._get_status(), "warn")


def test_tile_quality_status_fail():
    var ptq = PerTileQuality()
    ptq.max_deviation = 12.0
    assert_equal(ptq._get_status(), "fail")


def main():
    TestSuite.discover_tests[__functions_in_module()]().run()
