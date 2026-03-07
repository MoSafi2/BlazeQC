"""Unit tests for blazeqc.stats.tile_quality (pure Mojo)."""

from blazeseq import FastqRecord
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


def main():
    TestSuite.discover_tests[__functions_in_module()]().run()
