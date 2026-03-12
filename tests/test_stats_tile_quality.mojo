"""Unit tests for blazeqc.stats.tile_quality (pure Mojo)."""

from blazeseq import FastqRecord
from blazeqc.helpers import Matrix2D
from blazeqc.stats.tile_quality import (
    TileQualityCollector,
    TileQualityEntry,
    TileQualityModule,
    _subtract_group_averages,
    tile_quality_grade_from_deviation,
)
from blazeqc.stats.summary_utils import SummaryContext
from testing import assert_equal, assert_true, TestSuite


# ----- 1. Initialisation (collector) -----


def test_per_tile_init_n_zero():
    var c = TileQualityCollector()
    assert_equal(c.n, 0)


def test_per_tile_init_max_length_zero():
    var c = TileQualityCollector()
    assert_equal(c.max_length, 0)


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
    var c = TileQualityCollector()
    var rec = FastqRecord("readname", "ACGT", "IIII")
    var result = c._find_tile_info(rec)
    assert_equal(result, -1)


def test_find_tile_info_three_colons():
    var c = TileQualityCollector()
    var rec = FastqRecord("SIM:1:3:15", "ACGT", "IIII")
    var result = c._find_tile_info(rec)
    assert_equal(result, -1)


def test_find_tile_info_four_colons():
    var c = TileQualityCollector()
    var rec = FastqRecord("SIM:1:3:15:42", "ACGT", "IIII")
    var result = c._find_tile_info(rec)
    assert_equal(result, 2)


def test_find_tile_info_five_colons():
    var c = TileQualityCollector()
    var rec = FastqRecord("SIM:1:3:15:42:99", "ACGT", "IIII")
    var result = c._find_tile_info(rec)
    assert_equal(result, 2)


def test_find_tile_info_six_colons():
    var c = TileQualityCollector()
    var rec = FastqRecord("SIM:1:FCX:1:15:6329:1045", "ACGT", "IIII")
    var result = c._find_tile_info(rec)
    assert_equal(result, 4)


# ----- 4. _find_tile_value — numeric field extraction -----


def test_find_tile_value_four_colon_id():
    var c = TileQualityCollector()
    var rec = FastqRecord("SIM:1:3:15:42", "ACGT", "IIII")
    var result = c._find_tile_value(rec, 2)
    assert_equal(result, 3)


def test_find_tile_value_six_colon_id():
    var c = TileQualityCollector()
    var rec = FastqRecord("SIM:1:FCX:1:15:6329:1045", "ACGT", "IIII")
    var result = c._find_tile_value(rec, 4)
    assert_equal(result, 15)


def test_find_tile_value_returns_zero_on_non_numeric():
    var c = TileQualityCollector()
    var rec = FastqRecord("A:B:C:D:E", "ACGT", "IIII")
    var result = c._find_tile_value(rec, 2)
    assert_equal(result, 0)


# ----- 5. tally_read — basic state tracking -----


def test_per_tile_tally_increments_n():
    var c = TileQualityCollector()
    var rec = FastqRecord("SIM:1:FCX:1:15:6329:1045", "ACGT", "IIII")
    c.tally_read(rec)
    assert_equal(c.n, 1)


def test_per_tile_tally_updates_max_length():
    var c = TileQualityCollector()
    var rec = FastqRecord("SIM:1:FCX:1:15:6329:1045", "ACGT", "IIII")
    c.tally_read(rec)
    assert_equal(c.max_length, 4)


def test_per_tile_tally_n_increments_on_no_tile():
    var c = TileQualityCollector()
    var rec = FastqRecord("readname", "ACGT", "IIII")
    c.tally_read(rec)
    assert_equal(c.n, 1)


# ----- 6. FastQC alignment: normalization and plot -----


def test_tile_quality_normalization_column_sums_zero():
    var means = Matrix2D[DType.float64](2, 3)
    means.set(0, 0, 10.0)
    means.set(0, 1, 20.0)
    means.set(0, 2, 30.0)
    means.set(1, 0, 14.0)
    means.set(1, 1, 22.0)
    means.set(1, 2, 26.0)
    _ = _subtract_group_averages(means, 2, 3)
    assert_equal(means.col_sum(0), Float64(0.0))
    assert_equal(means.col_sum(1), Float64(0.0))
    assert_equal(means.col_sum(2), Float64(0.0))


def test_tile_quality_plot_runs_and_sets_max_deviation():
    var mod = TileQualityModule()
    var rec1 = FastqRecord("SIM:1:FCX:1:15:6329:1045", "ACGT", "IIII")
    var rec2 = FastqRecord("SIM:1:FCX:2:15:6329:1046", "ACGT", "IIII")
    mod.tally_read(rec1)
    mod.tally_read(rec2)
    var ctx = SummaryContext(2, 8, "test")
    mod.prepare_summarizers(ctx)
    _ = mod.plot_result()
    assert_true(mod.summarizer._prepared.max_deviation >= 0.0)
    assert_true(mod.collector.max_length > 0)


# ----- 6. Grade from deviation (pass/warn/fail) -----
# Limits: TILE_WARN=5, TILE_ERROR=10


def test_tile_quality_status_pass():
    assert_equal(tile_quality_grade_from_deviation(2.0), "pass")


def test_tile_quality_status_warn():
    assert_equal(tile_quality_grade_from_deviation(6.0), "warn")


def test_tile_quality_status_fail():
    assert_equal(tile_quality_grade_from_deviation(12.0), "fail")


def main():
    TestSuite.discover_tests[__functions_in_module()]().run()
