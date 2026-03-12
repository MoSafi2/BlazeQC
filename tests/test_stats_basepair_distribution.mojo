"""Unit tests for blazeqc.stats.basepair_distribution (pure Mojo)."""

from blazeseq import FastqRecord
from blazeqc.stats.basepair_distribution import (
    BasepairCollector,
    BasepairModule,
)
from blazeqc.stats.summary_utils import SummaryContext
from testing import assert_equal, assert_true, TestSuite


# ----- 1. Initialisation (collector) -----


def test_bp_dist_init_max_length_zero():
    var c = BasepairCollector()
    assert_equal(c.max_length, 0)


def test_bp_dist_init_min_length_max_int():
    var c = BasepairCollector()
    assert_equal(c.min_length, Int.MAX)


def test_bp_dist_init_matrix_shape():
    var c = BasepairCollector()
    var sh = c.bp_dist.shape()
    assert_equal(sh[0], 1)
    assert_equal(sh[1], 5)


def test_bp_dist_init_matrix_zeros():
    var c = BasepairCollector()
    assert_equal(c.bp_dist.get(0, 0), 0)
    assert_equal(c.bp_dist.get(0, 4), 0)


# ----- 2. tally_read — length tracking (collector) -----


def test_bp_dist_tally_updates_max_length():
    var c = BasepairCollector()
    var rec = FastqRecord("r1", "ACGT", "IIII")
    c.tally_read(rec)
    assert_equal(c.max_length, 4)


def test_bp_dist_tally_updates_min_length():
    var c = BasepairCollector()
    var rec = FastqRecord("r1", "ACGT", "IIII")
    c.tally_read(rec)
    assert_equal(c.min_length, 4)


def test_bp_dist_tally_max_length_grows():
    var c = BasepairCollector()
    var rec2 = FastqRecord("r1", "AC", "II")
    var rec6 = FastqRecord("r2", "ACGTAC", "IIIIII")
    c.tally_read(rec2)
    assert_equal(c.max_length, 2)
    c.tally_read(rec6)
    assert_equal(c.max_length, 6)


def test_bp_dist_tally_min_length_tracks_shorter():
    var c = BasepairCollector()
    var rec6 = FastqRecord("r1", "ACGTAC", "IIIIII")
    var rec2 = FastqRecord("r2", "AC", "II")
    c.tally_read(rec6)
    assert_equal(c.min_length, 6)
    c.tally_read(rec2)
    assert_equal(c.min_length, 2)


# ----- 3. tally_read — base-to-column mapping (collector) -----
# Column = (ASCII & 0b11111) % 5:
#   'T' (84) → & 31 = 20 → % 5 = 0
#   'A' (65) → & 31 =  1 → % 5 = 1
#   'G' (71) → & 31 =  7 → % 5 = 2
#   'C' (67) → & 31 =  3 → % 5 = 3
#   'N' (78) → & 31 = 14 → % 5 = 4


def test_bp_dist_tally_base_T_maps_to_col0():
    var c = BasepairCollector()
    var rec = FastqRecord("r1", "T", "I")
    c.tally_read(rec)
    assert_equal(c.bp_dist.get(0, 0), 1)   # col 0 (T)
    assert_equal(c.bp_dist.get(0, 1), 0)
    assert_equal(c.bp_dist.get(0, 2), 0)
    assert_equal(c.bp_dist.get(0, 3), 0)
    assert_equal(c.bp_dist.get(0, 4), 0)


def test_bp_dist_tally_base_A_maps_to_col1():
    var c = BasepairCollector()
    var rec = FastqRecord("r1", "A", "I")
    c.tally_read(rec)
    assert_equal(c.bp_dist.get(0, 1), 1)   # col 1 (A)


def test_bp_dist_tally_base_G_maps_to_col2():
    var c = BasepairCollector()
    var rec = FastqRecord("r1", "G", "I")
    c.tally_read(rec)
    assert_equal(c.bp_dist.get(0, 2), 1)   # col 2 (G)


def test_bp_dist_tally_base_C_maps_to_col3():
    var c = BasepairCollector()
    var rec = FastqRecord("r1", "C", "I")
    c.tally_read(rec)
    assert_equal(c.bp_dist.get(0, 3), 1)   # col 3 (C)


def test_bp_dist_tally_base_N_maps_to_col4():
    var c = BasepairCollector()
    var rec = FastqRecord("r1", "N", "I")
    c.tally_read(rec)
    assert_equal(c.bp_dist.get(0, 4), 1)   # col 4 (N)


def test_bp_dist_tally_multi_position_mapping():
    var c = BasepairCollector()
    var rec = FastqRecord("r1", "TAGCN", "IIIII")
    c.tally_read(rec)
    assert_equal(c.bp_dist.get(0, 0), 1)
    assert_equal(c.bp_dist.get(1, 1), 1)
    assert_equal(c.bp_dist.get(2, 2), 1)
    assert_equal(c.bp_dist.get(3, 3), 1)
    assert_equal(c.bp_dist.get(4, 4), 1)


def test_bp_dist_tally_accumulates_same_position():
    var c = BasepairCollector()
    var rec = FastqRecord("r1", "A", "I")
    c.tally_read(rec)
    c.tally_read(rec)
    assert_equal(c.bp_dist.get(0, 1), 2)


# ----- Status (N content, sequence content) via module + summarizers -----
# N: N_CONTENT_WARN=5, N_CONTENT_ERROR=20. Sequence: SEQUENCE_WARN=10, SEQUENCE_ERROR=20.


def test_bp_dist_status_n_pass():
    var mod = BasepairModule()
    for _ in range(24):
        mod.collector.tally_read(FastqRecord("r", "A", "I"))
    mod.collector.tally_read(FastqRecord("r", "N", "I"))
    var ctx = SummaryContext(26, 26, "test")
    mod.prepare_summarizers(ctx)
    assert_equal(mod.summarizer_n.grade().grade, "pass")


def test_bp_dist_status_n_fail():
    var mod = BasepairModule()
    for _ in range(3):
        mod.collector.tally_read(FastqRecord("r", "A", "I"))
    for _ in range(2):
        mod.collector.tally_read(FastqRecord("r", "N", "I"))
    var ctx = SummaryContext(5, 5, "test")
    mod.prepare_summarizers(ctx)
    assert_equal(mod.summarizer_n.grade().grade, "fail")


def test_bp_dist_status_sequence_fail():
    var mod = BasepairModule()
    for _ in range(6):
        mod.collector.tally_read(FastqRecord("r", "A", "I"))
    for _ in range(2):
        mod.collector.tally_read(FastqRecord("r", "T", "I"))
    for _ in range(2):
        mod.collector.tally_read(FastqRecord("r", "G", "I"))
    var ctx = SummaryContext(10, 10, "test")
    mod.prepare_summarizers(ctx)
    assert_equal(mod.summarizer_seq.grade().grade, "fail")


def main():
    TestSuite.discover_tests[__functions_in_module()]().run()
