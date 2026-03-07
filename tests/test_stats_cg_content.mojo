"""Unit tests for blazeqc.stats.cg_content (pure Mojo, no Python)."""

from blazeseq import FastqRecord
from blazeqc.stats.cg_content import CGContent
from testing import assert_equal, assert_true, TestSuite


# ----- 1. Initialisation -----


def test_cg_content_init_list_length():
    var cg = CGContent()
    assert_equal(len(cg.cg_content), 101)


def test_cg_content_init_all_zeros():
    var cg = CGContent()
    for i in range(101):
        assert_equal(cg.cg_content[i], 0)


def test_cg_content_init_theoretical_zeros():
    var cg = CGContent()
    assert_equal(len(cg.theoritical_distribution), 101)
    for i in range(101):
        assert_equal(cg.theoritical_distribution[i], 0)


# ----- 2. tally_read — GC detection logic -----
# Bitmask used: sequence[i] & 0b111 == 3  (C=67→3) or == 7 (G=71→7)
# A=65→1, T=84→4, N=78→6  — none of these are 3 or 7


def test_cg_content_tally_100_percent_gc():
    # "CCGG": all four bases are C or G → GC% = 100
    var cg = CGContent()
    var rec = FastqRecord("r1", "CCGG", "IIII")
    cg.tally_read(rec)
    assert_equal(cg.cg_content[100], 1)
    assert_equal(cg.cg_content[0], 0)


def test_cg_content_tally_0_percent_gc():
    # "AATT": no C or G → GC% = 0
    var cg = CGContent()
    var rec = FastqRecord("r1", "AATT", "IIII")
    cg.tally_read(rec)
    assert_equal(cg.cg_content[0], 1)
    assert_equal(cg.cg_content[100], 0)


def test_cg_content_tally_50_percent_gc():
    # "ACGT": C and G out of 4 → 50%
    var cg = CGContent()
    var rec = FastqRecord("r1", "ACGT", "IIII")
    cg.tally_read(rec)
    assert_equal(cg.cg_content[50], 1)


def test_cg_content_tally_25_percent_gc():
    # "ACAA": one G/C out of 4 → 25%
    var cg = CGContent()
    var rec = FastqRecord("r1", "ACAA", "IIII")
    cg.tally_read(rec)
    assert_equal(cg.cg_content[25], 1)


def test_cg_content_tally_empty_record_no_change():
    # Empty sequence → early return, cg_content unchanged
    var cg = CGContent()
    var rec = FastqRecord("r1", "", "")
    cg.tally_read(rec)
    for i in range(101):
        assert_equal(cg.cg_content[i], 0)


def test_cg_content_tally_accumulates_across_reads():
    var cg = CGContent()
    var rec1 = FastqRecord("r1", "CCGG", "IIII")  # 100%
    var rec2 = FastqRecord("r2", "CCGG", "IIII")  # 100%
    var rec3 = FastqRecord("r3", "AATT", "IIII")  # 0%
    cg.tally_read(rec1)
    cg.tally_read(rec2)
    cg.tally_read(rec3)
    assert_equal(cg.cg_content[100], 2)
    assert_equal(cg.cg_content[0], 1)


def test_cg_content_tally_single_base_c():
    # Single 'C' → 100% GC
    var cg = CGContent()
    var rec = FastqRecord("r1", "C", "I")
    cg.tally_read(rec)
    assert_equal(cg.cg_content[100], 1)


def test_cg_content_tally_single_base_a():
    # Single 'A' → 0% GC
    var cg = CGContent()
    var rec = FastqRecord("r1", "A", "I")
    cg.tally_read(rec)
    assert_equal(cg.cg_content[0], 1)


def main():
    TestSuite.discover_tests[__functions_in_module()]().run()
