"""Unit tests for blazeqc.stats.duplication and over_represented (pure Mojo)."""

from blazeseq import FastqRecord
from blazeqc.stats.duplication import DupReads
from blazeqc.stats.over_represented import OverRepresentedSequence
from testing import assert_equal, assert_true, TestSuite


# ----- 1. DupReads — initialisation -----


def test_dup_reads_init_n_zero():
    var dr = DupReads()
    assert_equal(dr.n, 0)


def test_dup_reads_init_unique_reads_zero():
    var dr = DupReads()
    assert_equal(dr.unique_reads, 0)


def test_dup_reads_init_count_at_max_zero():
    var dr = DupReads()
    assert_equal(dr.count_at_max, 0)


# ----- 2. DupReads — tally_read -----


def test_dup_reads_tally_n_increments():
    var dr = DupReads()
    var rec = FastqRecord("r1", "ACGT", "IIII")
    dr.tally_read(rec)
    assert_equal(dr.n, 1)


def test_dup_reads_tally_n_increments_on_duplicate():
    var dr = DupReads()
    var rec = FastqRecord("r1", "ACGT", "IIII")
    dr.tally_read(rec)
    dr.tally_read(rec)
    assert_equal(dr.n, 2)


def test_dup_reads_tally_unique_read_increments_unique_reads():
    var dr = DupReads()
    var rec = FastqRecord("r1", "ACGT", "IIII")
    dr.tally_read(rec)
    assert_equal(dr.unique_reads, 1)


def test_dup_reads_tally_duplicate_does_not_increment_unique_reads():
    # Same sequence twice → unique_reads stays 1
    var dr = DupReads()
    var rec = FastqRecord("r1", "ACGT", "IIII")
    dr.tally_read(rec)
    dr.tally_read(rec)
    assert_equal(dr.unique_reads, 1)


def test_dup_reads_tally_two_different_reads():
    var dr = DupReads()
    var rec1 = FastqRecord("r1", "ACGT", "IIII")
    var rec2 = FastqRecord("r2", "TTTT", "IIII")
    dr.tally_read(rec1)
    dr.tally_read(rec2)
    assert_equal(dr.unique_reads, 2)
    assert_equal(dr.n, 2)


def test_dup_reads_tally_count_at_max_updated():
    # count_at_max tracks n at the point the unique read was added
    var dr = DupReads()
    var rec = FastqRecord("r1", "ACGT", "IIII")
    dr.tally_read(rec)
    # count_at_max should equal n (1) after the first unique read
    assert_equal(dr.count_at_max, 1)


def test_dup_reads_tally_many_unique_reads():
    var dr = DupReads()
    for i in range(10):
        var rec = FastqRecord("r" + String(i), "ACGT" + String(i), "IIIII")
        dr.tally_read(rec)
    assert_equal(dr.unique_reads, 10)
    assert_equal(dr.n, 10)


# ----- 3. DupReads.correct_values (static method) -----


def test_correct_values_count_at_max_equals_total():
    # When count_at_max == total_count → return count_at_level unchanged
    var result = DupReads.correct_values(2, 5, 100, 100)
    assert_equal(result, 5.0)


def test_correct_values_remaining_less_than_count_at_max():
    # When total_count - count_at_level < count_at_max → return count_at_level
    # 50 - 40 = 10 < 20 → return 40
    var result = DupReads.correct_values(2, 40, 20, 50)
    assert_equal(result, 40.0)


def test_correct_values_returns_float64():
    # Just verify it runs and returns a non-negative value in the normal path
    var result = DupReads.correct_values(1, 1, 50, 1000)
    assert_true(result > 0.0)


# ----- 4. OverRepresentedSequence -----


def test_over_repr_init_string():
    var ors = OverRepresentedSequence(
        String("ACGT"), 42, 1.5, String("No Hit")
    )
    assert_equal(ors.seq, String("ACGT"))
    assert_equal(ors.count, 42)
    assert_equal(ors.percentage, 1.5)
    assert_equal(ors.hit, String("No Hit"))


def test_over_repr_init_string_literal():
    var ors = OverRepresentedSequence("ACGT", 10, 0.5, "No Hit")
    assert_equal(ors.seq, String("ACGT"))
    assert_equal(ors.count, 10)
    assert_equal(ors.percentage, 0.5)
    assert_equal(ors.hit, String("No Hit"))


def test_over_repr_zero_percentage():
    var ors = OverRepresentedSequence("AAA", 0, 0.0, "No Hit")
    assert_equal(ors.count, 0)
    assert_equal(ors.percentage, 0.0)


# ----- _get_status_duplication (pass/warn/fail) -----
# DUPLICATION_WARN=70, DUPLICATION_ERROR=50. Percent remaining after dedup.


def test_dup_status_pass():
    var dr = DupReads()
    for i in range(20):
        var rec = FastqRecord("r" + String(i), "ACGTACGTACGT" + String(i), "IIIIIIIIIIII")
        dr.tally_read(rec)
    assert_equal(dr._get_status_duplication(20), "pass")


def test_dup_status_warn():
    var dr = DupReads()
    for i in range(20):
        var rec = FastqRecord("r" + String(i), "ACGTACGTACGT" + String(i), "IIIIIIIIIIII")
        dr.tally_read(rec)
    for i in range(20, 60):
        for _ in range(2):
            var rec = FastqRecord("r" + String(i), "ACGTACGTACGT" + String(i), "IIIIIIIIIIII")
            dr.tally_read(rec)
    assert_equal(dr._get_status_duplication(100), "warn")


def test_dup_status_fail():
    var dr = DupReads()
    for i in range(50):
        for _ in range(4):
            var rec = FastqRecord("r" + String(i), "ACGTACGTACGT" + String(i), "IIIIIIIIIIII")
            dr.tally_read(rec)
    assert_equal(dr._get_status_duplication(200), "fail")


def main():
    TestSuite.discover_tests[__functions_in_module()]().run()
