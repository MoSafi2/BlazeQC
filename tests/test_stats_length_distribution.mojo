"""Unit tests for blazeqc.stats.length_distribution (pure Mojo)."""

from blazeseq import FastqRecord
from blazeqc.stats.length_distribution import LengthDistribution
from testing import assert_equal, assert_true, TestSuite


# ----- 1. Initialisation -----


def test_length_dist_init_empty():
    var ld = LengthDistribution()
    assert_equal(len(ld.length_vector), 0)


# ----- 2. tally_read -----
# length_vector is 0-indexed: a read of length N increments length_vector[N-1]


def test_length_dist_tally_single_read():
    var ld = LengthDistribution()
    var rec = FastqRecord("r1", "ACGT", "IIII")  # length 4
    ld.tally_read(rec)
    assert_equal(len(ld.length_vector), 4)
    assert_equal(ld.length_vector[3], 1)  # index N-1 = 3


def test_length_dist_tally_slot_zero():
    # A single-base read should increment slot 0 (length 1 → index 0)
    var ld = LengthDistribution()
    var rec = FastqRecord("r1", "A", "I")
    ld.tally_read(rec)
    assert_equal(ld.length_vector[0], 1)


def test_length_dist_tally_multiple_lengths():
    var ld = LengthDistribution()
    var rec2 = FastqRecord("r1", "AC", "II")        # length 2 → slot 1
    var rec4 = FastqRecord("r2", "ACGT", "IIII")   # length 4 → slot 3
    ld.tally_read(rec2)
    ld.tally_read(rec4)
    assert_equal(ld.length_vector[1], 1)
    assert_equal(ld.length_vector[3], 1)
    assert_equal(ld.length_vector[0], 0)
    assert_equal(ld.length_vector[2], 0)


def test_length_dist_tally_same_length_accumulates():
    var ld = LengthDistribution()
    var rec = FastqRecord("r1", "ACGT", "IIII")
    ld.tally_read(rec)
    ld.tally_read(rec)
    ld.tally_read(rec)
    assert_equal(ld.length_vector[3], 3)


def test_length_dist_tally_shorter_then_longer():
    # Tally a short read first, then a longer one — vector must grow
    var ld = LengthDistribution()
    var short_rec = FastqRecord("r1", "AC", "II")
    var long_rec = FastqRecord("r2", "ACGTACGT", "IIIIIIII")
    ld.tally_read(short_rec)
    ld.tally_read(long_rec)
    assert_equal(len(ld.length_vector), 8)
    assert_equal(ld.length_vector[1], 1)   # length 2 → slot 1
    assert_equal(ld.length_vector[7], 1)   # length 8 → slot 7


# ----- 3. length_average -----
# cum = sum(length_vector[i] * (i+1)); average = cum / num_reads


def test_length_average_uniform():
    var ld = LengthDistribution()
    var rec = FastqRecord("r1", "ACGT", "IIII")   # length 4
    ld.tally_read(rec)
    ld.tally_read(rec)
    # cum = 4*2 = 8, num_reads = 2, average = 4.0
    var avg = ld.length_average(2)
    assert_equal(avg, 4.0)


def test_length_average_mixed():
    var ld = LengthDistribution()
    var rec2 = FastqRecord("r1", "AC", "II")       # length 2
    var rec4 = FastqRecord("r2", "ACGT", "IIII")  # length 4
    ld.tally_read(rec2)
    ld.tally_read(rec4)
    # cum = 2 + 4 = 6, num_reads = 2, average = 3.0
    var avg = ld.length_average(2)
    assert_equal(avg, 3.0)


def test_length_average_single_read():
    var ld = LengthDistribution()
    var rec = FastqRecord("r1", "ACGTACGT", "IIIIIIII")  # length 8
    ld.tally_read(rec)
    var avg = ld.length_average(1)
    assert_equal(avg, 8.0)


# ----- 4. get_size_distribution -----
# Returns (starting, interval) such that (max_val - min_val) / interval <= 50


def test_get_size_distribution_small_range():
    var ld = LengthDistribution()
    starting, interval = ld.get_size_distribution(0, 10)
    # range=10, interval should be 1 (10/1=10 ≤ 50)
    assert_equal(interval, 1)
    assert_equal(starting, 0)


def test_get_size_distribution_medium_range():
    var ld = LengthDistribution()
    starting, interval = ld.get_size_distribution(0, 51)
    # range=51, interval 1 gives 51 > 50, interval 2 gives 25.5 ≤ 50
    assert_equal(interval, 2)
    assert_equal(starting, 0)


def test_get_size_distribution_nonzero_start():
    var ld = LengthDistribution()
    starting, interval = ld.get_size_distribution(100, 200)
    # range=100, interval 2 gives 50 ≤ 50
    assert_equal(interval, 2)
    assert_equal(starting, 100)


def test_get_size_distribution_interval_nonnegative():
    var ld = LengthDistribution()
    starting, interval = ld.get_size_distribution(0, 1000)
    assert_true(interval > 0)
    # range/interval should be ≤ 50
    assert_true((1000 - 0) / interval <= 50)


def main():
    TestSuite.discover_tests[__functions_in_module()]().run()
