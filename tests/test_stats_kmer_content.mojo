"""Unit tests for blazeqc.stats.kmer_content (pure Mojo).

NOTE: kmer_to_index() is NOT tested directly because it iterates
      range(len(kmer), -1, -1) which accesses kmer[len(kmer)] — an
      out-of-bounds read that panics under -D ASSERT=all.  tally_read
      tests are kept to sequences shorter than KMERSIZE so the inner
      kmer loop never executes, avoiding the bug.
"""

from blazeseq import FastqRecord
from blazeqc.stats.kmer_content import KmerContent
from testing import assert_equal, assert_true, TestSuite


# ----- 1. Initialisation -----
# KmerContent[K] allocates 4^K kmer slots, each with one zero element.


def test_kmer_content_init_slot_count_k1():
    # KMERSIZE=1 → 4^1 = 4 slots
    var kc = KmerContent[1]()
    assert_equal(len(kc.kmers), 4)


def test_kmer_content_init_slot_count_k2():
    # KMERSIZE=2 → 4^2 = 16 slots
    var kc = KmerContent[2]()
    assert_equal(len(kc.kmers), 16)


def test_kmer_content_init_slot_count_k3():
    # KMERSIZE=3 → 4^3 = 64 slots
    var kc = KmerContent[3]()
    assert_equal(len(kc.kmers), 64)


def test_kmer_content_init_max_length_zero():
    var kc = KmerContent[2]()
    assert_equal(kc.max_length, 0)


def test_kmer_content_init_each_slot_has_one_element():
    var kc = KmerContent[2]()
    for i in range(len(kc.kmers)):
        assert_equal(len(kc.kmers[i]), 1)


def test_kmer_content_init_slot_values_zero():
    var kc = KmerContent[2]()
    for i in range(len(kc.kmers)):
        assert_equal(kc.kmers[i][0], 0)


# ----- 2. tally_read — sampling (1-in-50) -----
# Only reads where read_num % 50 == 0 are processed.


def test_kmer_tally_skips_non_multiples_of_50():
    # read_num=1 → skipped → max_length unchanged
    var kc = KmerContent[2]()
    var rec = FastqRecord("r1", "ACGTACGT", "IIIIIIII")
    kc.tally_read(rec, 1)
    assert_equal(kc.max_length, 0)


def test_kmer_tally_skips_read_num_49():
    var kc = KmerContent[2]()
    var rec = FastqRecord("r1", "ACGTACGT", "IIIIIIII")
    kc.tally_read(rec, 49)
    assert_equal(kc.max_length, 0)


def test_kmer_tally_processes_read_num_0():
    # read_num=0 → processed; sequence shorter than KMERSIZE → no kmer loop
    var kc = KmerContent[2]()
    # length 1 < KMERSIZE 2 → inner loop range(-1) → empty → kmer_to_index not called
    var rec = FastqRecord("r1", "A", "I")
    kc.tally_read(rec, 0)
    assert_equal(kc.max_length, 1)


def test_kmer_tally_processes_multiples_of_50():
    # read_num=50 → 50 % 50 == 0 → processed
    var kc = KmerContent[2]()
    var rec = FastqRecord("r1", "A", "I")
    kc.tally_read(rec, 50)
    assert_equal(kc.max_length, 1)


def test_kmer_tally_processes_read_num_100():
    var kc = KmerContent[2]()
    var rec = FastqRecord("r1", "AC", "II")  # length 2 = KMERSIZE → inner range(0) = empty
    kc.tally_read(rec, 100)
    assert_equal(kc.max_length, 2)


def test_kmer_tally_max_length_updates_on_longer_read():
    var kc = KmerContent[2]()
    var short_rec = FastqRecord("r1", "A", "I")
    var long_rec = FastqRecord("r2", "AC", "II")
    kc.tally_read(short_rec, 0)
    assert_equal(kc.max_length, 1)
    kc.tally_read(long_rec, 50)
    assert_equal(kc.max_length, 2)


def main():
    TestSuite.discover_tests[__functions_in_module()]().run()
