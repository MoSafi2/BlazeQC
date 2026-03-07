"""Unit tests for blazeqc.stats.adapter_content (pure Mojo)."""

from blazeseq import FastqRecord
from blazeqc.stats.adapter_content import AdapterContent
from blazeqc.config import hash_list as get_hash_list
from blazeqc.helpers import Matrix2D
from testing import assert_equal, assert_true, TestSuite


# ----- 1. Initialisation -----


def test_adapter_content_init_max_length_zero():
    var hashes = List[UInt64]()
    hashes.append(12345)
    var ac = AdapterContent(hashes^, 12)
    assert_equal(ac.max_length, 0)


def test_adapter_content_init_hash_counts_shape():
    # hash_counts should be Matrix2D(num_hashes, 1)
    var hashes = List[UInt64]()
    hashes.append(1)
    hashes.append(2)
    var ac = AdapterContent(hashes^, 5)
    var sh = ac.hash_counts.shape()
    assert_equal(sh[0], 2)   # one row per hash
    assert_equal(sh[1], 1)   # single initial column


def test_adapter_content_init_hash_list_stored():
    var hashes = List[UInt64]()
    hashes.append(0xDEADBEEF)
    var ac = AdapterContent(hashes^, 5)
    assert_equal(len(ac.hash_list), 1)
    assert_equal(ac.hash_list[0], 0xDEADBEEF)


def test_adapter_content_init_kmer_len_stored():
    # bits=3 (default), so max kmer_len = 64 // 3 = 21
    var hashes = List[UInt64]()
    hashes.append(1)
    var ac = AdapterContent(hashes^, 12)
    assert_equal(ac.kmer_len, 12)


def test_adapter_content_init_kmer_len_capped():
    # kmer_len larger than 64 // bits (21) is capped
    var hashes = List[UInt64]()
    hashes.append(1)
    var ac = AdapterContent(hashes^, 100)
    # kmer_len = min(100, 64 // 3) = min(100, 21) = 21
    assert_equal(ac.kmer_len, 21)


def test_adapter_content_init_kmer_len_zero():
    var hashes = List[UInt64]()
    hashes.append(1)
    var ac = AdapterContent(hashes^)  # default kmer_len=0
    assert_equal(ac.kmer_len, 0)


# ----- 2. tally_read — basic state update -----


def test_adapter_content_tally_updates_max_length():
    var hashes = List[UInt64]()
    hashes.append(9999)
    var ac = AdapterContent(hashes^, 12)
    var rec = FastqRecord("r1", "ACGTACGT", "IIIIIIII")
    ac.tally_read(rec)
    assert_equal(ac.max_length, 8)


def test_adapter_content_tally_max_length_grows():
    var hashes = List[UInt64]()
    hashes.append(9999)
    var ac = AdapterContent(hashes^, 12)
    var rec4 = FastqRecord("r1", "ACGT", "IIII")
    var rec8 = FastqRecord("r2", "ACGTACGT", "IIIIIIII")
    ac.tally_read(rec4)
    assert_equal(ac.max_length, 4)
    ac.tally_read(rec8)
    assert_equal(ac.max_length, 8)


def test_adapter_content_tally_no_match_counts_zero():
    # A hash that will never match (0xDEADBEEF) → counts stay zero
    var hashes = List[UInt64]()
    hashes.append(0xDEADBEEF)
    var ac = AdapterContent(hashes^, 4)
    var rec = FastqRecord("r1", "ACGTACGT", "IIIIIIII")
    ac.tally_read(rec)
    # Verify the counts for hash 0 across all positions are 0
    var total: Int64 = 0
    for col in range(ac.hash_counts.shape()[1]):
        total += ac.hash_counts.get(0, col)
    assert_equal(total, 0)


def test_adapter_content_empty_hash_list_no_crash():
    # With no hashes, tally_read should not crash (the _check_hashes loop is skipped)
    var hashes = List[UInt64]()
    var ac = AdapterContent(hashes^, 4)
    var rec = FastqRecord("r1", "ACGT", "IIII")
    ac.tally_read(rec)
    assert_equal(ac.max_length, 4)


# ----- 3. tally_read — canonical adapter match -----
# Uses hash_list() from config (same hashes used in the real QC pipeline).
# The Illumina Universal Adapter is "AGATCGGAAGAG" (12 chars).
# Feeding a record that starts with this sequence should increment row 0
# of hash_counts at position 12 (after the kmer window fills).


def test_adapter_content_illumina_universal_match():
    var hashes = get_hash_list()
    # kmer_len=12 matches the 12-char adapter sequences in config
    var ac = AdapterContent(hashes^, 12)
    # Record whose first 12 bases are the Illumina Universal Adapter
    var rec = FastqRecord(
        "r1", "AGATCGGAAGAGACGT", "IIIIIIIIIIIIIIII"
    )
    ac.tally_read(rec)
    # At least one count in row 0 (Illumina Universal Adapter) must be > 0
    var hit: Int64 = 0
    for col in range(ac.hash_counts.shape()[1]):
        hit += ac.hash_counts.get(0, col)
    assert_true(hit > 0)


def main():
    TestSuite.discover_tests[__functions_in_module()]().run()
