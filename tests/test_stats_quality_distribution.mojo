"""Unit tests for blazeqc.stats.quality_distribution (pure Mojo)."""

from blazeseq import FastqRecord
from blazeqc.stats.quality_distribution import QualityDistribution
from testing import assert_equal, assert_true, TestSuite


# ----- 1. Initialisation -----


def test_quality_dist_init_max_length_zero():
    var qd = QualityDistribution()
    assert_equal(qd.max_length, 0)


def test_quality_dist_init_min_qu():
    # min_qu starts at 128 (sentinel high value)
    var qd = QualityDistribution()
    assert_equal(Int(qd.min_qu), 128)


def test_quality_dist_init_max_qu():
    # max_qu starts at 0
    var qd = QualityDistribution()
    assert_equal(Int(qd.max_qu), 0)


def test_quality_dist_init_seq_length():
    var qd = QualityDistribution()
    assert_equal(len(qd.qu_dist_seq), 128)


def test_quality_dist_init_seq_zeros():
    var qd = QualityDistribution()
    for i in range(len(qd.qu_dist_seq)):
        assert_equal(qd.qu_dist_seq[i], 0)


# ----- 2. tally_read — state updates -----
# Quality string "IIII" → ASCII 73 per byte
# Quality "!" → ASCII 33, quality "A" → ASCII 65, quality "B" → ASCII 66


def test_quality_dist_tally_updates_max_length():
    var qd = QualityDistribution()
    var rec = FastqRecord("r1", "ACGT", "IIII")
    qd.tally_read(rec)
    assert_equal(qd.max_length, 4)


def test_quality_dist_tally_max_length_grows():
    var qd = QualityDistribution()
    var rec2 = FastqRecord("r1", "AC", "II")
    var rec4 = FastqRecord("r2", "ACGT", "IIII")
    qd.tally_read(rec2)
    assert_equal(qd.max_length, 2)
    qd.tally_read(rec4)
    assert_equal(qd.max_length, 4)


def test_quality_dist_tally_updates_max_qu():
    # 'I' = 73
    var qd = QualityDistribution()
    var rec = FastqRecord("r1", "A", "I")
    qd.tally_read(rec)
    assert_equal(Int(qd.max_qu), 73)


def test_quality_dist_tally_updates_min_qu():
    # '!' = 33 — should become min_qu
    var qd = QualityDistribution()
    var rec = FastqRecord("r1", "A", "!")
    qd.tally_read(rec)
    assert_equal(Int(qd.min_qu), 33)


def test_quality_dist_tally_min_max_mixed():
    # Two quality values: '!' (33) and 'I' (73)
    var qd = QualityDistribution()
    var rec = FastqRecord("r1", "AC", "!I")
    qd.tally_read(rec)
    assert_equal(Int(qd.min_qu), 33)
    assert_equal(Int(qd.max_qu), 73)


def test_quality_dist_tally_matrix_increments():
    # After one read with quality 'I' (73), qu_dist.get(0, 73) == 1
    var qd = QualityDistribution()
    var rec = FastqRecord("r1", "ACGT", "IIII")
    qd.tally_read(rec)
    assert_equal(qd.qu_dist.get(0, 73), 1)
    assert_equal(qd.qu_dist.get(1, 73), 1)
    assert_equal(qd.qu_dist.get(2, 73), 1)
    assert_equal(qd.qu_dist.get(3, 73), 1)


def test_quality_dist_tally_matrix_accumulates():
    # Two identical reads → each cell incremented twice
    var qd = QualityDistribution()
    var rec = FastqRecord("r1", "AC", "II")
    qd.tally_read(rec)
    qd.tally_read(rec)
    assert_equal(qd.qu_dist.get(0, 73), 2)
    assert_equal(qd.qu_dist.get(1, 73), 2)


def test_quality_dist_tally_seq_average():
    # Uniform quality 'I' (73) across 4 bases → average = 73 → qu_dist_seq[73] += 1
    var qd = QualityDistribution()
    var rec = FastqRecord("r1", "ACGT", "IIII")
    qd.tally_read(rec)
    assert_equal(qd.qu_dist_seq[73], 1)


def test_quality_dist_tally_seq_average_two_reads():
    var qd = QualityDistribution()
    var rec = FastqRecord("r1", "ACGT", "IIII")
    qd.tally_read(rec)
    qd.tally_read(rec)
    assert_equal(qd.qu_dist_seq[73], 2)


# ----- 3. _guess_schema -----
# Discriminated by min_qu:
#   min_qu < 64       → Illumina v1.8  (LOWER=33, OFFSET=33)
#   min_qu == 65      → Illumina v1.3  (LOWER=64, OFFSET=64)
#   65 < min_qu ≤ 126 → Illumina v1.5  (LOWER=66, OFFSET=64)


def test_guess_schema_illumina_18():
    var qd = QualityDistribution()
    # '!' = 33 → min_qu = 33 < 64 → Illumina v1.8
    var rec = FastqRecord("r1", "A", "!")
    qd.tally_read(rec)
    var schema = qd._guess_schema()
    assert_equal(Int(schema.OFFSET), 33)
    assert_equal(Int(schema.LOWER), 33)


def test_guess_schema_illumina_13():
    var qd = QualityDistribution()
    # 'A' = 65 → min_qu = 65 == 64+1 → Illumina v1.3
    var rec = FastqRecord("r1", "A", "A")
    qd.tally_read(rec)
    var schema = qd._guess_schema()
    assert_equal(Int(schema.OFFSET), 64)
    assert_equal(Int(schema.LOWER), 64)


def test_guess_schema_illumina_15():
    var qd = QualityDistribution()
    # 'B' = 66 → min_qu = 66; 66 != 65 and 66 <= 126 → Illumina v1.5
    var rec = FastqRecord("r1", "A", "B")
    qd.tally_read(rec)
    var schema = qd._guess_schema()
    assert_equal(Int(schema.LOWER), 66)
    assert_equal(Int(schema.OFFSET), 64)


# ----- _get_status_per_base, _get_status_per_sequence -----


def test_quality_dist_status_per_base_pass():
    var qd = QualityDistribution()
    for _ in range(100):
        var rec = FastqRecord("r", "ACGTACGTACGT", "IIIIIIIIIIII")
        qd.tally_read(rec)
    assert_equal(qd._get_status_per_base(), "pass")


def test_quality_dist_status_per_sequence_pass():
    var qd = QualityDistribution()
    for _ in range(100):
        var rec = FastqRecord("r", "ACGTACGTACGT", "IIIIIIIIIIII")
        qd.tally_read(rec)
    assert_equal(qd._get_status_per_sequence(), "pass")


def main():
    TestSuite.discover_tests[__functions_in_module()]().run()
