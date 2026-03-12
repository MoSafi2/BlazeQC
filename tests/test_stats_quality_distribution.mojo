"""Unit tests for blazeqc.stats.quality_distribution (pure Mojo)."""

from blazeseq import FastqRecord
from blazeqc.stats.quality_distribution import QualityCollector, QualityModule
from blazeqc.stats.summary_utils import SummaryContext
from testing import assert_equal, assert_true, TestSuite


# ----- 1. Initialisation (QualityCollector) -----


def test_quality_dist_init_max_length_zero():
    var qc = QualityCollector()
    assert_equal(qc.max_length, 0)


def test_quality_dist_init_min_qu():
    # min_qu starts at 128 (sentinel high value)
    var qc = QualityCollector()
    assert_equal(Int(qc.min_qu), 128)


def test_quality_dist_init_max_qu():
    # max_qu starts at 0
    var qc = QualityCollector()
    assert_equal(Int(qc.max_qu), 0)


def test_quality_dist_init_seq_length():
    var qc = QualityCollector()
    assert_equal(len(qc.qu_dist_seq), 128)


def test_quality_dist_init_seq_zeros():
    var qc = QualityCollector()
    for i in range(len(qc.qu_dist_seq)):
        assert_equal(qc.qu_dist_seq[i], 0)


# ----- 2. tally_read — state updates (QualityCollector) -----
# Quality string "IIII" → ASCII 73 per byte
# Quality "!" → ASCII 33, quality "A" → ASCII 65, quality "B" → ASCII 66


def test_quality_dist_tally_updates_max_length():
    var qc = QualityCollector()
    var rec = FastqRecord("r1", "ACGT", "IIII")
    qc.tally_read(rec)
    assert_equal(qc.max_length, 4)


def test_quality_dist_tally_max_length_grows():
    var qc = QualityCollector()
    var rec2 = FastqRecord("r1", "AC", "II")
    var rec4 = FastqRecord("r2", "ACGT", "IIII")
    qc.tally_read(rec2)
    assert_equal(qc.max_length, 2)
    qc.tally_read(rec4)
    assert_equal(qc.max_length, 4)


def test_quality_dist_tally_updates_max_qu():
    # 'I' = 73
    var qc = QualityCollector()
    var rec = FastqRecord("r1", "A", "I")
    qc.tally_read(rec)
    assert_equal(Int(qc.max_qu), 73)


def test_quality_dist_tally_updates_min_qu():
    # '!' = 33 — should become min_qu
    var qc = QualityCollector()
    var rec = FastqRecord("r1", "A", "!")
    qc.tally_read(rec)
    assert_equal(Int(qc.min_qu), 33)


def test_quality_dist_tally_min_max_mixed():
    # Two quality values: '!' (33) and 'I' (73)
    var qc = QualityCollector()
    var rec = FastqRecord("r1", "AC", "!I")
    qc.tally_read(rec)
    assert_equal(Int(qc.min_qu), 33)
    assert_equal(Int(qc.max_qu), 73)


def test_quality_dist_tally_matrix_increments():
    # After one read with quality 'I' (73), qu_dist.get(0, 73) == 1
    var qc = QualityCollector()
    var rec = FastqRecord("r1", "ACGT", "IIII")
    qc.tally_read(rec)
    assert_equal(qc.qu_dist.get(0, 73), 1)
    assert_equal(qc.qu_dist.get(1, 73), 1)
    assert_equal(qc.qu_dist.get(2, 73), 1)
    assert_equal(qc.qu_dist.get(3, 73), 1)


def test_quality_dist_tally_matrix_accumulates():
    # Two identical reads → each cell incremented twice
    var qc = QualityCollector()
    var rec = FastqRecord("r1", "AC", "II")
    qc.tally_read(rec)
    qc.tally_read(rec)
    assert_equal(qc.qu_dist.get(0, 73), 2)
    assert_equal(qc.qu_dist.get(1, 73), 2)


def test_quality_dist_tally_seq_average():
    # Uniform quality 'I' (73) across 4 bases → average = 73 → qu_dist_seq[73] += 1
    var qc = QualityCollector()
    var rec = FastqRecord("r1", "ACGT", "IIII")
    qc.tally_read(rec)
    assert_equal(qc.qu_dist_seq[73], 1)


def test_quality_dist_tally_seq_average_two_reads():
    var qc = QualityCollector()
    var rec = FastqRecord("r1", "ACGT", "IIII")
    qc.tally_read(rec)
    qc.tally_read(rec)
    assert_equal(qc.qu_dist_seq[73], 2)


# ----- 3. _guess_schema (QualityCollector) -----
# Discriminated by min_qu:
#   min_qu < 64       → Illumina v1.8  (LOWER=33, OFFSET=33)
#   min_qu == 65      → Illumina v1.3  (LOWER=64, OFFSET=64)
#   65 < min_qu ≤ 126 → Illumina v1.5  (LOWER=66, OFFSET=64)


def test_guess_schema_illumina_18():
    var qc = QualityCollector()
    # '!' = 33 → min_qu = 33 < 64 → Illumina v1.8
    var rec = FastqRecord("r1", "A", "!")
    qc.tally_read(rec)
    var schema = qc._guess_schema()
    assert_equal(Int(schema.OFFSET), 33)
    assert_equal(Int(schema.LOWER), 33)


def test_guess_schema_illumina_13():
    var qc = QualityCollector()
    # 'A' = 65 → min_qu = 65 == 64+1 → Illumina v1.3
    var rec = FastqRecord("r1", "A", "A")
    qc.tally_read(rec)
    var schema = qc._guess_schema()
    assert_equal(Int(schema.OFFSET), 64)
    assert_equal(Int(schema.LOWER), 64)


def test_guess_schema_illumina_15():
    var qc = QualityCollector()
    # 'B' = 66 → min_qu = 66; 66 != 65 and 66 <= 126 → Illumina v1.5
    var rec = FastqRecord("r1", "A", "B")
    qc.tally_read(rec)
    var schema = qc._guess_schema()
    assert_equal(Int(schema.LOWER), 66)
    assert_equal(Int(schema.OFFSET), 64)


# ----- Status (QualityModule: feed + summerize then grade) -----


def test_quality_dist_status_per_base_pass():
    var module = QualityModule()
    # Use high Phred quality so quartile/median pass (e.g. '?' = 63 in ASCII -> Phred 30 in Sanger)
    for _ in range(100):
        var rec = FastqRecord("r", "ACGTACGTACGT", "????????????")
        module.collector.tally_read(rec)
    var ctx = SummaryContext(100, 1200, "test")
    module.summarizer_base.feed(module.collector)
    module.summarizer_seq.feed(module.collector)
    module.summarizer_base.summerize(ctx)
    module.summarizer_seq.summerize(ctx)
    assert_equal(module.summarizer_base.grade().grade, "pass")


def test_quality_dist_status_per_sequence_pass():
    var module = QualityModule()
    for _ in range(100):
        var rec = FastqRecord("r", "ACGTACGTACGT", "????????????")
        module.collector.tally_read(rec)
    var ctx = SummaryContext(100, 1200, "test")
    module.summarizer_base.feed(module.collector)
    module.summarizer_seq.feed(module.collector)
    module.summarizer_base.summerize(ctx)
    module.summarizer_seq.summerize(ctx)
    assert_equal(module.summarizer_seq.grade().grade, "pass")


def main():
    TestSuite.discover_tests[__functions_in_module()]().run()
