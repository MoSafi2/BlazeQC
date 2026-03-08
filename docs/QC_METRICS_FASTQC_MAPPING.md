# Quality Control Metrics: BlazeQC ↔ FastQC

Mapping of BlazeSeq (BlazeQC) quality control metrics to [FastQC modules](https://github.com/s-andrews/FastQC/tree/master/uk/ac/babraham/FastQC/Modules).

## Table


| **BlazeQC (BlazeSeq)**                                                | **FastQC Module**                   | **Notes**                                                                                                         |
| --------------------------------------------------------------------- | ----------------------------------- | ----------------------------------------------------------------------------------------------------------------- |
| **Base statistics** (`make_base_stats()` in `full_stats.mojo`)        | **BasicStats.java**                 | Filename, encoding, total sequences, total bases, sequence length, %GC.                                           |
| **QualityDistribution** (`blazeqc/stats/quality_distribution.mojo`)   | **PerBaseQualityScores.java**       | Per-base quality (quality by position).                                                                           |
| **QualityDistribution** (same struct)                                 | **PerSequenceQualityScores.java**   | Per-sequence quality (quality distribution per read).                                                             |
| **PerTileQuality** (`blazeqc/stats/tile_quality.mojo`)                | **PerTileQualityScores.java**       | Quality by tile (flow cell position).                                                                             |
| **BasepairDistribution** (`blazeqc/stats/basepair_distribution.mojo`) | **PerBaseSequenceContent.java**     | Per-base A/T/G/C content (base pair distribution).                                                                |
| **BasepairDistribution** (same struct)                                | **NContent.java**                   | Per-base N content (reported as "Base pair N percentage" in BlazeQC).                                             |
| **CGContent** (`blazeqc/stats/cg_content.mojo`)                       | **PerSequenceGCContent.java**       | Per-sequence GC content distribution.                                                                             |
| **LengthDistribution** (`blazeqc/stats/length_distribution.mojo`)     | **SequenceLengthDistribution.java** | Distribution of read lengths.                                                                                     |
| **DupReads** (`blazeqc/stats/duplication.mojo`)                       | **DuplicationLevel.java**           | Duplication levels (sequence duplication levels).                                                                 |
| **DupReads** (same struct, via `over_represented`)                    | **OverRepresentedSeqs.java**        | Overrepresented sequences.                                                                                        |
| **AdapterContent** (`blazeqc/stats/adapter_content.mojo`)             | **AdapterContent.java**             | Adapter sequence content along reads.                                                                             |
| **KmerContent** (`blazeqc/stats/kmer_content.mojo`)                   | **KmerContent.java**                | K-mer overrepresentation; **implemented in BlazeQC but not in the default report** (`FullStats` does not use it). |


## Summary

- **In the default BlazeQC report:** Base statistics, per-base and per-sequence quality, per-tile quality, per-base sequence content, N content, per-sequence GC content, length distribution, duplication levels, overrepresented sequences, and adapter content.
- **FastQC modules with a direct BlazeQC counterpart:** All of the above; BlazeQC’s `BasepairDistribution` covers both PerBaseSequenceContent and NContent; `QualityDistribution` covers both per-base and per-sequence quality; `DupReads` covers both DuplicationLevel and OverRepresentedSeqs.
- **Present in BlazeQC but not in the default report:** `KmerContent` (has tests and implementation; not wired into `FullStats`).

