"""Central pass/warn/fail limits for BlazeQC modules (FastQC-aligned)."""

# Per-tile quality: max absolute deviation (Phred) vs tile average
comptime TILE_WARN = 5.0
comptime TILE_ERROR = 10.0

# Adapter: % of reads with adapter k-mer at any position
comptime ADAPTER_WARN = 5.0
comptime ADAPTER_ERROR = 10.0

# Duplication: % remaining after deduplication (below = warn/fail)
comptime DUPLICATION_WARN = 70.0
comptime DUPLICATION_ERROR = 50.0

# Per-base quality: lower quartile and median (below = warn/fail)
comptime QUALITY_BASE_LOWER_WARN = 10.0
comptime QUALITY_BASE_LOWER_ERROR = 5.0
comptime QUALITY_BASE_MEDIAN_WARN = 25.0
comptime QUALITY_BASE_MEDIAN_ERROR = 20.0

# Per-sequence quality: most frequent mean Phred (below = warn/fail)
comptime QUALITY_SEQUENCE_WARN = 27.0
comptime QUALITY_SEQUENCE_ERROR = 20.0

# Per-base sequence content: max deviation A/T or C/G (%)
comptime SEQUENCE_WARN = 10.0
comptime SEQUENCE_ERROR = 20.0

# Per-sequence GC: max deviation from theoretical (%)
comptime GC_SEQUENCE_WARN = 15.0
comptime GC_SEQUENCE_ERROR = 30.0

# N content: % N at any position
comptime N_CONTENT_WARN = 5.0
comptime N_CONTENT_ERROR = 20.0
