from blazeseq.quality_schema import (
    QualitySchema,
    generic_schema,
    sanger_schema,
    solexa_schema,
    illumina_1_3_schema,
    illumina_1_5_schema,
    illumina_1_8_schema,
)
from sys.info import simd_width_of

comptime KB = 1024
comptime MB = 1024 * KB
comptime GB = 1024 * MB


comptime USE_SIMD = True
comptime read_header = 64
comptime quality_header = 43
comptime new_line = 10
comptime carriage_return = 13

comptime simd_width: Int = simd_width_of[UInt8]()

comptime DEFAULT_CAPACITY = 64 * KB
comptime MAX_SHIFT = 30
comptime MAX_CAPACITY = 2**MAX_SHIFT
