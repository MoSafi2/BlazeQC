from pathlib import Path
from blazeseq import FastqParser, FileReader, RapidgzipReader, ParserConfig
from blazeqc.stats import FullStats
from time import perf_counter_ns
from sys import argv
from os import abort, stat
from utils import Variant

# Report progress every this many bytes of stream position
comptime PROGRESS_BYTE_INTERVAL: Int = 1_000_000  # 1 MB

comptime config = ParserConfig(check_ascii=False, check_quality=False, buffer_capacity=64*1024)

fn main() raises:
    args = argv()

    if len(args) != 2:
        print("Usage: QC.mojo <fastq_file>")
        abort()
    if not Path(args[1]).exists():
        print("File does not exist")
        abort()

    var stats = FullStats()
    var n = 0
    t0 = perf_counter_ns()

    var input_path = Path(args[1])
    var total_size = stat(input_path).st_size

    var fastq_file = String(args[1])
    var next_threshold: Int64
    if fastq_file.endswith(".gz"):
        next_threshold = PROGRESS_BYTE_INTERVAL
        var parser = FastqParser[config=config](RapidgzipReader(fastq_file, parallelism=0), "generic")
        for record in parser.ref_records():
            n += 1
            stats.tally(record)
            var pos = parser.buffer.source._file.tell()
            if pos >= next_threshold:
                if total_size > 0:
                    print("\rParsing & tallying: ", Int(100 * pos // total_size), "% (", n, " reads)", end="", flush=True)
                else:
                    print("\rParsing & tallying: ", n, " reads", end="", flush=True)
                next_threshold = pos + PROGRESS_BYTE_INTERVAL
        if total_size > 0:
            print("\rParsing & tallying: 100% (", n, " reads)", end="", flush=True)
        else:
            print("\rParsing & tallying: ", n, " reads", end="", flush=True)
        print()
        stats.make_html(String(args[1]))

    elif fastq_file.endswith(".fastq") or fastq_file.endswith(".fq"):
        next_threshold = PROGRESS_BYTE_INTERVAL
        var parser = FastqParser[config=config](FileReader(Path(fastq_file)), "generic")
        for record in parser.ref_records():
            n += 1
            stats.tally(record)
            var pos = parser.get_file_position()
            if pos >= next_threshold:
                if total_size > 0:
                    print("\rParsing & tallying: ", Int(100 * pos // total_size), "% (", n, " reads)", end="", flush=True)
                else:
                    print("\rParsing & tallying: ", n, " reads", end="", flush=True)
                next_threshold = pos + PROGRESS_BYTE_INTERVAL
        if total_size > 0:
            print("\rParsing & tallying: 100% (", n, " reads)", end="", flush=True)
        else:
            print("\rParsing & tallying: ", n, " reads", end="", flush=True)
        print()
        stats.make_html(String(args[1]))
    else:
        print("Unsupported file type")
        abort()

    var t1 = perf_counter_ns()
    print("Time taken: ", t1 - t0, " ns")
    print("Number of reads: ", n)
    print("Total bases: ", stats.total_bases)
    print("Number of reads: ", n)
