from pathlib import Path
from blazeseq import FastqParser, FileReader, GZFile
from blazeqc.stats import FullStats
from time import perf_counter_ns
from sys import argv
from os import abort
from utils import Variant


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

    var fastq_file = String(args[1])
    if fastq_file.endswith(".gz"):
        var parser = FastqParser(GZFile(String(fastq_file), "r"), "generic")
        for record in parser.ref_records():
            n += 1
            stats.tally(record)
        stats.make_html(String(args[1]))

    elif fastq_file.endswith(".fastq"):
        var parser = FastqParser(FileReader(Path(fastq_file)), "generic")
        for record in parser.ref_records():
            n += 1
            stats.tally(record)
        stats.make_html(String(args[1]))
    else:
        print("Unsupported file type")
        abort()

    var t1 = perf_counter_ns()
    print("Time taken: ", t1 - t0, " ns")
    print("Number of reads: ", n)
    print("Total bases: ", stats.total_bases)
    print("Number of reads: ", n)
