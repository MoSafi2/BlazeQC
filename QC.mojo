from pathlib import Path
from blazeseq import FastqParser, FileReader
from blazeqc.stats import FullStats
from time import perf_counter_ns
from sys import argv
from os import abort


fn main() raises:
    args = argv()

    var parser = FastqParser(FileReader(Path(args[1])), "generic")
    var stats = FullStats()

    var n = 0

    t0 = perf_counter_ns()
    for record in parser.ref_records():
        n += 1
        stats.tally(record)
    print(n)
    stats.make_html(String(args[1]))
