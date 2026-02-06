from blazeseq.parser import RecordParser
from blazeseq.iostream import FileReader
from blazeqc.stats import FullStats
from time import perf_counter_ns
from sys import argv
from os import abort


fn main() raises:
    args = argv()

    var parser = RecordParser[check_ascii=False, check_quality=False](
        FileReader(String(args[1]))
    )
    var stats = FullStats()

    var n = 0

    t0 = perf_counter_ns()
    while True:
        try:
            var record = parser.next()
            n += 1
            if record is not None:
                stats.tally(record.take())
            else:
                break
        except Error:
            print(Error)
            t1 = perf_counter_ns()
            stats.make_html(String(args[1]))
            t2 = perf_counter_ns()
            print(n)
            print("Total tally time: ", Float64(t1 - t0) / 1e9, "s")
            print("Total Plot time", Float64(t2 - t1) / 1e9, "s")
            break


    stats.make_html(String(args[1]))