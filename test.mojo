import blazeseq
from blazeseq import RapidgzipReader

fn main():
    var reader = RapidgzipReader("test.fastq.gz")
    for record in reader:
        print(record)