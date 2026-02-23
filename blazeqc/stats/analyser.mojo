"""Shared analyser trait and constants (split from stats_.mojo)."""

from blazeseq import FastqRecord, RefRecord

trait Analyser(Copyable):
    fn tally_read(mut self, record: FastqRecord):
        ...

    fn tally_read(mut self, record: RefRecord):
        ...
