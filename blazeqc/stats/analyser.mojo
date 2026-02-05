"""Shared analyser trait and constants (split from stats_.mojo)."""

from blazeseq import FastqRecord

trait Analyser(Copyable):
    fn tally_read(mut self, record: FastqRecord):
        ...
