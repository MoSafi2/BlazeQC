"""Shared analyser trait and constants (split from stats_.mojo)."""

from blazeseq import FastqRecord

# TODO: Make this dynamic
comptime py_lib: String = ".pixi/envs/default/lib/python3.12/site-packages/"


trait Analyser(Copyable):
    fn tally_read(mut self, record: FastqRecord):
        ...
