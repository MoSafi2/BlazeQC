"""Stats traits: summarization and output interfaces (separate from collection implementations)."""

from python import PythonObject
from blazeqc.html_maker import result_panel
from blazeseq import FastqRecord, RefRecord
from blazeqc.stats.summary_utils import SummaryContext, GradeEntry, DefaultOutputter


trait Collector(Copyable):
    fn tally_read(mut self, record: FastqRecord):
        ...

    fn tally_read(mut self, record: RefRecord):
        ...


trait Summarizer(Copyable, PlotOutput):
    """Main interface for summarization: prepare derived data and expose grades (pass/warn/fail).
    Also provide different types of Output interfaces.
    """
    fn summerize(mut self, ctx: SummaryContext) raises:
        ...

    fn grade(self) raises -> GradeEntry:
        ...


trait PlotOutput(Copyable):
    """Output: plot results."""
    fn plot_result(self) raises -> PythonObject:
        ...


trait FastqcHtmlOutput(Copyable):
    """Output: result_panel for HTML report."""
    fn to_html(self) raises -> result_panel:
        ...



trait FastqcDataOutput(Copyable):
    """Output: FastQC-style data block text."""
    fn to_data_text(self, ctx: SummaryContext) raises -> String:
        ...


