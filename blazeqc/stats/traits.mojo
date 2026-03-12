"""Stats traits: collection and summarization (core). FastQC output traits live in reporting_traits.mojo."""

from python import PythonObject
from blazeseq import FastqRecord, RefRecord
from blazeqc.stats.summary_utils import SummaryContext, GradeEntry


trait Collector(Copyable):
    fn tally_read(mut self, record: FastqRecord):
        ...

    fn tally_read(mut self, record: RefRecord):
        ...


trait Summarizer(Copyable, PlotOutput):
    """Main interface for summarization: prepare derived data and expose grades (pass/warn/fail)."""
    fn summerize(mut self, ctx: SummaryContext) raises:
        ...

    fn grade(self) raises -> GradeEntry:
        ...

    fn data_block_body(self) raises -> String:
        """Raw data lines (header + rows) for the data file block."""
        ...

    fn module_legend(self) -> String:
        """Display name for the module (e.g. for HTML report)."""
        ...

    fn panel_id(self) -> String:
        """HTML panel id (e.g. for result_panel)."""
        ...


trait PlotOutput(Copyable):
    """Output: plot results."""
    fn plot_result(self) raises -> PythonObject:
        ...


