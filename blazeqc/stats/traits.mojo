"""Stats traits: summarization and output interfaces (separate from collection implementations)."""

from collections.dict import Dict
from collections.list import List
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


trait FastqcHtmlOutput(Copyable):
    """Output: result_panel for HTML report. Caller calls plot_result() and passes figures to to_html_panels."""
    fn to_html(self) raises -> result_panel:
        ...

    fn to_html_panels(self, figures: List[PythonObject]) raises -> Dict[String, result_panel]:
        """Build panels from pre-computed figures; returns dict keyed by module_legend."""
        ...


trait FastqcDataOutput(Copyable):
    """Output: FastQC-style data block text."""
    fn to_data_text(self, ctx: SummaryContext) raises -> String:
        ...


trait ModuleReport(FastqcDataOutput, FastqcHtmlOutput, Copyable):
    """Container module: to_data_text, to_html_panels, and plot_result (delegates to summarizer)."""
    fn plot_result(self) raises -> PythonObject:
        ...


