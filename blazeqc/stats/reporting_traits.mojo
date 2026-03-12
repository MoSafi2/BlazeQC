"""FastQC-specific output traits: data file and HTML report. Depends on summary_utils and html_maker."""

from collections.dict import Dict
from collections.list import List
from python import PythonObject
from blazeqc.html_maker import result_panel
from blazeqc.stats.summary_utils import (
    SummaryContext,
    DataEntry,
    PanelEntry,
    wrap_data_blocks,
    make_panels,
)


trait FastqcDataOutput(Copyable):
    """Output: FastQC-style data block text."""
    fn data_entries(self, ctx: SummaryContext) raises -> List[DataEntry]:
        """Return list of data entries for default to_data_text."""
        ...

    fn to_data_text(self, ctx: SummaryContext) raises -> String:
        """Default implementation builds data block from data_entries."""
        return wrap_data_blocks(self.data_entries(ctx))


trait FastqcHtmlOutput(Copyable):
    """Output: result_panel for HTML report. Caller calls plot_result() and passes figures to to_html_panels."""
    fn to_html(self) raises -> result_panel:
        """Default: not implemented; use to_html_panels(figures) for report modules."""
        raise Error("to_html not implemented; use to_html_panels")

    fn panel_entries(self, figures: List[PythonObject]) raises -> List[PanelEntry]:
        """Return list of panel entries (image or table) for default to_html_panels."""
        ...

    fn to_html_panels(self, figures: List[PythonObject]) raises -> Dict[String, result_panel]:
        """Build panels from pre-computed figures; returns dict keyed by module_legend. Default uses panel_entries + make_panels."""
        return make_panels(self.panel_entries(figures))


trait ModuleReport(FastqcDataOutput, FastqcHtmlOutput, Copyable):
    """Container module: to_data_text, to_html_panels, and plot_result (delegates to summarizer)."""
    fn plot_result(self) raises -> PythonObject:
        ...
