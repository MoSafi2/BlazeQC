"""Summarization utilities: shared structs used by stats modules."""

from python import PythonObject
from blazeqc.html_maker import result_panel
from blazeqc.helpers import encode_img_b64


struct SummaryContext(Copyable):
    """Context passed into prepare() and output methods (e.g. num_reads for modules that need it)."""
    var num_reads: Int64
    var total_bases: Int64
    var file_name: String

    fn __init__(out self, num_reads: Int64, total_bases: Int64, file_name: String):
        self.num_reads = num_reads
        self.total_bases = total_bases
        self.file_name = file_name


struct GradeEntry(Copyable):
    """One (panel_legend, grade) entry from a module; multi-panel modules return multiple."""
    var panel_legend: String
    var grade: String

    fn __init__(out self, panel_legend: String, grade: String):
        self.panel_legend = panel_legend
        self.grade = grade


struct DefaultOutputter(Copyable):
    """Stateless helper: wrap (module_name, grade, body) into data-block text and (panel_id, grade, legend, figure) into result_panel."""

    fn __init__(out self):
        pass

    fn wrap_data_block(self, module_name: String, grade: String, body: String) -> String:
        """Returns >>module_name\tgrade\nbody>>END_MODULE\n (body should include its header line and data lines)."""
        return ">>{}\t{}\n{}>>END_MODULE\n".format(module_name, grade, body)

    fn make_panel(
        self, panel_id: String, grade: String, legend: String, figure: PythonObject
    ) raises -> result_panel:
        """Encodes figure to base64, returns result_panel for HTML."""
        var encoded = encode_img_b64(figure)
        return result_panel(panel_id, grade, legend, encoded, panel_type="image")

