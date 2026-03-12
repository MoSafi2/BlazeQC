"""Summarization utilities: shared structs used by stats modules."""

from collections.dict import Dict
from collections.list import List
from python import PythonObject
from blazeqc.html_maker import result_panel
from blazeqc.helpers import encode_img_b64


struct DataEntry(Copyable, ImplicitlyCopyable):
    """One module's data block: (module_legend, grade, body) for FastQC data file."""
    var module_legend: String
    var grade: String
    var body: String

    fn __init__(out self, module_legend: String, grade: String, body: String):
        self.module_legend = module_legend
        self.grade = grade
        self.body = body


struct PanelEntry(Copyable, ImplicitlyCopyable):
    """One panel for HTML report: image (figure) or table (html_content)."""
    var panel_id: String
    var grade: String
    var legend: String
    var panel_type: String  # "image" or "table"
    var figure: PythonObject  # for image; use Python.none() for table
    var html_content: String  # for table; empty for image

    fn __init__(
        out self,
        panel_id: String,
        grade: String,
        legend: String,
        panel_type: String,
        figure: PythonObject,
        html_content: String,
    ):
        self.panel_id = panel_id
        self.grade = grade
        self.legend = legend
        self.panel_type = panel_type
        self.figure = figure
        self.html_content = html_content


fn wrap_data_blocks(entries: List[DataEntry]) -> String:
    """Build FastQC-style data block text from a list of data entries."""
    var out = DefaultOutputter()
    var result = String("")
    for i in range(len(entries)):
        var e = entries[i].copy()
        result += out.wrap_data_block(e.module_legend, e.grade, e.body)
    return result


fn make_panels(entries: List[PanelEntry]) raises -> Dict[String, result_panel]:
    """Build dict of result_panel keyed by legend from panel entries (image or table)."""
    var out = DefaultOutputter()
    var d = Dict[String, result_panel]()
    for i in range(len(entries)):
        var e = entries[i].copy()
        var p: result_panel
        if e.panel_type == "image":
            p = out.make_panel(e.panel_id, e.grade, e.legend, e.figure)
        else:
            p = result_panel(e.panel_id, e.grade, e.legend, e.html_content, panel_type="table")
        d[p.legand] = p^
    return d^


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

