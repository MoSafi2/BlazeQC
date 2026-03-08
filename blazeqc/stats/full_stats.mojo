"""FullStats aggregate (split from stats_.mojo)."""

from collections.dict import Dict
from collections.list import List
from python import Python, PythonObject
from blazeseq import FastqRecord, RefRecord
from blazeqc.stats.analyser import Analyser
from blazeqc.stats.basepair_distribution import BasepairDistribution
from blazeqc.stats.cg_content import CGContent
from blazeqc.stats.duplication import DupReads
from blazeqc.stats.length_distribution import LengthDistribution
from blazeqc.stats.quality_distribution import QualityDistribution
from blazeqc.stats.tile_quality import PerTileQuality
from blazeqc.stats.adapter_content import AdapterContent
from blazeqc.config import hash_list
from blazeqc.html_maker import (
    result_panel,
    insert_result_panel,
    html_template,
    format_length,
)

# ----- Report data (FastQC-style data file) -----
# FullStats.write_data(file_name) writes ##BlazeQC, basic stats, then each module's write_module_data(f).

@fieldwise_init
struct FullStats(Copyable):
    var num_reads: Int64
    var total_bases: Int64
    var bp_dist: BasepairDistribution
    var len_dist: LengthDistribution
    var qu_dist: QualityDistribution
    var cg_content: CGContent
    var dup_reads: DupReads
    var tile_qual: PerTileQuality
    var adpt_cont: AdapterContent[3]

    fn __init__(out self) raises:
        self.num_reads = 0
        self.total_bases = 0
        self.len_dist = LengthDistribution()
        self.bp_dist = BasepairDistribution()
        self.cg_content = CGContent()
        self.qu_dist = QualityDistribution()
        self.dup_reads = DupReads()
        self.tile_qual = PerTileQuality()
        self.adpt_cont = AdapterContent[bits=3](hash_list(), 12)

    @always_inline
    fn tally(mut self, record: FastqRecord):
        self.num_reads += 1
        self.total_bases += len(record)
        self.bp_dist.tally_read(record)
        self.len_dist.tally_read(record)
        self.cg_content.tally_read(record)  # Almost Free
        self.dup_reads.tally_read(record)
        self.qu_dist.tally_read(record)
        self.adpt_cont.tally_read(record, self.num_reads)
        self.tile_qual.tally_read(record)

    @always_inline
    fn tally(mut self, record: RefRecord):
        self.num_reads += 1
        self.total_bases += len(record)
        self.bp_dist.tally_read(record)
        self.len_dist.tally_read(record)
        self.cg_content.tally_read(record)
        self.dup_reads.tally_read(record)
        self.qu_dist.tally_read(record)
        self.adpt_cont.tally_read(record, self.num_reads)
        self.tile_qual.tally_read(record)

    @always_inline
    fn make_base_stats(self) raises -> result_panel:
        var sum: Int64 = 0
        for i in range(len(self.cg_content.cg_content)):
            sum += self.cg_content.cg_content[i] * i
        var avg_cg = sum / self.num_reads
        var schema = self.qu_dist._guess_schema()

        var total_bases = format_length(Float64(self.total_bases))

        var lengths: String
        if self.bp_dist.max_length == self.bp_dist.min_length:
            lengths = String(self.bp_dist.max_length)
        else:
            lengths = "{}-{}".format(
                self.bp_dist.min_length, self.bp_dist.max_length
            )

        var table_template = """
                        <table>
                        <thead>
                            <tr>
                                <th>Measure</th>
                                <th>Value</th>
                            </tr>
                        </thead>
                        <tbody>
                        <tr>
                            <td>Filename</td>
                            <td><<filename>></td>
                        <tr>
                            <td>Encoding</td>
                            <td>{}</td>
                        </tr>
                        <tr>
                            <td>Total Sequences</td>
                            <td>{}</td>
                        </tr>
                        <tr>
                            <td>Total Bases	</td>
                            <td>{}</td>
                        </tr>
                        <tr>
                            <td>Sequence length</td>
                            <td>{}
                        </td>
                        <tr>
                            <td>%GC	</td>
                            <td>{}</td>
                        </tr>
                        </tbody>
                        </table>
                        """.format(
            schema.SCHEMA,
            self.num_reads,
            total_bases,
            lengths.split("/")[-1],
            avg_cg,
        )
        var res = result_panel(
            "base_stats",
            "pass",
            "Basic Statistics",
            table_template,
            panel_type="table",
        )
        return res^

    @always_inline
    fn plot(mut self) raises -> List[PythonObject]:
        var plots = List[PythonObject]()

        var bp_plots = self.bp_dist.plot(self.num_reads)
        plots.append(bp_plots[0])
        plots.append(bp_plots[1])
        plots.append(self.cg_content.plot())
        plots.append(self.len_dist.plot())
        var dup_plot_result = self.dup_reads.plot(Int(self.num_reads))
        plots.append(dup_plot_result[0])
        var qu_plots = self.qu_dist.plot()
        plots.append(qu_plots[0])
        plots.append(qu_plots[1])
        plots.append(self.tile_qual.plot())
        plots.append(self.adpt_cont.plot(self.num_reads))

        return plots^

    fn prepare_data(mut self, file_name: String) raises:
        """Fill all module caches. Call before write_data and build_panels."""
        self.qu_dist.prepare_data()
        self.tile_qual.prepare_data()
        self.bp_dist.prepare_data(self.num_reads)
        self.cg_content.prepare_data()
        self.len_dist.prepare_data()
        self.dup_reads.prepare_data(Int(self.num_reads))
        self.adpt_cont.prepare_data(self.num_reads)

    fn write_data(mut self, file_name: String) raises:
        """Write FastQC-style data file. Calls prepare_data then writes each module block."""
        self.prepare_data(file_name)

        with open("{}_blazeseq_data.txt".format(file_name), "w") as f:
            f.write("##BlazeQC\t0.12.1\n")

            # Basic Statistics
            var sum: Int64 = 0
            for i in range(len(self.cg_content.cg_content)):
                sum += self.cg_content.cg_content[i] * i
            var avg_cg = sum / self.num_reads
            var schema = self.qu_dist._guess_schema()
            var total_bases_str = format_length(Float64(self.total_bases))
            var lengths: String
            if self.bp_dist.max_length == self.bp_dist.min_length:
                lengths = String(self.bp_dist.max_length)
            else:
                lengths = "{}-{}".format(
                    self.bp_dist.min_length, self.bp_dist.max_length
                )
            var filename_display = file_name
            if file_name.find("/") >= 0:
                filename_display = String(file_name.split("/")[-1])
            f.write(">>Basic Statistics\tpass\n")
            f.write("#Measure\tValue\n")
            f.write("Filename\t{}\n".format(filename_display))
            f.write("File type\tConventional base calls\n")
            f.write("Encoding\t{}\n".format(schema.SCHEMA))
            f.write("Total Sequences\t{}\n".format(self.num_reads))
            f.write("Total Bases\t{}\n".format(total_bases_str))
            f.write("Sequences flagged as poor quality\t0\n")
            f.write("Sequence length\t{}\n".format(lengths))
            f.write("%GC\t{}\n".format(avg_cg))
            f.write(">>END_MODULE\n")

            # Module blocks in panel order (each returns its block text)
            f.write(self.qu_dist.get_module_data())
            f.write(self.tile_qual.get_module_data())
            f.write(self.bp_dist.get_module_data(self.num_reads))
            f.write(self.cg_content.get_module_data())
            f.write(self.len_dist.get_module_data())
            f.write(self.dup_reads.get_module_data(Int(self.num_reads)))
            f.write(self.adpt_cont.get_module_data(self.num_reads))

    fn build_panels(mut self) raises -> Dict[String, result_panel]:
        var panels = Dict[String, result_panel]()
        var base_stats = self.make_base_stats()
        panels[base_stats.legand] = base_stats^

        var qu_html = self.qu_dist.make_html()
        var per_base_quality_panel = qu_html[0].copy()
        var per_sequence_quality_panel = qu_html[1].copy()
        panels[per_sequence_quality_panel.legand] = per_sequence_quality_panel^
        panels[per_base_quality_panel.legand] = per_base_quality_panel^

        var tile_quality = self.tile_qual.make_html()
        panels[tile_quality.legand] = tile_quality^

        var bp_html = self.bp_dist.make_html(self.num_reads)
        var base_pair_N_percentage = bp_html[0].copy()
        var base_pair_distribution = bp_html[1].copy()
        panels[base_pair_distribution.legand] = base_pair_distribution^
        panels[base_pair_N_percentage.legand] = base_pair_N_percentage^

        var per_sequence_cg_content = self.cg_content.make_html()
        panels[per_sequence_cg_content.legand] = per_sequence_cg_content^

        var sequence_length_distribution = self.len_dist.make_html()
        panels[sequence_length_distribution.legand] = sequence_length_distribution^

        var dup_html = self.dup_reads.make_html(Int(self.num_reads))
        var sequence_duplication_levels = dup_html[0].copy()
        var overrepresented_sequences = dup_html[1].copy()
        panels[sequence_duplication_levels.legand] = sequence_duplication_levels^
        panels[overrepresented_sequences.legand] = overrepresented_sequences^

        var adapter_content = self.adpt_cont.make_html(self.num_reads)
        panels[adapter_content.legand] = adapter_content^

        return panels^

    fn orchestrate(mut self, file_name: String) raises:
        self.write_data(file_name)  # write data file before any plotting
        var panels = self.build_panels()
        var html = render_html(panels.copy(), file_name)
        write_html_file(html, file_name)
        write_summary_file(panels, file_name)

    fn make_html(mut self, file_name: String) raises:
        self.orchestrate(file_name)


fn render_html(
    panels: Dict[String, result_panel], file_name: String
) raises -> String:
    var order = _panel_order()
    var html: String = html_template
    for name in order:
        html = insert_result_panel(html, panels[name])
    while html.find("<<filename>>") > -1:
        html = html.replace("<<filename>>", file_name)
    var py_dt = Python.import_module("datetime")
    var dt_now = py_dt.datetime.now().strftime("%a %d %b %Y")
    while html.find("<<date>>") > -1:
        html = html.replace("<<date>>", String(dt_now))
    return html


fn write_html_file(html: String, file_name: String) raises:
    with open("{}_blazeseq.html".format(file_name), "w") as f:
        print("{}_blazeseq.html".format(file_name))
        f.write(html)


fn write_summary_file(
    panels: Dict[String, result_panel], file_name: String
) raises:
    var order = _panel_order()
    with open("{}_blazeseq_summary.txt".format(file_name), "w") as f:
        for name in order:
            var line = _summary_line(panels[name].grade, name, file_name)
            f.write(line)
            f.write("\n")


fn _panel_order() raises -> List[String]:
    var order = List[String]()
    order.append("Basic Statistics")
    order.append("Per Sequence Quality Scores")
    order.append("Per Tile Sequence Quality")
    order.append("Per Base Sequence Quality")
    order.append("Per Base Sequence Content")
    order.append("Per Sequence GC Content")
    order.append("Per Base N Content")
    order.append("Sequence Length Distribution")
    order.append("Duplicate Sequences")
    order.append("Overrepresented Sequences")
    order.append("Adapter Content")
    return order^


fn _summary_line(grade: String, module_name: String, file_name: String) -> String:
    var status: String
    if grade == "pass":
        status = "PASS"
    elif grade == "warn":
        status = "WARN"
    else:
        status = "FAIL"
    return "{}\t{}\t{}".format(status, module_name, file_name)
