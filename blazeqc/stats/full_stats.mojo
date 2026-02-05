"""FullStats aggregate (split from stats_.mojo)."""

from python import Python, PythonObject
from blazeseq import FastqRecord, RecordCoord
from blazeqc.stats.analyser import py_lib
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


@value
struct FullStats(Movable & Copyable):
    var num_reads: Int64
    var total_bases: Int64
    var bp_dist: BasepairDistribution
    var len_dist: LengthDistribution
    var qu_dist: QualityDistribution
    var cg_content: CGContent
    var dup_reads: DupReads
    var tile_qual: PerTileQuality
    var adpt_cont: AdapterContent

    fn __init__(out self):
        self.num_reads = 0
        self.total_bases = 0
        self.len_dist = LengthDistribution()
        self.bp_dist = BasepairDistribution()
        self.cg_content = CGContent()
        self.qu_dist = QualityDistribution()
        self.dup_reads = DupReads()
        self.tile_qual = PerTileQuality()
        self.adpt_cont = AdapterContent(hash_list(), 12)

    @always_inline
    fn tally(mut self, record: FastqRecord):
        self.num_reads += 1
        self.total_bases += record.len_record()
        self.bp_dist.tally_read(record)
        self.len_dist.tally_read(record)
        self.cg_content.tally_read(record)  # Almost Free
        self.dup_reads.tally_read(record)
        self.qu_dist.tally_read(record)
        self.adpt_cont.tally_read(record, self.num_reads)
        self.tile_qual.tally_read(record)

    @always_inline
    fn tally(mut self, record: RecordCoord):
        self.num_reads += 1
        self.total_bases += Int(record.seq_len())
        self.bp_dist.tally_read(record)
        self.len_dist.tally_read(record)
        self.cg_content.tally_read(record)
        self.qu_dist.tally_read(record)
        pass

    @always_inline
    fn make_base_stats(self) raises -> result_panel:
        var sum: Int64 = 0
        for i in range(len(self.cg_content.cg_content)):
            sum += self.cg_content.cg_content[i] * i
        var avg_cg = (sum / self.num_reads)
        var schema = self.qu_dist._guess_schema()

        total_bases = format_length(Float64(self.total_bases))

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
            "Base Statistics",
            table_template,
            panel_type="table",
        )
        return res

    @always_inline
    fn plot(mut self) raises -> List[PythonObject]:
        plots = List[PythonObject]()

        img1, img2 = self.bp_dist.plot(self.num_reads)
        plots.append(img1)
        plots.append(img2)
        plots.append(self.cg_content.plot())
        plots.append(self.len_dist.plot())
        img, _ = self.dup_reads.plot(Int(self.num_reads))
        plots.append(img)
        img1, img2 = self.qu_dist.plot()
        plots.append(img1)
        plots.append(img2)
        plots.append(self.tile_qual.plot())
        plots.append(self.adpt_cont.plot(self.num_reads))

        return plots

    fn make_html(mut self, file_name: String) raises:
        Python.add_to_path(py_lib.as_string_slice())

        py_dt = Python.import_module("datetime")
        dt_now = py_dt.datetime.now().strftime("%a %d %b %Y")

        var results = List[result_panel]()
        results.append(self.make_base_stats())
        res1, res2 = self.qu_dist.make_html()
        results.append(res2)
        results.append(res1)

        n_dist, bp_dist = self.bp_dist.make_html(self.num_reads)

        results.append(bp_dist)
        results.append(self.cg_content.make_html())
        results.append(n_dist)
        results.append(self.len_dist.make_html())
        dup_reads, over_represented = self.dup_reads.make_html(
            Int(self.num_reads)
        )
        results.append(dup_reads)
        results.append(over_represented)
        results.append(self.tile_qual.make_html())
        results.append(self.adpt_cont.make_html(self.num_reads))

        var html: String = html_template
        for entry in results:
            html = insert_result_panel(html, entry[])

        while html.find("<<filename>>") > -1:
            html = html.replace("<<filename>>", file_name)

        while html.find("<<date>>") > -1:
            html = html.replace("<<date>>", String(dt_now))

        with open("{}_blazeseq.html".format(file_name), "w") as f:
            print("{}_blazeseq.html".format(file_name))
            f.write(html)
