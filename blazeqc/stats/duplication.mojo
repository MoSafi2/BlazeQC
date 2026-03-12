"""Duplication and over-represented sequences: Collector + two Summarizers + DupModule."""

from collections.dict import Dict
from collections.list import List
from python import Python, PythonObject
from blazeseq import FastqRecord, RefRecord
from blazeqc.stats.traits import Collector, Summarizer, PlotOutput
from blazeqc.stats.reporting_traits import FastqcDataOutput, FastqcHtmlOutput
from blazeqc.stats.summary_utils import SummaryContext, GradeEntry, DefaultOutputter, DataEntry, PanelEntry
from blazeqc.stats.over_represented import OverRepresentedSequence
from blazeqc.helpers import list_float64_to_numpy, encode_img_b64
from blazeqc.html_maker import result_panel, _make_row, _make_table
from blazeqc.limits import DUPLICATION_WARN, DUPLICATION_ERROR


comptime MAX_READS = 100_000

# ----- Helper: correct_values (from original DupReads) -----

fn _correct_values(
    dup_level: Int, count_at_level: Int, count_at_max: Int, total_count: Int
) -> Float64:
    if count_at_max == total_count:
        return count_at_level
    if total_count - count_at_level < count_at_max:
        return count_at_level
    var pNotSeeingAtLimit: Float64 = 1
    var limitOfCaring = Float64(1) - (
        count_at_level / (count_at_level + 0.01)
    )
    for i in range(count_at_max):
        pNotSeeingAtLimit *= ((total_count - i) - dup_level) / (
            total_count - i
        )
        if pNotSeeingAtLimit < limitOfCaring:
            pNotSeeingAtLimit = 0
            break
    var pSeeingAtLimit: Float64 = 1 - pNotSeeingAtLimit
    return count_at_level / pSeeingAtLimit


fn _dup_slot(dup_level: Int) -> Int:
    var dup_slot = min(max(dup_level - 1, 0), 15)
    if dup_slot > 9999 or dup_slot < 0:
        return 15
    if dup_slot > 4999:
        return 14
    if dup_slot > 999:
        return 13
    if dup_slot > 499:
        return 12
    if dup_slot > 99:
        return 11
    if dup_slot > 49:
        return 10
    if dup_slot > 9:
        return 9
    return dup_slot


fn cmp_over_repr(
    a: OverRepresentedSequence,
    b: OverRepresentedSequence,
) capturing -> Bool:
    if a.percentage > b.percentage:
        return True
    elif a.percentage < b.percentage:
        return False
    return False


# ----- Collector: tally only -----

struct DupCollector(Collector, Copyable, Movable):
    var unique_dict: Dict[String, Int]
    var unique_reads: Int
    var count_at_max: Int
    var n: Int

    fn __init__(out self):
        self.unique_dict = Dict[String, Int](
            power_of_two_initial_capacity=2**18
        )
        self.unique_reads = 0
        self.count_at_max = 0
        self.n = 0

    fn tally_read(mut self, record: FastqRecord):
        self.n += 1
        var read_len = min(len(record), 50)
        var s: String
        try:
            s = String(unsafe_from_utf8=record.sequence()[0:read_len])
        except:
            s = ""
        if s in self.unique_dict:
            try:
                self.unique_dict[s] += 1
                return
            except error:
                print(error)
                pass
        if self.unique_reads <= MAX_READS:
            self.unique_dict[s] = 1
            self.unique_reads += 1
            if self.unique_reads <= MAX_READS:
                self.count_at_max = self.n
        else:
            return

    fn tally_read(mut self, record: RefRecord):
        self.n += 1
        var read_len = min(len(record), 50)
        var s: String
        try:
            s = String(unsafe_from_utf8=record.sequence().as_bytes()[0:read_len])
        except:
            s = ""
        if s in self.unique_dict:
            try:
                self.unique_dict[s] += 1
                return
            except error:
                print(error)
                pass
        if self.unique_reads <= MAX_READS:
            self.unique_dict[s] = 1
            self.unique_reads += 1
            if self.unique_reads <= MAX_READS:
                self.count_at_max = self.n
        else:
            return


# ----- Prepared data (computed once from collector) -----

fn _compute_corrected_counts(collector: DupCollector) -> Dict[Int, Float64]:
    var dup_dict = Dict[Int, Int]()
    for entry in collector.unique_dict.items():
        if Int(entry.value) in dup_dict:
            try:
                dup_dict[Int(entry.value)] += 1
            except error:
                print(error)
        else:
            dup_dict[Int(entry.value)] = 1
    var corrected_reads = Dict[Int, Float64]()
    for entry in dup_dict:
        try:
            var level = entry
            var count = dup_dict[level]
            var corrected_count = _correct_values(
                level, count, collector.count_at_max, collector.n
            )
            corrected_reads[level] = corrected_count
        except:
            print("Error")
    return corrected_reads^


fn _percent_remaining_after_dedup(
    corrected_counts: Dict[Int, Float64], n: Int
) -> Float64:
    var dedup_total: Float64 = 0
    var raw_total: Float64 = 0
    for entry in corrected_counts.items():
        dedup_total += entry.value
        raw_total += entry.value * Float64(entry.key)
    if raw_total <= 0:
        return 100.0
    return (dedup_total / raw_total) * 100.0


struct DupPreparedData(Copyable, Movable):
    var dup_percentages: List[Float64]
    var dup_grade: String
    var overrepresented: List[OverRepresentedSequence]

    fn __init__(out self):
        self.dup_percentages = List[Float64]()
        self.dup_grade = ""
        self.overrepresented = List[OverRepresentedSequence]()

fn _prepare_dup_data(
    collector: DupCollector, total_reads: Int
) -> DupPreparedData:
    var corrected_counts = _compute_corrected_counts(collector)
    var total_percentages = List[Float64](capacity=16)
    for _ in range(16):
        total_percentages.append(0)
    for entry in corrected_counts.items():
        var count = entry.value
        var dup_level = entry.key
        total_percentages[_dup_slot(dup_level)] += count * dup_level
    var result = DupPreparedData()
    for i in range(16):
        result.dup_percentages.append(
            (total_percentages[i] / Float64(total_reads)) * Float64(100)
        )
    var pct = _percent_remaining_after_dedup(corrected_counts, collector.n)
    if pct < DUPLICATION_ERROR:
        result.dup_grade = "fail"
    elif pct < DUPLICATION_WARN:
        result.dup_grade = "warn"
    else:
        result.dup_grade = "pass"

    for key in collector.unique_dict.items():
        var seq_pct = (Float64(key.value) / Float64(collector.n)) * 100.0
        if seq_pct > 0.1:
            result.overrepresented.append(
                OverRepresentedSequence(
                    String(key.key), key.value, seq_pct, String("No Hit")
                )
            )
    sort[cmp_fn=cmp_over_repr](result.overrepresented)
    return result^


# ----- Summarizer: Duplicate Sequences -----

struct DuplicateSequencesSummarizer(Summarizer, PlotOutput, Copyable, Movable):
    var _cache_dup_percentages: List[Float64]
    var _cache_grade: String
    var _cache_ready: Bool

    fn __init__(out self):
        self._cache_dup_percentages = List[Float64]()
        self._cache_grade = ""
        self._cache_ready = False

    fn feed_prepared(mut self, dup_percentages: List[Float64], grade: String):
        self._cache_dup_percentages = dup_percentages.copy()
        self._cache_grade = grade
        self._cache_ready = True

    fn summerize(mut self, ctx: SummaryContext) raises:
        pass

    fn grade(self) raises -> GradeEntry:
        return GradeEntry("Duplicate Sequences", self._cache_grade)

    fn data_block_body(self) raises -> String:
        if not self._cache_ready:
            return ""
        var out = "#Duplication Level\tPercentage of total\n"
        var tick_labels = List[String]()
        tick_labels.append("1"); tick_labels.append("2"); tick_labels.append("3"); tick_labels.append("4"); tick_labels.append("5"); tick_labels.append("6"); tick_labels.append("7"); tick_labels.append("8"); tick_labels.append("9")
        tick_labels.append(">10"); tick_labels.append(">50"); tick_labels.append(">100"); tick_labels.append(">500"); tick_labels.append(">1k"); tick_labels.append(">5k"); tick_labels.append(">10k+")
        for i in range(len(self._cache_dup_percentages)):
            out += "{}\t{}\n".format(tick_labels[i], self._cache_dup_percentages[i])
        return out

    fn module_legend(self) -> String:
        return "Duplicate Sequences"

    fn panel_id(self) -> String:
        return "dup_reads"

    fn plot_result(self) raises -> PythonObject:
        var plt = Python.import_module("matplotlib.pyplot")
        var final_arr = list_float64_to_numpy(self._cache_dup_percentages)
        var f = plt.subplots()
        var fig = f[0]
        var ax = f[1]
        ax.plot(final_arr)
        var tick_positions = Python.list(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
        ax.set_xticks(tick_positions)
        var tick_labels = Python.list(
            "1", "2", "3", "4", "5", "6", "7", "8", "9",
            ">10", ">50", ">100", ">500", ">1k", ">5k", ">10k+",
        )
        ax.set_xticklabels(tick_labels)
        ax.set_xlabel("Sequence Duplication Level")
        ax.set_title("Sequence duplication levels")
        return fig


# ----- Summarizer: Overrepresented Sequences (table panel) -----

struct OverrepresentedSequencesSummarizer(Summarizer, PlotOutput, Copyable, Movable):
    var _cache_overrepresented: List[OverRepresentedSequence]
    var _cache_ready: Bool

    fn __init__(out self):
        self._cache_overrepresented = List[OverRepresentedSequence]()
        self._cache_ready = False

    fn feed_prepared(mut self, overrepresented: List[OverRepresentedSequence]):
        self._cache_overrepresented = overrepresented.copy()
        self._cache_ready = True

    fn summerize(mut self, ctx: SummaryContext) raises:
        pass

    fn grade(self) raises -> GradeEntry:
        return GradeEntry("Overrepresented Sequences", "pass")

    fn data_block_body(self) raises -> String:
        if not self._cache_ready:
            return ""
        var out = "#Sequence\tCount\tPercentage\tPossible Source\n"
        for entry in self._cache_overrepresented:
            out += "{}\t{}\t{}\t{}\n".format(entry.seq, entry.count, entry.percentage, entry.hit)
        return out

    fn module_legend(self) -> String:
        return "Overrepresented Sequences"

    fn panel_id(self) -> String:
        return "over_represented_seqs"

    fn plot_result(self) raises -> PythonObject:
        return Python.evaluate("None")

    fn table_html(self) raises -> String:
        var rows: String = ""
        for entry in self._cache_overrepresented:
            rows += _make_row(
                entry.seq, entry.count, entry.percentage, entry.hit
            )
        return _make_table(rows)


# ----- Assembled module -----

struct DupModule(FastqcDataOutput, FastqcHtmlOutput, Copyable, Movable):
    var collector: DupCollector
    var summarizer_dup: DuplicateSequencesSummarizer
    var summarizer_overrepr: OverrepresentedSequencesSummarizer

    fn __init__(out self):
        self.collector = DupCollector()
        self.summarizer_dup = DuplicateSequencesSummarizer()
        self.summarizer_overrepr = OverrepresentedSequencesSummarizer()

    fn tally_read(mut self, record: FastqRecord):
        self.collector.tally_read(record)

    fn tally_read(mut self, record: RefRecord):
        self.collector.tally_read(record)

    fn prepare_summarizers(mut self, ctx: SummaryContext) raises:
        var total_reads = Int(ctx.num_reads)
        var prepared = _prepare_dup_data(self.collector, total_reads)
        self.summarizer_dup.feed_prepared(
            prepared.dup_percentages, prepared.dup_grade
        )
        self.summarizer_dup.summerize(ctx)
        self.summarizer_overrepr.feed_prepared(prepared.overrepresented)
        self.summarizer_overrepr.summerize(ctx)

    fn data_entries(self, ctx: SummaryContext) raises -> List[DataEntry]:
        var body_dup = self.summarizer_dup.data_block_body()
        var g_dup = self.summarizer_dup.grade()
        var body_over = self.summarizer_overrepr.data_block_body()
        var g_over = self.summarizer_overrepr.grade()
        var entries = List[DataEntry]()
        entries.append(DataEntry(self.summarizer_dup.module_legend(), g_dup.grade, body_dup))
        entries.append(DataEntry("Overrepresented sequences", g_over.grade, body_over))
        return entries^

    fn panel_entries(self, figures: List[PythonObject]) raises -> List[PanelEntry]:
        var table_html = self.summarizer_overrepr.table_html()
        var py_none = Python.evaluate("None")
        var entries = List[PanelEntry]()
        entries.append(PanelEntry(
            self.summarizer_dup.panel_id(),
            self.summarizer_dup.grade().grade,
            self.summarizer_dup.module_legend(),
            "image",
            figures[0],
            "",
        ))
        entries.append(PanelEntry(
            self.summarizer_overrepr.panel_id(),
            self.summarizer_overrepr.grade().grade,
            self.summarizer_overrepr.module_legend(),
            "table",
            py_none,
            table_html,
        ))
        return entries^
