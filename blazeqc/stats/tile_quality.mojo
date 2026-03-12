"""Per-tile quality: Collector + Summarizer + TileQualityModule."""

from collections.dict import DictEntry, Dict, default_hasher
from collections.list import List
from python import Python, PythonObject
from blazeseq import FastqRecord, RefRecord
from blazeqc.stats.traits import Collector, Summarizer, PlotOutput
from blazeqc.stats.summary_utils import SummaryContext, GradeEntry, DefaultOutputter
from blazeqc.helpers import (
    Matrix2D,
    grow_tensor,
    make_linear_base_groups,
    matrix_to_numpy,
)
from blazeqc.html_maker import result_panel
from blazeqc.limits import TILE_WARN, TILE_ERROR

# Helper for tests: grade string from max deviation.
fn tile_quality_grade_from_deviation(max_deviation: Float64) -> String:
    if max_deviation > TILE_ERROR:
        return "fail"
    if max_deviation > TILE_WARN:
        return "warn"
    return "pass"


struct TileQualityEntry(Copyable, Movable):
    var tile: Int
    var count: Int
    var quality: List[Int64]

    fn __init__(out self, tile: Int, count: Int, length: Int):
        self.tile = tile
        self.count = count
        self.quality = List[Int64](capacity=length)
        for _ in range(length):
            self.quality.append(0)

    fn __hash__(self) -> UInt64:
        return hash(self.tile)

    fn __add__(self, other: Int) -> Int:
        return self.count + other

    fn __iadd__(mut self, other: Int):
        self.count += other


# ----- Collector: tally only -----

struct TileQualityCollector(Collector, Copyable, Movable):
    var n: Int
    var map: Dict[Int, TileQualityEntry]
    var max_length: Int
    var enabled: Bool

    fn __init__(out self):
        self.map = Dict[Int, TileQualityEntry](
            power_of_two_initial_capacity=2**14
        )
        self.n = 0
        self.max_length = 0
        self.enabled = True

    fn tally_read(mut self, record: FastqRecord):
        if not self.enabled:
            return
        self.n += 1
        if self.n >= 10_000:
            if self.n % 10 != 0:
                return
        var x = self._find_tile_info(record)
        if x == -1:
            self.enabled = False
            return
        var val = self._find_tile_value(record, x)
        index = self.map._find_index(hash(val), val)
        if index[0]:
            pos = index[2]
            entry = self.map._entries[pos]
            var deref_value = entry.unsafe_value().value.copy()
            deref_value.count += 1
            if len(deref_value.quality) < len(record):
                deref_value.quality = grow_tensor(
                    deref_value.quality, len(record)
                )
            var qu_span = record.quality()
            for i in range(len(record)):
                deref_value.quality[i] += Int(qu_span[i])
            self.map._entries[pos] = DictEntry[
                Int, TileQualityEntry, default_hasher
            ](val, deref_value^)
        else:
            if len(self.map) >= 2500:
                self.enabled = False
                return
            self.map[val] = TileQualityEntry(val, 1, len(record))
        if self.max_length < len(record):
            self.max_length = len(record)

    fn tally_read(mut self, record: RefRecord):
        if not self.enabled:
            return
        self.n += 1
        if self.n >= 10_000:
            if self.n % 10 != 0:
                return
        var x = self._find_tile_info(record)
        if x == -1:
            self.enabled = False
            return
        var val = self._find_tile_value(record, x)
        var index = self.map._find_index(hash(val), val)
        if index[0]:
            var pos = index[2]
            entry = self.map._entries[pos]
            var deref_value = entry.unsafe_value().value.copy()
            deref_value.count += 1
            if len(deref_value.quality) < len(record):
                deref_value.quality = grow_tensor(
                    deref_value.quality, len(record)
                )
            var qu_span = record.quality().as_bytes()
            for i in range(len(record)):
                deref_value.quality[i] += Int(qu_span[i])
            self.map._entries[pos] = DictEntry[
                Int, TileQualityEntry, default_hasher
            ](val, deref_value^)
        else:
            if len(self.map) >= 2500:
                self.enabled = False
                return
            self.map[val] = TileQualityEntry(val, 1, len(record))
        if self.max_length < len(record):
            self.max_length = len(record)

    @always_inline
    fn _find_tile_info(self, record: FastqRecord) -> Int:
        comptime sep: UInt8 = ord(":")
        var id_bytes = record.id()
        count = 0
        for i in range(len(id_bytes)):
            if id_bytes[i] == sep:
                count += 1
        var split_position: Int
        if count >= 6:
            split_position = 4
        elif count >= 4:
            split_position = 2
        else:
            return -1
        return split_position

    @always_inline
    fn _find_tile_value(self, record: FastqRecord, pos: Int) -> Int:
        comptime sep: UInt8 = ord(":")
        var index_1 = 0
        var index_2 = 0
        var count = 0
        var id_bytes = record.id()
        for i in range(len(id_bytes)):
            if id_bytes[i] == sep:
                count += 1
                if count == pos:
                    index_1 = i + 1
                if count == pos + 1:
                    index_2 = i
                    break
        var s = String(unsafe_from_utf8=record.id())
        try:
            return atol(s[index_1:index_2])
        except:
            return 0

    @always_inline
    fn _find_tile_info(self, record: RefRecord) -> Int:
        comptime sep: UInt8 = ord(":")
        var id_bytes = record.id().as_bytes()
        var count = 0
        for i in range(len(id_bytes)):
            if id_bytes[i] == sep:
                count += 1
        var split_position: Int
        if count >= 6:
            split_position = 4
        elif count >= 4:
            split_position = 2
        else:
            return -1
        return split_position

    @always_inline
    fn _find_tile_value(self, record: RefRecord, pos: Int) -> Int:
        comptime sep: UInt8 = ord(":")
        var index_1 = 0
        var index_2 = 0
        var count = 0
        var id_bytes = record.id().as_bytes()
        for i in range(len(id_bytes)):
            if id_bytes[i] == sep:
                count += 1
                if count == pos:
                    index_1 = i + 1
                if count == pos + 1:
                    index_2 = i
                    break
        var s = String(unsafe_from_utf8=record.id().as_bytes())
        try:
            return atol(s[index_1:index_2])
        except:
            return 0


# ----- Helpers (used by summarizer / prepare) -----

fn _sorted_tile_ids_from_map(map: Dict[Int, TileQualityEntry]) -> List[Int]:
    var ids = List[Int]()
    for k in map.keys():
        ids.append(k)
    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            if ids[i] > ids[j]:
                var tmp = ids[i]
                ids[i] = ids[j]
                ids[j] = tmp
    return ids^


fn _compute_group_means(
    map: Dict[Int, TileQualityEntry],
    max_length: Int,
    tile_ids: List[Int],
    groups: List[Int],
) raises -> Matrix2D[DType.float64]:
    var n_tiles = len(tile_ids)
    var n_groups = len(groups)
    var means = Matrix2D[DType.float64](n_tiles, n_groups)
    for t in range(n_tiles):
        ref entry = map[tile_ids[t]]
        var count = Float64(entry.count)
        for g in range(n_groups):
            var g_start = groups[g] - 1
            var g_end = groups[g + 1] - 1 if g + 1 < n_groups else max_length
            var qual_sum: Float64 = 0.0
            var width = g_end - g_start
            for p in range(g_start, g_end):
                if p < len(entry.quality):
                    qual_sum += Float64(entry.quality[p])
            if count > 0 and width > 0:
                means.set(t, g, qual_sum / (count * Float64(width)))
    return means^


fn _subtract_group_averages(
    mut means: Matrix2D[DType.float64],
    n_tiles: Int,
    n_groups: Int,
) -> Float64:
    for g in range(n_groups):
        var avg = means.col_sum(g) / Float64(n_tiles) if n_tiles > 0 else Float64(0.0)
        for t in range(n_tiles):
            means.set(t, g, means.get(t, g) - avg)
    var max_dev: Float64 = 0.0
    for t in range(n_tiles):
        for g in range(n_groups):
            var v = means.get(t, g)
            var absval = v if v >= 0.0 else -v
            if absval > max_dev:
                max_dev = absval
    return max_dev


fn _draw_heatmap(
    means: Matrix2D[DType.float64],
    tile_ids: List[Int],
    groups: List[Int],
    max_length: Int,
    max_deviation: Float64,
) raises -> PythonObject:
    var sns = Python.import_module("seaborn")
    var plt = Python.import_module("matplotlib.pyplot")
    var n_tiles = len(tile_ids)
    var n_groups = len(groups)
    var x_labels = Python.list()
    for g in range(n_groups):
        var start_bp = groups[g]
        var end_bp: Int
        if g + 1 < n_groups:
            end_bp = groups[g + 1] - 1
        else:
            end_bp = max_length
        var width = end_bp - start_bp + 1
        if width == 1:
            x_labels.append(String(start_bp))
        else:
            x_labels.append(String(start_bp) + "-" + String(end_bp))
    var y_labels = Python.list()
    for t in range(n_tiles):
        y_labels.append(tile_ids[t])
    var vmax = max_deviation if max_deviation > 0.0 else 1.0
    var z = plt.subplots(figsize=Python.tuple(10, max(4, n_tiles // 4 + 2)))
    var fig = z[0]
    var ax = z[1]
    sns.heatmap(
        matrix_to_numpy(means),
        cmap="RdYlBu",
        center=0.0,
        vmin=-vmax,
        vmax=vmax,
        yticklabels=y_labels,
        xticklabels=x_labels,
        ax=ax,
    )
    ax.set_title("Per tile quality")
    ax.set_xlabel("Position in read (bp)")
    ax.set_ylabel("Tile")
    return fig


# ----- Prepared data -----

struct TileQualityPreparedData(Copyable, Movable):
    var tile_ids: List[Int]
    var groups: List[Int]
    var means: Matrix2D[DType.float64]
    var max_length: Int
    var max_deviation: Float64

    fn __init__(out self):
        self.tile_ids = List[Int]()
        self.groups = List[Int]()
        self.means = Matrix2D[DType.float64](0, 0)
        self.max_length = 0
        self.max_deviation = 0.0

    fn __copyinit__(out self, other: Self):
        self.tile_ids = other.tile_ids.copy()
        self.groups = other.groups.copy()
        self.means = other.means
        self.max_length = other.max_length
        self.max_deviation = other.max_deviation


fn _prepare_tile_quality(collector: TileQualityCollector) raises -> TileQualityPreparedData:
    var tile_ids = _sorted_tile_ids_from_map(collector.map)
    var groups = make_linear_base_groups(collector.max_length)
    var means = _compute_group_means(
        collector.map, collector.max_length, tile_ids, groups
    )
    var n_tiles = len(tile_ids)
    var n_groups = len(groups)
    var max_dev = _subtract_group_averages(means, n_tiles, n_groups)
    var result = TileQualityPreparedData()
    result.tile_ids = tile_ids^
    result.groups = groups^
    result.means = means^
    result.max_length = collector.max_length
    result.max_deviation = max_dev
    return result^


# ----- Summarizer -----

struct TileQualitySummarizer(Summarizer, PlotOutput, Copyable, Movable):
    var _prepared: TileQualityPreparedData
    var _cache_ready: Bool

    fn __init__(out self):
        self._prepared = TileQualityPreparedData()
        self._cache_ready = False

    fn feed_prepared(mut self, prepared: TileQualityPreparedData):
        self._prepared = prepared.copy()

    fn summerize(mut self, ctx: SummaryContext) raises:
        self._cache_ready = True

    fn grade(self) raises -> GradeEntry:
        var status = tile_quality_grade_from_deviation(self._prepared.max_deviation)
        return GradeEntry("Per Tile Sequence Quality", status)

    fn data_block_body(self) raises -> String:
        if not self._cache_ready:
            return ""
        var out = "#Tile\tBase\tMean\n"
        var n_tiles = len(self._prepared.tile_ids)
        var n_groups = len(self._prepared.groups)
        for t in range(n_tiles):
            var tile_id = self._prepared.tile_ids[t]
            for g in range(n_groups):
                var start_bp = self._prepared.groups[g]
                var end_bp: Int
                if g + 1 < n_groups:
                    end_bp = self._prepared.groups[g + 1] - 1
                else:
                    end_bp = self._prepared.max_length
                var base_label: String
                if end_bp - start_bp + 1 == 1:
                    base_label = String(start_bp)
                else:
                    base_label = "{}-{}".format(start_bp, end_bp)
                var mean_val = self._prepared.means.get(t, g)
                out += "{}\t{}\t{}\n".format(tile_id, base_label, mean_val)
        return out

    fn module_legend(self) -> String:
        return "Per Tile Sequence Quality"

    fn panel_id(self) -> String:
        return "tile_quality"

    fn plot_result(self) raises -> PythonObject:
        return _draw_heatmap(
            self._prepared.means,
            self._prepared.tile_ids,
            self._prepared.groups,
            self._prepared.max_length,
            self._prepared.max_deviation,
        )


# ----- Assembled module -----

struct TileQualityModule(Copyable, Movable):
    var collector: TileQualityCollector
    var summarizer: TileQualitySummarizer

    fn __init__(out self):
        self.collector = TileQualityCollector()
        self.summarizer = TileQualitySummarizer()

    fn tally_read(mut self, record: FastqRecord):
        self.collector.tally_read(record)

    fn tally_read(mut self, record: RefRecord):
        self.collector.tally_read(record)

    fn prepare_summarizers(mut self, ctx: SummaryContext) raises:
        var prepared = _prepare_tile_quality(self.collector)
        self.summarizer.feed_prepared(prepared)
        self.summarizer.summerize(ctx)

    fn to_data_text(self, ctx: SummaryContext) raises -> String:
        var out = DefaultOutputter()
        var body = self.summarizer.data_block_body()
        var g = self.summarizer.grade()
        return out.wrap_data_block(
            self.summarizer.module_legend(), g.grade, body
        )

    fn to_html_panels(self) raises -> List[result_panel]:
        var panels = List[result_panel]()
        var out = DefaultOutputter()
        var fig = self.summarizer.plot_result()
        panels.append(
            out.make_panel(
                self.summarizer.panel_id(),
                self.summarizer.grade().grade,
                self.summarizer.module_legend(),
                fig,
            )
        )
        return panels^

    fn plot_result(self) raises -> PythonObject:
        return self.summarizer.plot_result()
