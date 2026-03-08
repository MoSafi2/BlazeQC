"""Per-tile quality (split from stats_.mojo)."""

from collections.dict import DictEntry, Dict, default_hasher
from python import Python, PythonObject
from blazeseq import FastqRecord, RefRecord
from blazeqc.stats.analyser import Analyser
from blazeqc.helpers import (
    Matrix2D,
    encode_img_b64,
    grow_tensor,
    make_linear_base_groups,
    matrix_to_numpy,
)
from blazeqc.html_maker import result_panel


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


struct PerTileQuality(Analyser, Copyable, Movable):
    var n: Int
    var map: Dict[Int, TileQualityEntry]
    var max_length: Int
    var enabled: Bool
    var max_deviation: Float64

    fn __init__(out self):
        self.map = Dict[Int, TileQualityEntry](
            power_of_two_initial_capacity=2**14
        )
        self.n = 0
        self.max_length = 0
        self.enabled = True
        self.max_deviation = 0.0

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

        # Low-level access to the hashmap to avoid the overhead of calling `_find_index` multiple times.
        # Should be replaced with a cleaner version once Mojo dict is more performant.

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

            for i in range(len(record)):
                deref_value.quality[i] += Int(record.quality[i])

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

            for i in range(len(record)):
                deref_value.quality[i] += Int(record.quality[i])

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

    fn plot(mut self) raises -> PythonObject:
        """Produce a per-tile quality deviation heatmap matching FastQC output.

        Algorithm (mirrors FastQC PerTileQualityScores.getPercentages):
          1. Sort tile IDs and bin read positions into exponential groups
             (e.g. 1,2,...,9,10-14,15-19,...) via make_linear_base_groups().
          2. For each tile, compute the mean quality per position group:
               mean[t][g] = sum(quality[g_start:g_end]) / (count * group_width)
          3. Compute the cross-tile average per group:
               avg[g] = mean over all tiles of mean[t][g]
          4. Subtract to get deviations:
               deviation[t][g] = mean[t][g] - avg[g]
          5. Plot a diverging RdYlBu heatmap centered at 0, scaled to
             ±maxDeviation. Blue = tile worse than average; red = better.
        """
        var tile_ids = self._sorted_tile_ids()
        var n_tiles = len(tile_ids)
        var groups = make_linear_base_groups(self.max_length)

        var means = self._compute_group_means(tile_ids, groups)
        var n_groups = len(groups)
        self.max_deviation = self._subtract_group_averages(means, n_tiles, n_groups)

        return self._draw_heatmap(means, tile_ids, groups)

    fn _sorted_tile_ids(self) -> List[Int]:
        """Return tile IDs from the map sorted in ascending order."""
        var ids = List[Int]()
        for k in self.map.keys():
            ids.append(k)
        for i in range(len(ids)):
            for j in range(i + 1, len(ids)):
                if ids[i] > ids[j]:
                    var tmp = ids[i]
                    ids[i] = ids[j]
                    ids[j] = tmp
        return ids^

    fn _compute_group_means(
        self, tile_ids: List[Int], groups: List[Int]
    ) raises -> Matrix2D[DType.float64]:
        """Compute mean quality per tile per position group.

        Returns a [n_tiles x n_groups] matrix where:
          means[t][g] = sum(raw_quality_bytes in group) / (count * group_width)
        Raw quality bytes are Phred+offset sums; the offset cancels in the
        subsequent deviation subtraction so no explicit decoding is needed.
        """
        var n_tiles = len(tile_ids)
        var n_groups = len(groups)
        var means = Matrix2D[DType.float64](n_tiles, n_groups)
        for t in range(n_tiles):
            ref entry = self.map[tile_ids[t]]
            var count = Float64(entry.count)
            for g in range(n_groups):
                var g_start = groups[g] - 1
                var g_end = groups[g + 1] - 1 if g + 1 < n_groups else self.max_length
                var qual_sum: Float64 = 0.0
                var width = g_end - g_start
                for p in range(g_start, g_end):
                    if p < len(entry.quality):
                        qual_sum += Float64(entry.quality[p])
                if count > 0 and width > 0:
                    means.set(t, g, qual_sum / (count * Float64(width)))
        return means^

    fn _subtract_group_averages(
        mut self,
        mut means: Matrix2D[DType.float64],
        n_tiles: Int,
        n_groups: Int,
    ) -> Float64:
        """Subtract per-group cross-tile average from means in-place.

        Returns the maximum absolute deviation (used for status and colormap scaling).
        """
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
        self,
        means: Matrix2D[DType.float64],
        tile_ids: List[Int],
        groups: List[Int],
    ) raises -> PythonObject:
        """Render the deviation matrix as a seaborn heatmap and return the figure."""
        var sns = Python.import_module("seaborn")
        var plt = Python.import_module("matplotlib.pyplot")

        var n_tiles = len(tile_ids)
        var n_groups = len(groups)

        var x_labels = Python.list()
        for g in range(n_groups):
            x_labels.append(groups[g])

        var y_labels = Python.list()
        for t in range(n_tiles):
            y_labels.append(tile_ids[t])

        var vmax = self.max_deviation if self.max_deviation > 0.0 else 1.0
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
        ax.set_title("Quality per tile")
        ax.set_xlabel("Position in read (bp)")
        ax.set_ylabel("Tile")
        return fig

    fn _get_status(self) -> String:
        if self.max_deviation > 5.0:
            return "fail"
        if self.max_deviation > 2.0:
            return "warn"
        return "pass"

    fn make_html(mut self) raises -> result_panel:
        fig1 = self.plot()
        var encoded_fig1 = encode_img_b64(fig1)
        var result_1 = result_panel(
            "tile_quality",
            self._get_status(),
            "Quality per tile",
            encoded_fig1,
        )
        return result_1^

    # TODO: This function should return tile information
    @always_inline
    fn _find_tile_info(self, record: FastqRecord) -> Int:
        comptime sep: UInt8 = ord(":")
        count = 0
        for i in range(len(record.id)):
            if record.id[i] == sep:
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
        var id_slice = record.id_slice()
        # TODO: Add Error Handling
        for i in range(len(record.id)):
            if record.id[i] == sep:
                count += 1
                if count == pos:
                    index_1 = i + 1
                if count == pos + 1:
                    index_2 = i
                    break

        var s = String(id_slice)

        try:
            return atol(s[index_1:index_2])
        except:
            return 0

    @always_inline
    fn _find_tile_info(self, record: RefRecord) -> Int:
        comptime sep: UInt8 = ord(":")
        var count = 0
        for i in range(len(record.id)):
            if record.id[i] == sep:
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
        var id_slice = record.id_slice()
        for i in range(len(record.id)):
            if record.id[i] == sep:
                count += 1
                if count == pos:
                    index_1 = i + 1
                if count == pos + 1:
                    index_2 = i
                    break

        var s = String(id_slice)

        try:
            return atol(s[index_1:index_2])
        except:
            return 0
