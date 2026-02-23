"""Per-tile quality (split from stats_.mojo)."""

from collections.dict import DictEntry, Dict, default_hasher
from python import Python, PythonObject
from blazeseq import FastqRecord
from blazeqc.stats.analyser import Analyser
from blazeqc.helpers import tensor_to_numpy_1d, encode_img_b64, grow_tensor
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

    fn __init__(out self):
        self.map = Dict[Int, TileQualityEntry](
            power_of_two_initial_capacity=2**14
        )
        self.n = 0
        self.max_length = 0

    # TODO: Add tracking for the number of items inside the hashmaps to limit it to 2_500 items.
    fn tally_read(mut self, record: FastqRecord):
        self.n += 1
        if self.n >= 10_000:
            if self.n % 10 != 0:
                return

        var x = self._find_tile_info(record)
        var val = self._find_tile_value(record, x)

        # Low-level access to the hashmap to avoid the overhead of calling `_find_index` multiple times.
        # Should be replcaed with a cleaner version once Mojo dict is more performant.

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
            self.map[val] = TileQualityEntry(val, 1, len(record))

        if self.max_length < len(record):
            self.max_length = len(record)

    # TODO: Construct a n_keys*max_length array to hold all information.
    fn plot(self) raises -> PythonObject:
        var np = Python.import_module("numpy")
        var sns = Python.import_module("seaborn")
        var plt = Python.import_module("matplotlib.pyplot")

        var z = plt.subplots()
        var fig = z[0]
        var ax = z[1]

        var arr = np.zeros(self.max_length)
        for i in self.map.keys():
            temp_arr = tensor_to_numpy_1d(self.map[i].quality)
            arr = np.vstack(Python.tuple(arr, temp_arr))

        var ks = Python.list()
        for i in self.map.keys():
            ks.append(i)
        sns.heatmap(arr[1:,], cmap="Blues_r", yticklabels=ks, ax=ax)
        ax.set_title("Quality per tile")
        ax.set_xlabel("Position in read (bp)")

        return fig

    fn make_html(self) raises -> result_panel:
        fig1 = self.plot()
        var encoded_fig1 = encode_img_b64(fig1)
        var result_1 = result_panel(
            "tile_quality",
            "pass",
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
