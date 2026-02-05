"""Duplication and over-represented sequences (split from stats_.mojo)."""

from collections.dict import Dict
from python import Python, PythonObject
from blazeseq import FastqRecord
from blazeqc.stats.analyser import Analyser
from blazeqc.stats.over_represented import OverRepresentedSequence
from blazeqc.helpers import list_float64_to_numpy, encode_img_b64
from blazeqc.html_maker import result_panel, _make_row, _make_table


struct DupReads(Analyser, Copyable, Movable):
    var unique_dict: Dict[String, Int]
    var unique_reads: Int
    var count_at_max: Int
    var n: Int
    var corrected_counts: Dict[Int, Float64]
    comptime MAX_READS = 100_000

    fn __init__(out self):
        self.unique_dict = Dict[String, Int](
            power_of_two_initial_capacity=2**18
        )
        self.unique_reads = 0
        self.count_at_max = 0
        self.n = 0
        self.corrected_counts = Dict[Int, Float64]()

    # TODO: Check if the Stringslice to String Conversion is right
    fn tally_read(mut self, record: FastqRecord):
        self.n += 1
        var read_len = min(len(record), 50)
        var s: String
        try:
            s = String(record.get_seq()[0:read_len])
        except:
            s = ""
        if s in self.unique_dict:
            try:
                self.unique_dict[s] += 1
                return
            except error:
                print(error)
                pass

        if self.unique_reads <= self.MAX_READS:
            self.unique_dict[s] = 1
            self.unique_reads += 1

            if self.unique_reads <= self.MAX_READS:
                self.count_at_max = self.n
        else:
            return

    fn predict_reads(mut self):
        # Construct Duplication levels dict
        var dup_dict = Dict[Int, Int]()
        for entry in self.unique_dict.items():
            if Int(entry.value) in dup_dict:
                try:
                    dup_dict[Int(entry.value)] += 1
                except error:
                    print(error)
            else:
                dup_dict[Int(entry.value)] = 1

        # Correct reads levels
        var corrected_reads = Dict[Int, Float64]()
        for entry in dup_dict:
            try:
                var level = entry
                var count = dup_dict[level]
                var corrected_count = self.correct_values(
                    level, count, self.count_at_max, self.n
                )
                corrected_reads[level] = corrected_count
            except:
                print("Error")

        self.corrected_counts = corrected_reads^

    @staticmethod
    fn correct_values(
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
        var trueCount = count_at_level / pSeeingAtLimit
        return trueCount

    fn plot(
        mut self, total_reads: Int
    ) raises -> Tuple[PythonObject, List[OverRepresentedSequence]]:
        ###################################################################
        ###                     Duplicate Reads                         ###
        ###################################################################
        # Correct if we didn't hit the number of unique reads, thus count_at_max stays at 0:
        self.predict_reads()
        var total_percentages = List[Float64](capacity=16)
        for _ in range(16):
            total_percentages.append(0)
        var dedup_total: Float64 = 0
        var raw_total: Float64 = 0
        for entry in self.corrected_counts.items():
            var count = entry.value
            var dup_level = entry.key
            dedup_total += count
            raw_total += count * dup_level
            var dup_slot = min(max(dup_level - 1, 0), 15)
            # Handle edge cases for duplication levels
            if dup_slot > 9999 or dup_slot < 0:
                dup_slot = 15
            elif dup_slot > 4999:
                dup_slot = 14
            elif dup_slot > 999:
                dup_slot = 13
            elif dup_slot > 499:
                dup_slot = 12
            elif dup_slot > 99:
                dup_slot = 11
            elif dup_slot > 49:
                dup_slot = 10
            elif dup_slot > 9:
                dup_slot = 9
            total_percentages[dup_slot] += count * dup_level

        var plt = Python.import_module("matplotlib.pyplot")
        var new_arr = List[Float64](capacity=len(total_percentages))
        for i in range(len(total_percentages)):
            new_arr.append(
                (total_percentages[i] / Float64(total_reads)) * Float64(100)
            )
        var final_arr = list_float64_to_numpy(new_arr)
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
        ax.set_title("Sequences Duplication levels")
        # ax.set_ylim(0, 100)

        ################################################################
        ####               Over-Represented Sequences                ###
        ################################################################

        # TODO: Check also those over-representing stuff against the contaimination list.
        var overrepresented_seqs = List[OverRepresentedSequence]()
        for key in self.unique_dict.items():
            var seq_precent = (Float64(key.value) / Float64(self.n)) * 100.0
            if seq_precent > 0.1:
                overrepresented_seqs.append(
                    OverRepresentedSequence(
                        String(key.key), key.value, seq_precent, String("No Hit")
                    )
                )

        sort[cmp_fn=cmp_over_repr](overrepresented_seqs)

        return (fig^, overrepresented_seqs^)

    fn make_html(
        mut self, total_reads: Int
    ) raises -> Tuple[result_panel, result_panel]:
        var plot_result = self.plot(total_reads)
        var fig = plot_result[0]
        var encoded_fig1 = encode_img_b64(fig)
        var result_1 = result_panel(
            "dup_reads",
            "pass",
            "Sequence Duplication Levels",
            encoded_fig1,
        )

        var rows: String = ""
        for entry in plot_result[1]:
            var row = _make_row(
                entry.seq, entry.count, entry.percentage, entry.hit
            )
            rows += row
        var over_repr_table = _make_table(rows)

        var result_2 = result_panel(
            "over_represented_seqs",
            "pass",
            "Over-represented Sequences",
            over_repr_table,
            panel_type="table",
        )

        return (result_1^, result_2^)


fn cmp_over_repr(
    a: OverRepresentedSequence,
    b: OverRepresentedSequence,
) capturing -> Bool:
    if a.percentage > b.percentage:
        return True
    elif a.percentage < b.percentage:
        return False
    else:
        return False
