"""K-mer content (split from stats_.mojo)."""

from utils import Index
from memory import Span
from python import Python, PythonObject
from blazeseq import FastqRecord
from blazeqc.stats.analyser import Analyser
from blazeqc.helpers import base2int, grow_tensor, Matrix2D, matrix_to_numpy


# TODO: Test if rolling hash function with power of two modulu's would work.
struct KmerContent[KMERSIZE: Int](Copyable, Movable):
    var kmers: List[List[Int64]]
    var max_length: Int

    fn __init__(out self):
        self.kmers = List[List[Int64]](capacity=pow(4, Self.KMERSIZE))
        for _ in range(pow(4, Self.KMERSIZE)):
            var one = List[Int64](capacity=1)
            one.append(0)
            self.kmers.append(one^)
        self.max_length = 0

    @always_inline
    fn tally_read(mut self, record: FastqRecord, read_num: Int64):
        comptime N_b = ord("N")
        comptime n_b = ord("n")

        if read_num % 50 != 0:
            return

        if len(record) > self.max_length:
            self.max_length = len(record)
            var new_kmers = List[List[Int64]](capacity=pow(4, Self.KMERSIZE))
            for i in range(pow(4, Self.KMERSIZE)):
                new_kmers.append(grow_tensor(self.kmers[i], self.max_length))

            self.kmers = new_kmers^

        var s = record.get_seq().as_bytes()
        # INFO: As per FastQC: limit the Kmers to the first 500 BP for long reads
        for i in range(min(len(record), 500) - Self.KMERSIZE):
            var kmer = s[i : i + Self.KMERSIZE]
            var contains_n = False
            # TODO: Check how this optimized in Cpp
            for l in kmer:
                if l == N_b or l == n_b:
                    contains_n = True
            if contains_n:
                continue
            self.kmers[self.kmer_to_index(kmer)][i] += 1

    # From: https://github.com/smithlabcode/falco/blob/f4f0e6ca35e262cbeffc81fdfc620b3413ecfe2c/src/smithlab_utils.hpp#L357
    @always_inline
    fn kmer_to_index(self, kmer: Span[T=Byte]) -> Int:
        var index: UInt = 0
        var multiplier = 1
        for i in range(len(kmer), -1, -1):
            index += Int(base2int(kmer[i]))
            multiplier *= 4
        return Int(index)

    # TODO: Figure out how the enrichment calculation is carried out.
    # Check: https://github.com/smithlabcode/falco/blob/f4f0e6ca35e262cbeffc81fdfc620b3413ecfe2c/src/Module.cpp#L2068
    fn plot(self) raises -> PythonObject:
        var agg = Matrix2D(pow(4, Self.KMERSIZE), self.max_length)
        for i in range(len(self.kmers)):
            for j in range(len(self.kmers[i])):
                agg[Index(i, j)] = self.kmers[i][j]

        var mat = matrix_to_numpy(agg)
        return mat

    # TODO: Sort the Kmers to report
    # Ported from Falco C++ implementation: https://github.com/smithlabcode/falco/blob/f4f0e6ca35e262cbeffc81fdfc620b3413ecfe2c/src/Module.cpp#L2057
    fn get_kmer_stats(
        self, kmer_count: PythonObject, num_kmers: Int
    ) raises -> PythonObject:
        var np = Python.import_module("numpy")

        var min_obs_exp_to_report = 1e-2

        var num_kmer_bases = min(self.max_length, 500)
        var obs_exp_max = np.zeros(num_kmers, dtype=np.float64)
        var where_obs_exp_is_max = np.zeros(num_kmers, dtype=np.int32)
        var total_kmer_counts = np.zeros(num_kmers, dtype=np.int64)

        total_kmer_counts = kmer_count[:num_kmer_bases, :].sum(axis=0)
        var num_seen_kmers = np.count_nonzero(total_kmer_counts)

        var dividend: PythonObject
        if Int(py=num_seen_kmers) > 0:
            dividend = np.float64(num_seen_kmers)
        else:
            dividend = np.finfo(np.float64).eps

        for pos in range(num_kmer_bases):
            var observed_counts = kmer_count[pos : pos + 1, :]
            var expected_counts = kmer_count[pos : pos + 1, :] / dividend
            var obs_exp_ratios = np.divide(
                observed_counts,
                expected_counts,
                out=np.zeros_like(observed_counts, dtype=np.float64),
                where=expected_counts > 0,
            )

            # Update maximum obs/exp ratios and positions
            var mask = obs_exp_ratios > obs_exp_max
            obs_exp_max = np.where(mask, obs_exp_ratios, obs_exp_max)
            where_obs_exp_is_max = np.where(
                mask, pos + 1 - 7, where_obs_exp_is_max
            )

        var kmers_to_report = Python.list()

        # Filter and sort k-mers with significant obs/exp ratios
        var significant_mask = obs_exp_max > min_obs_exp_to_report
        var significant_kmers = np.where(significant_mask)[0]
        for kmer in significant_kmers:
            kmers_to_report.append(Python.tuple(kmer, obs_exp_max[kmer]))

        return kmers_to_report
