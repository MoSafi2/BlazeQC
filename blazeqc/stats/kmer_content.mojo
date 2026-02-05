"""K-mer content (split from stats_.mojo)."""

from utils import Index
from memory import Span
from python import Python, PythonObject
from blazeseq import FastqRecord
from blazeqc.stats.analyser import py_lib
from blazeqc.helpers import base2int, grow_tensor, Matrix2D, matrix_to_numpy


# TODO: Test if rolling hash function with power of two modulu's would work.
@value
struct KmerContent[KMERSIZE: Int]:
    var kmers: List[List[Int64]]
    var max_length: Int

    fn __init__(out self):
        self.kmers = List[List[Int64]](capacity=pow(4, KMERSIZE))
        for _ in range(pow(4, KMERSIZE)):
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
            var new_kmers = List[List[Int64]](capacity=pow(4, KMERSIZE))
            for i in range(pow(4, KMERSIZE)):
                new_kmers.append(grow_tensor(self.kmers[i], self.max_length))

            self.kmers = new_kmers^

        var s = record.get_seq().as_bytes()
        # INFO: As per FastQC: limit the Kmers to the first 500 BP for long reads
        for i in range(min(record.len_record(), 500) - KMERSIZE):
            var kmer = s[i : i + KMERSIZE]
            var contains_n = False
            # TODO: Check how this optimized in Cpp
            for l in kmer:
                if l[] == N_b or l[] == n_b:
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
        return index

    # TODO: Figure out how the enrichment calculation is carried out.
    # Check: https://github.com/smithlabcode/falco/blob/f4f0e6ca35e262cbeffc81fdfc620b3413ecfe2c/src/Module.cpp#L2068
    fn plot(self) raises -> PythonObject:
        Python.add_to_path(py_lib.as_string_slice())
        var agg = Matrix2D(pow(4, KMERSIZE), self.max_length)
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
        Python.add_to_path(py_lib.as_string_slice())
        var np = Python.import_module("numpy")

        min_obs_exp_to_report = 1e-2

        num_kmer_bases = min(self.max_length, 500)
        obs_exp_max = np.zeros(num_kmers, dtype=np.float64)
        where_obs_exp_is_max = np.zeros(num_kmers, dtype=np.int32)
        total_kmer_counts = np.zeros(num_kmers, dtype=np.int64)

        total_kmer_counts = kmer_count[:num_kmer_bases, :].sum(axis=0)
        num_seen_kmers = np.count_nonzero(total_kmer_counts)

        dividend = (
            Float64(num_seen_kmers) if num_seen_kmers
            > 0 else np.finfo(np.float64).eps
        )

        for pos in range(num_kmer_bases):
            observed_counts = kmer_count[pos : pos + 1, :]
            expected_counts = kmer_count[pos : pos + 1, :] / dividend
            obs_exp_ratios = np.divide(
                observed_counts,
                expected_counts,
                out=np.zeros_like(observed_counts, dtype=np.float64),
                where=expected_counts > 0,
            )

            # Update maximum obs/exp ratios and positions
            mask = obs_exp_ratios > obs_exp_max
            obs_exp_max = np.where(mask, obs_exp_ratios, obs_exp_max)
            where_obs_exp_is_max = np.where(
                mask, pos + 1 - 7, where_obs_exp_is_max
            )

        kmers_to_report = PythonObject([])

        # Filter and sort k-mers with significant obs/exp ratios
        significant_mask = obs_exp_max > min_obs_exp_to_report
        significant_kmers = np.where(significant_mask)[0]
        for kmer in significant_kmers:
            kmers_to_report.append((kmer, obs_exp_max[kmer]))

        return kmers_to_report
