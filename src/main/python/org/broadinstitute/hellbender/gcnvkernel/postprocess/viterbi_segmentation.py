import itertools
import logging
import os
import shutil
from typing import List, Tuple, Dict, TypeVar, Generator

import numpy as np

import vcf

from .segment_quality_utils import HMMSegmentationQualityCalculator
from .. import types
from ..io import io_consts, io_commons, io_denoising_calling, io_intervals_and_counts, io_vcf_parsing
from ..models.model_denoising_calling import DenoisingModelConfig, CopyNumberCallingConfig, \
    HHMMClassAndCopyNumberBasicCaller
from ..models.theano_hmm import TheanoForwardBackward, TheanoViterbi
from ..structs.interval import Interval
from ..structs.metadata import IntervalListMetadata
from ..structs.metadata import SampleMetadataCollection
from ..structs.segment import IntegerCopyNumberSegment

_logger = logging.getLogger(__name__)


class ViterbiSegmentationEngine:
    """This class runs the forward-backward and Viterbi algorithm on gCNV model/calls shards for a single sample,
     obtains constant copy-number segments, calculates various quality metrics, and saves the result to disk.

    Note:
        It is assumed that the model and calls shards are provided in order according to the SAM sequence dictionary.
        It is not checked or enforced here.
    """
    def __init__(self,
                 model_shards_paths: List[str],
                 calls_shards_paths: List[str],
                 sample_metadata_collection: SampleMetadataCollection,
                 sample_index: int,
                 output_path: str,
                 intervals_vcf: str = None,
                 clustered_vcf: str = None):
        """Initializer.
        
        Args:
            model_shards_paths: list of paths to model shards
            calls_shards_paths: list of paths to calls shards
            sample_metadata_collection: sample metadata collection (must contain sample being analyzed)
            sample_index: index of the sample in the callset
            output_path: output path for writing segmentation results
            intervals_vcf: file with single-sample copy number calls for all intervals
            clustered_vcf: file with clustered breakpoints and calls for each sample
        """
        try:
            self._validate_args(model_shards_paths, calls_shards_paths, sample_metadata_collection, sample_index,
                                clustered_vcf)
        except AssertionError as ex:
            raise AssertionError("Inconsistency detected in the provided model and calls shards.") from ex

        self.sample_index = sample_index
        self.output_path = output_path
        self.calls_shards_paths = calls_shards_paths
        self.sample_metadata_collection = sample_metadata_collection
        self.denoising_config = self._get_denoising_config(model_shards_paths[0])
        self.calling_config = self._get_calling_config(model_shards_paths[0])
        self.intervals_vcf = intervals_vcf
        self.clustered_vcf = clustered_vcf

        # assemble scattered global entities (interval list, log_q_tau_tk)
        _logger.info("Assembling interval list and copy-number class posterior from model shards...")
        self.interval_list: List[Interval] = []
        log_q_tau_tk_shards: Tuple[np.ndarray] = ()
        for model_path in model_shards_paths:
            self.interval_list += self._get_interval_list_from_model_shard(model_path)
            log_q_tau_tk_shards += (self._get_log_q_tau_tk_from_model_shard(model_path),)
        self.log_q_tau_tk: np.ndarray = np.concatenate(log_q_tau_tk_shards, axis=0)

        # extract SAM header lines from one of the interval lists
        self.interval_list_sam_header_lines = io_intervals_and_counts.extract_sam_header_from_file(
            os.path.join(model_shards_paths[0], io_consts.default_interval_list_filename))

        # sample names
        self.sample_name = self._get_sample_name_from_calls_shard(calls_shards_paths[0], sample_index)

        # interval list metadata
        interval_list_metadata: IntervalListMetadata = IntervalListMetadata(self.interval_list)
        self.ordered_contig_list = interval_list_metadata.ordered_contig_list
        self.contig_interval_indices = interval_list_metadata.contig_interval_indices
        self.contig_interval_lists: Dict[str, List[Interval]] = {
            contig: [self.interval_list[ti] for ti in self.contig_interval_indices[contig]]
            for contig in self.ordered_contig_list}

        # cnv stay probability for each contig
        self.cnv_stay_prob_t_j: Dict[str, np.ndarray] = dict()
        for contig in self.ordered_contig_list:
            contig_interval_list = self.contig_interval_lists[contig]
            dist_t = np.asarray([contig_interval_list[ti + 1].distance(contig_interval_list[ti])
                                 for ti in range(len(contig_interval_list) - 1)], dtype=types.floatX)
            self.cnv_stay_prob_t_j[contig] = np.exp(-dist_t / self.calling_config.cnv_coherence_length)

        # forward-backward algorithm
        _logger.info("Compiling theano forward-backward function...")
        self.theano_forward_backward = TheanoForwardBackward(
            log_posterior_probs_output_tc=None,
            resolve_nans=False,
            do_thermalization=False,
            do_admixing=False,
            include_update_size_output=False,
            include_alpha_beta_output=True)

        # viterbi algorithm
        _logger.info("Compiling theano Viterbi function...")
        self.theano_viterbi = TheanoViterbi()

        # copy-number HMM specs generator
        _logger.info("Compiling theano variational HHMM...")
        self.get_copy_number_hmm_specs = HHMMClassAndCopyNumberBasicCaller\
            .get_compiled_copy_number_hmm_specs_theano_func()

        # initialize log likelihood
        self.log_likelihood = 0.

    def _viterbi_segments_generator(self) -> Generator[IntegerCopyNumberSegment, None, None]:
        """Performs Viterbi segmentation and segment quality calculation for a single sample in
        the call-set and returns a generator for segments.

        Returns:
            a generator for segments
        """

        # load copy number log emission for the sample
        copy_number_log_emission_tc_shards = ()
        for calls_path in self.calls_shards_paths:
            copy_number_log_emission_tc_shards += (self._get_log_copy_number_emission_tc_from_calls_shard(
                calls_path, self.sample_index),)
        copy_number_log_emission_tc = np.concatenate(copy_number_log_emission_tc_shards, axis=0)

        # iterate over contigs and perform segmentation
        sample_name = self.sample_name
        for contig_index, contig in enumerate(self.ordered_contig_list):
            _logger.info("Segmenting contig ({0}/{1}) (contig name: {2})...".format(
                contig_index + 1, len(self.ordered_contig_list), contig))

            # copy-number prior probabilities for each class
            contig_baseline_copy_number = self.sample_metadata_collection\
                .get_sample_ploidy_metadata(sample_name)\
                .get_contig_ploidy(contig)
            pi_jkc = HHMMClassAndCopyNumberBasicCaller.get_copy_number_prior_for_sample_jkc(
                self.calling_config.num_copy_number_states,
                self.calling_config.p_alt,
                np.asarray([contig_baseline_copy_number], dtype=types.med_uint))

            # contig interval list and indices
            contig_interval_list = self.contig_interval_lists[contig]
            contig_interval_indices = self.contig_interval_indices[contig]

            # mapping from intervals to contig index (since we have a single contig, all intervals map to index=0)
            t_to_j_map = np.zeros((len(contig_interval_list),), dtype=types.med_uint)

            # copy-number class log probability
            log_q_tau_tk = self.log_q_tau_tk[contig_interval_indices, :]

            # copy-number log emission probability for contig intervals
            copy_number_log_emission_contig_tc = copy_number_log_emission_tc[contig_interval_indices, :]

            # get HMM specs
            hmm_specs = self.get_copy_number_hmm_specs(
                pi_jkc, self.cnv_stay_prob_t_j[contig], log_q_tau_tk, t_to_j_map)
            log_prior_c = hmm_specs[0]
            log_trans_contig_tcc = hmm_specs[1]

            # run forward-back algorithm
            fb_result = self.theano_forward_backward.perform_forward_backward(
                log_prior_c, log_trans_contig_tcc, copy_number_log_emission_contig_tc)
            log_posterior_prob_tc = fb_result.log_posterior_probs_tc
            log_data_likelihood = fb_result.log_data_likelihood
            alpha_tc = fb_result.alpha_tc
            beta_tc = fb_result.beta_tc
            self.log_likelihood += log_data_likelihood

            # initialize the segment quality calculator
            segment_quality_calculator: HMMSegmentationQualityCalculator = HMMSegmentationQualityCalculator(
                copy_number_log_emission_contig_tc, log_trans_contig_tcc,
                alpha_tc, beta_tc, log_posterior_prob_tc, log_data_likelihood)

            if self.clustered_vcf is None or self.intervals_vcf is None:
                # validate args -- should be both none or neither none
                if bool(self.clustered_vcf is None) != bool(self.intervals_vcf is None):
                    raise Exception("If clustered_vcf is provided, then intervals_vcf must be provided.")
                # run viterbi algorithm
                viterbi_path_t_contig = self.theano_viterbi.get_viterbi_path(
                    log_prior_c, log_trans_contig_tcc, copy_number_log_emission_contig_tc)

                # coalesce into piecewise constant copy-number segments
                segments = self._coalesce_seq_into_segments(viterbi_path_t_contig)
            else:
                # use events from clustered_vcf
                segments = io_vcf_parsing.read_sample_segments_and_calls(self.intervals_vcf, self.clustered_vcf, self.sample_name, contig)

            # calculate qualities
            for call_copy_number, start_index, end_index in segments:
                num_points = end_index - start_index + 1
                try:
                    segment = IntegerCopyNumberSegment(contig,
                                                       contig_interval_list[start_index].start,
                                                       contig_interval_list[end_index].end,
                                                       num_points,
                                                       call_copy_number,
                                                       contig_baseline_copy_number)
                except IndexError:
                    print("end index out of bounds: {0} requested, max is {1}".format(end_index, len(contig_interval_list)))
                if num_points > 1:
                    segment.quality_some_called = segment_quality_calculator.get_segment_quality_some_called(
                        start_index, end_index, call_copy_number)
                    segment.quality_all_called = segment_quality_calculator.get_segment_quality_all_called(
                        start_index, end_index, call_copy_number)
                    segment.quality_start = segment_quality_calculator.get_segment_quality_start(
                        start_index, call_copy_number)
                    segment.quality_end = segment_quality_calculator.get_segment_quality_end(
                        end_index, call_copy_number)

                else:  # for single-interval segments, all qualities must be the same
                    segment.quality_some_called = segment_quality_calculator.get_segment_quality_some_called(
                        start_index, end_index, call_copy_number)
                    segment.quality_all_called = segment.quality_some_called
                    segment.quality_start = segment.quality_some_called
                    segment.quality_end = segment.quality_some_called

                yield segment

    def write_results(self):
        """Performs Viterbi segmentation and segment quality calculation for a single sample in
        the call-set and saves the results to disk.

        """
        sample_name = self.sample_name
        _logger.info("Processing sample index: {0}, sample name: {1}...".format(self.sample_index, sample_name))
        sample_output_path = os.path.join(self.output_path, io_consts.sample_folder_prefix + repr(self.sample_index))
        io_commons.assert_output_path_writable(sample_output_path, try_creating_output_path=True)

        # write configs, gcnvkernel version and sample name to output path
        shutil.copy(os.path.join(self.calls_shards_paths[0], io_consts.default_denoising_config_json_filename),
                    sample_output_path)
        shutil.copy(os.path.join(self.calls_shards_paths[0], io_consts.default_calling_config_json_filename),
                    sample_output_path)
        io_commons.write_gcnvkernel_version(sample_output_path)
        io_commons.write_sample_name_to_txt_file(sample_output_path, sample_name)

        seg_file = os.path.join(sample_output_path, io_consts.default_copy_number_segments_tsv_filename)
        with open(seg_file, 'w') as of:
            # copy SAM header lines from model/calls interval list
            for sam_header_line in self.interval_list_sam_header_lines:
                of.write(sam_header_line + '\n')

            # add sample name header
            of.write('@' + io_consts.sample_name_sam_header_prefix + sample_name + '\n')

            # add table column headers
            of.write(IntegerCopyNumberSegment.get_header_column_string() + '\n')

            # add segments
            for segment in self._viterbi_segments_generator():
                of.write(repr(segment) + '\n')

        log_likelihood_file = os.path.join(sample_output_path, io_consts.default_log_likelihood_txt_filename)
        with open(log_likelihood_file, 'w') as of:
            of.write(repr(self.log_likelihood) + '\n')

    @staticmethod
    def _validate_args(model_shards_paths: List[str],
                       calls_shards_paths: List[str],
                       sample_metadata_collection: SampleMetadataCollection,
                       sample_index: int,
                       clustered_vcf: str):
        assert len(model_shards_paths) > 0, "At least one model shard must be provided."
        assert len(calls_shards_paths) == len(model_shards_paths),\
            "The number of model shards ({0}) and calls shards ({1}) must match.".format(
                len(model_shards_paths), len(calls_shards_paths))
        assert sample_index >= 0, "Sample index must be an integer non-negative number"

        scattered_sample_names: List[str] = []
        for model_path, calls_path in zip(model_shards_paths, calls_shards_paths):
            # assert interval lists are identical
            model_interval_list_file = os.path.join(model_path, io_consts.default_interval_list_filename)
            calls_interval_list_file = os.path.join(calls_path, io_consts.default_interval_list_filename)
            io_commons.assert_files_are_identical(model_interval_list_file, calls_interval_list_file)

            # assert gcnvkernel versions are identical
            model_gcnvkernel_version_file = os.path.join(model_path, io_consts.default_gcnvkernel_version_json_filename)
            calls_gcnvkernel_version_file = os.path.join(calls_path, io_consts.default_gcnvkernel_version_json_filename)
            try:
                io_commons.assert_files_are_identical(model_gcnvkernel_version_file, calls_gcnvkernel_version_file)
            except AssertionError:
                _logger.warning("Different gcnvkernel versions between model and calls -- proceeding at your own risk!")

            # assert denoising configs are identical
            model_denoising_config_file = os.path.join(model_path, io_consts.default_denoising_config_json_filename)
            calls_denoising_config_file = os.path.join(calls_path, io_consts.default_denoising_config_json_filename)
            try:
                io_commons.assert_files_are_identical(model_denoising_config_file, calls_denoising_config_file)
            except AssertionError:
                _logger.warning("Different denoising configuration between model and calls -- "
                                "proceeding at your own risk!")

            # assert callings configs are identical
            model_calling_config_file = os.path.join(model_path, io_consts.default_calling_config_json_filename)
            calls_calling_config_file = os.path.join(calls_path, io_consts.default_calling_config_json_filename)
            try:
                io_commons.assert_files_are_identical(model_calling_config_file, calls_calling_config_file)
            except AssertionError:
                _logger.warning("Different calling configuration between model and calls -- "
                                "proceeding at your own risk!")

            # extract and store sample names for the current shard
            scattered_sample_names.append(
                ViterbiSegmentationEngine._get_sample_name_from_calls_shard(calls_path, sample_index))

        # all scattered calls have the same set of samples and in the same order
        assert len(set(scattered_sample_names)) == 1,\
            "The calls shards contain different sample names and/or different number of samples."

        if clustered_vcf is not None:
            clustered_reader = vcf.Reader(filename=clustered_vcf)
            assert set(clustered_reader.samples).issuperset(set(scattered_sample_names)), \
                "The clustered VCF does not contain all samples in the calls shard."

        # all samples have ploidy calls in the metadata collection
        sample_names = list(scattered_sample_names[0])
        sample_metadata_collection.all_samples_have_ploidy_metadata(sample_names)

    @staticmethod
    def _get_sample_name_from_calls_shard(calls_path: str, sample_index: int) -> str:
        sample_posteriors_path = io_denoising_calling.get_sample_posterior_path(calls_path, sample_index)
        if not os.path.isdir(sample_posteriors_path):
            raise Exception("Could not find any sample posterior calls in {0} for sample with index {1}.".
                            format(calls_path, sample_index))
        sample_name = io_commons.get_sample_name_from_txt_file(sample_posteriors_path)
        return sample_name

    @staticmethod
    def _get_denoising_config(input_path: str) -> DenoisingModelConfig:
        return DenoisingModelConfig.from_json_file(os.path.join(
            input_path, io_consts.default_denoising_config_json_filename))

    @staticmethod
    def _get_calling_config(input_path: str) -> CopyNumberCallingConfig:
        return CopyNumberCallingConfig.from_json_file(os.path.join(
            input_path, io_consts.default_calling_config_json_filename))

    @staticmethod
    def _get_interval_list_from_model_shard(model_path: str) -> List[Interval]:
        interval_list_file = os.path.join(model_path, io_consts.default_interval_list_filename)
        return io_intervals_and_counts.load_interval_list_tsv_file(interval_list_file)

    @staticmethod
    def _get_log_q_tau_tk_from_model_shard(model_path: str) -> np.ndarray:
        return io_commons.read_ndarray_from_tsv(os.path.join(
            model_path, io_consts.default_class_log_posterior_tsv_filename))

    @staticmethod
    def _get_log_copy_number_emission_tc_from_calls_shard(calls_path: str, sample_index: int):
        return io_denoising_calling.SampleDenoisingAndCallingPosteriorsReader.\
            read_ndarray_tc_with_copy_number_header(
                io_denoising_calling.get_sample_posterior_path(calls_path, sample_index),
                io_consts.default_copy_number_log_emission_tsv_filename)

    @staticmethod
    def _coalesce_seq_into_segments(seq: List[TypeVar('_T')]) -> List[Tuple[TypeVar('_T'), int, int]]:
        """Coalesces a sequence of objects into piecewise constant segments, along with start and end indices
        for each constant segment.

        Example:
            seq = ['a', 'a', 'a', 'a', 'b', 'c', 'c', 'a', 'a', 'a']
            result = [('a', 0, 3), ('b', 4, 4), ('c', 5, 6), ('a', 7, 9)]

        Args:
            seq: a sequence of objects that implement __equals__

        Returns:
            a generator for (object, start_index, end_index)
        """
        for seg in itertools.groupby(enumerate(seq), key=lambda elem: elem[1]):
            seg_const = seg[0]
            grouper = seg[1]
            start_index = grouper.__next__()[0]
            end_index = start_index
            try:
                while True:
                    end_index = grouper.__next__()[0]
            except StopIteration:
                pass
            yield (seg_const, start_index, end_index)
