import argparse
import collections
import inspect
import json
import logging
from abc import abstractmethod
from typing import List, Tuple, Set, Dict, Optional

import numpy as np
import pymc3 as pm
import scipy.sparse as sp
import theano as th
import theano.sparse as tst
import theano.tensor as tt
from pymc3 import Normal, Deterministic, DensityDist, Lognormal, Exponential

from . import commons
from .dists import HalfFlat
from .fancy_model import GeneralizedContinuousModel
from .theano_hmm import TheanoForwardBackward
from .. import config, types
from ..structs.interval import Interval, GCContentAnnotation
from ..structs.metadata import SampleMetadataCollection
from ..tasks.inference_task_base import HybridInferenceParameters

_logger = logging.getLogger(__name__)

_eps = commons.eps


class DenoisingModelConfig:
    """Configuration for the coverage denoising model, including hyper-parameters, model feature selection,
    and choice of approximation schemes."""

    # approximation schemes for calculating expectations with respect to copy number posteriors
    _q_c_expectation_modes = ['map', 'exact', 'hybrid']

    def __init__(self,
                 max_bias_factors: int = 5,
                 mapping_error_rate: float = 0.01,
                 psi_t_scale: float = 0.001,
                 psi_s_scale: float = 0.0001,
                 depth_correction_tau: float = 10000.0,
                 log_mean_bias_std: float = 0.1,
                 init_ard_rel_unexplained_variance: float = 0.1,
                 num_gc_bins: int = 20,
                 gc_curve_sd: float = 1.0,
                 q_c_expectation_mode: str = 'hybrid',
                 active_class_padding_hybrid_mode: int = 50000,
                 enable_bias_factors: bool = True,
                 enable_explicit_gc_bias_modeling: bool = False,
                 disable_bias_factors_in_active_class: bool = False):
        """See `expose_args` for the description of arguments"""
        self.max_bias_factors = max_bias_factors
        self.mapping_error_rate = mapping_error_rate
        self.psi_t_scale = psi_t_scale
        self.psi_s_scale = psi_s_scale
        self.depth_correction_tau = depth_correction_tau
        self.log_mean_bias_std = log_mean_bias_std
        self.init_ard_rel_unexplained_variance = init_ard_rel_unexplained_variance
        self.num_gc_bins = num_gc_bins
        self.gc_curve_sd = gc_curve_sd
        self.q_c_expectation_mode = q_c_expectation_mode
        self.active_class_padding_hybrid_mode = active_class_padding_hybrid_mode
        self.enable_bias_factors = enable_bias_factors
        self.enable_explicit_gc_bias_modeling = enable_explicit_gc_bias_modeling
        self.disable_bias_factors_in_active_class = disable_bias_factors_in_active_class

    @staticmethod
    def expose_args(args: argparse.ArgumentParser,
                    hide: Set[str] = None):
        """Exposes arguments of `__init__` to a given instance of `ArgumentParser`.

        Args:
            args: an instance of `ArgumentParser`
            hide: a set of arguments not to expose

        Returns:
            None
        """
        group = args.add_argument_group(title="Coverage denoising model parameters")
        if hide is None:
            hide = set()

        initializer_params = inspect.signature(DenoisingModelConfig.__init__).parameters
        valid_args = {"--" + arg for arg in initializer_params.keys()}
        for hidden_arg in hide:
            assert hidden_arg in valid_args, \
                "Initializer argument to be hidden {0} is not a valid initializer arguments; possible " \
                "choices are: {1}".format(hidden_arg, valid_args)

        def process_and_maybe_add(arg, **kwargs):
            full_arg = "--" + arg
            if full_arg in hide:
                return
            kwargs['default'] = initializer_params[arg].default
            group.add_argument(full_arg, **kwargs)

        def str_to_bool(value: str):
            if value.lower() in ('yes', 'true', 't', 'y', '1'):
                return True
            elif value.lower() in ('no', 'false', 'f', 'n', '0'):
                return False
            else:
                raise argparse.ArgumentTypeError('Boolean value expected.')

        process_and_maybe_add("max_bias_factors",
                              type=int,
                              help="Maximum number of bias factors")

        process_and_maybe_add("mapping_error_rate",
                              type=float,
                              help="Typical mapping error rate")

        process_and_maybe_add("psi_t_scale",
                              type=float,
                              help="Typical scale of interval-specific unexplained variance")

        process_and_maybe_add("psi_s_scale",
                              type=float,
                              help="Typical scale of sample-specific unexplained variance")

        process_and_maybe_add("depth_correction_tau",
                              type=float,
                              help="Precision of pinning read-depth in the coverage denoising model "
                                   "to its globally determined value")

        process_and_maybe_add("log_mean_bias_std",
                              type=float,
                              help="Standard deviation of mean bias in log space")

        process_and_maybe_add("init_ard_rel_unexplained_variance",
                              type=float,
                              help="Initial value of automatic relevance determination (ARD) precisions relative "
                                   "to the typical interval-specific unexplained variance scale")

        process_and_maybe_add("num_gc_bins",
                              type=int,
                              help="Number of knobs on the GC curve")

        process_and_maybe_add("gc_curve_sd",
                              type=float,
                              help="Prior standard deviation of the GC curve from a flat curve")

        process_and_maybe_add("q_c_expectation_mode",
                              type=str,
                              choices=DenoisingModelConfig._q_c_expectation_modes,
                              help="The strategy for calculating copy number posterior expectations in the denoising "
                                   "model. Choices: \"exact\": summation over all states, \"map\": drop all terms "
                                   "except for the maximum a posteriori (MAP) copy number estimate, \"hybrid\": "
                                   "use MAP strategy in silent regions and exact strategy in active regions.")

        process_and_maybe_add("active_class_padding_hybrid_mode",
                              type=int,
                              help="If q_c_expectation_mode is set to \"hybrid\", the active intervals "
                                   "will be further padded by this value (in the units of bp) in order to achieve "
                                   "higher sensitivity in detecting common CNVs and to avoid boundary artifacts")

        process_and_maybe_add("enable_bias_factors",
                              type=str_to_bool,
                              help="Enable discovery of novel bias factors")

        process_and_maybe_add("enable_explicit_gc_bias_modeling",
                              type=str_to_bool,
                              help="Enable explicit modeling of GC bias (if enabled, the provided modeling interval "
                                   "list of contain a column for {0} values)".format(GCContentAnnotation.get_key()))

        process_and_maybe_add("disable_bias_factors_in_active_class",
                              type=str_to_bool,
                              help="Disable novel bias factor discovery CNV-active regions")

    @staticmethod
    def from_args_dict(args_dict: Dict) -> 'DenoisingModelConfig':
        """Initialize an instance of `DenoisingModelConfig` from a dictionary of arguments.

        Args:
            args_dict: a dictionary of arguments; the keys must match argument names in
                `DenoisingModelConfig.__init__`

        Returns:
            an instance of `DenoisingModelConfig`
        """
        relevant_keys = set(inspect.getfullargspec(DenoisingModelConfig.__init__).args)
        relevant_kwargs = {k: v for k, v in args_dict.items() if k in relevant_keys}
        return DenoisingModelConfig(**relevant_kwargs)

    @staticmethod
    def from_json_file(json_file: str) -> 'DenoisingModelConfig':
        with open(json_file, 'r') as fp:
            imported_denoising_config_dict = json.load(fp)
        return DenoisingModelConfig.from_args_dict(imported_denoising_config_dict)


class CopyNumberCallingConfig:
    """Configuration of the copy number caller."""
    def __init__(self,
                 p_alt: float = 1e-6,
                 p_active: float = 1e-3,
                 cnv_coherence_length: float = 10000.0,
                 class_coherence_length: float = 10000.0,
                 max_copy_number: int = 5,
                 num_calling_processes: int = 1):
        """See `expose_args` for the description of arguments"""
        assert 0.0 <= p_alt <= 1.0
        assert 0.0 <= p_active <= 1.0
        assert cnv_coherence_length > 0.0
        assert class_coherence_length > 0.0
        assert max_copy_number > 0
        assert max_copy_number * p_alt < 1.0
        assert num_calling_processes > 0

        self.p_alt = p_alt
        self.p_active = p_active
        self.cnv_coherence_length = cnv_coherence_length
        self.class_coherence_length = class_coherence_length
        self.max_copy_number = max_copy_number
        self.num_calling_processes = num_calling_processes

        self.num_copy_number_states = max_copy_number + 1
        self.num_copy_number_classes = 2

    @staticmethod
    def expose_args(args: argparse.ArgumentParser, hide: Set[str] = None):
        """Exposes arguments of `__init__` to a given instance of `ArgumentParser`.

        Args:
            args: an instance of `ArgumentParser`
            hide: a set of arguments not to expose

        Returns:
            None
        """
        group = args.add_argument_group(title="Copy number calling parameters")
        if hide is None:
            hide = set()

        initializer_params = inspect.signature(CopyNumberCallingConfig.__init__).parameters
        valid_args = {"--" + arg for arg in initializer_params.keys()}
        for hidden_arg in hide:
            assert hidden_arg in valid_args, \
                "Initializer argument to be hidden {0} is not a valid initializer arguments; possible " \
                "choices are: {1}".format(hidden_arg, valid_args)

        def process_and_maybe_add(arg, **kwargs):
            full_arg = "--" + arg
            if full_arg in hide:
                return
            kwargs['default'] = initializer_params[arg].default
            group.add_argument(full_arg, **kwargs)

        def str_to_bool(value: str):
            if value.lower() in ('yes', 'true', 't', 'y', '1'):
                return True
            elif value.lower() in ('no', 'false', 'f', 'n', '0'):
                return False
            else:
                raise argparse.ArgumentTypeError('Boolean value expected.')

        process_and_maybe_add("p_alt",
                              type=float,
                              help="Prior probability of alternate copy number with respect to contig baseline "
                                   "state in CNV-silent intervals")

        process_and_maybe_add("p_active",
                              type=float,
                              help="Prior probability of treating an interval as CNV-active")

        process_and_maybe_add("cnv_coherence_length",
                              type=float,
                              help="Coherence length of CNV events (in the units of bp)")

        process_and_maybe_add("class_coherence_length",
                              type=float,
                              help="Coherence length of CNV-silent and CNV-active domains (in the units of bp)")

        process_and_maybe_add("max_copy_number",
                              type=int,
                              help="Highest called copy number state")

        process_and_maybe_add("num_calling_processes",
                              type=int,
                              help="Number of concurrent forward-backward threads (not implemented yet)")

    @staticmethod
    def from_args_dict(args_dict: Dict):
        """Initialize an instance of `CopyNumberCallingConfig` from a dictionary of arguments.

        Args:
            args_dict: a dictionary of arguments; the keys must match argument names in
                `CopyNumberCallingConfig.__init__`

        Returns:
            an instance of `CopyNumberCallingConfig`
        """
        relevant_keys = set(inspect.getfullargspec(CopyNumberCallingConfig.__init__).args)
        relevant_kwargs = {k: v for k, v in args_dict.items() if k in relevant_keys}
        return CopyNumberCallingConfig(**relevant_kwargs)

    @staticmethod
    def from_json_file(json_file: str) -> 'CopyNumberCallingConfig':
        with open(json_file, 'r') as fp:
            imported_calling_config_dict = json.load(fp)
        return CopyNumberCallingConfig.from_args_dict(imported_calling_config_dict)


class PosteriorInitializer:
    """Base class for posterior initializers."""
    @staticmethod
    @abstractmethod
    def initialize_posterior(denoising_config: DenoisingModelConfig,
                             calling_config: CopyNumberCallingConfig,
                             shared_workspace: 'DenoisingCallingWorkspace') -> None:
        raise NotImplementedError


class TrivialPosteriorInitializer(PosteriorInitializer):
    """Initialize posteriors to reasonable values based on priors."""
    @staticmethod
    def initialize_posterior(denoising_config: DenoisingModelConfig,
                             calling_config: CopyNumberCallingConfig,
                             shared_workspace: 'DenoisingCallingWorkspace'):
        # interval class log posterior probs
        class_probs_k = np.asarray([1.0 - calling_config.p_active, calling_config.p_active], dtype=types.floatX)
        log_q_tau_tk = np.tile(np.log(class_probs_k), (shared_workspace.num_intervals, 1))
        shared_workspace.log_q_tau_tk = th.shared(log_q_tau_tk, name="log_q_tau_tk", borrow=config.borrow_numpy)

        # copy number log posterior probs
        log_q_c_stc = np.zeros((shared_workspace.num_samples, shared_workspace.num_intervals,
                                calling_config.num_copy_number_states), dtype=types.floatX)
        t_to_j_map = shared_workspace.t_to_j_map.get_value(borrow=True)
        for si in range(shared_workspace.num_samples):
            sample_baseline_copy_number_j = shared_workspace.baseline_copy_number_sj[si, :]
            sample_pi_jkc = HHMMClassAndCopyNumberBasicCaller.get_copy_number_prior_for_sample_jkc(
                calling_config.num_copy_number_states,
                calling_config.p_alt,
                sample_baseline_copy_number_j)
            sample_log_pi_jc = np.log(np.sum(sample_pi_jkc * class_probs_k[np.newaxis, :, np.newaxis], axis=1))
            for ti in range(shared_workspace.num_intervals):
                log_q_c_stc[si, ti, :] = sample_log_pi_jc[t_to_j_map[ti], :]
        shared_workspace.log_q_c_stc = th.shared(log_q_c_stc, name="log_q_c_stc", borrow=config.borrow_numpy)


class DenoisingCallingWorkspace:
    """This class contains objects (numpy arrays, theano tensors, etc) shared between the denoising model
    and the copy number caller."""
    def __init__(self,
                 denoising_config: DenoisingModelConfig,
                 calling_config: CopyNumberCallingConfig,
                 interval_list: List[Interval],
                 n_st: np.ndarray,
                 sample_names: List[str],
                 sample_metadata_collection: SampleMetadataCollection,
                 posterior_initializer: Optional[PosteriorInitializer] = TrivialPosteriorInitializer):
        self.denoising_config = denoising_config
        self.calling_config = calling_config
        self.interval_list = interval_list
        self.sample_names = sample_names

        assert n_st.ndim == 2, "Read counts matrix must be a 2-dim ndarray with shape (num_samples, num_intervals)"

        self.num_samples: int = n_st.shape[0]
        self.num_intervals: int = n_st.shape[1]

        assert self.num_intervals >= 2, "At least two intervals must be provided"
        assert len(interval_list) == self.num_intervals,\
            "The length of the interval list is incompatible with the shape of the read counts matrix"

        # a list of unique contigs appearing in the interval list; the ordering is arbitrary and
        # is only used internally
        # Note: j is the index subscript used for contig index hereafter
        self.contig_list = list({interval.contig for interval in interval_list})
        self.num_contigs = len(self.contig_list)
        contig_to_j_map = {contig: self.contig_list.index(contig) for contig in self.contig_list}
        t_to_j_map = np.asarray([contig_to_j_map[interval.contig] for interval in interval_list],
                                dtype=types.small_uint)
        self.t_to_j_map: types.TensorSharedVariable = th.shared(
            t_to_j_map, name="t_to_j_map", borrow=config.borrow_numpy)

        self.global_read_depth_s, average_ploidy_s, self.baseline_copy_number_sj = \
            DenoisingCallingWorkspace._get_baseline_copy_number_and_read_depth(
                sample_metadata_collection, sample_names, self.contig_list)

        max_baseline_copy_number = np.max(self.baseline_copy_number_sj)
        assert max_baseline_copy_number <= calling_config.max_copy_number, \
            "The highest contig ploidy ({0}) must be smaller or equal to the highest copy number state ({1})".format(
                max_baseline_copy_number, calling_config.max_copy_number)

        # shared theano tensors from the input data
        self.n_st: types.TensorSharedVariable = th.shared(
            n_st.astype(types.med_uint), name="n_st", borrow=config.borrow_numpy)
        self.average_ploidy_s: types.TensorSharedVariable = th.shared(
            average_ploidy_s.astype(types.floatX), name="average_ploidy_s", borrow=config.borrow_numpy)

        # copy-number event stay probability
        self.dist_t = np.asarray([self.interval_list[ti + 1].distance(self.interval_list[ti])
                                  for ti in range(self.num_intervals - 1)])
        cnv_stay_prob_t = np.exp(-self.dist_t / calling_config.cnv_coherence_length)
        self.cnv_stay_prob_t = th.shared(cnv_stay_prob_t, name='cnv_stay_prob_t', borrow=config.borrow_numpy)

        # copy number values for each copy number state
        copy_number_values_c = np.arange(0, calling_config.num_copy_number_states, dtype=types.small_uint)
        self.copy_number_values_c = th.shared(copy_number_values_c, name='copy_number_values_c',
                                              borrow=config.borrow_numpy)

        # copy number log posterior and derived quantities (to be initialized by `PosteriorInitializer`)
        self.log_q_c_stc: Optional[types.TensorSharedVariable] = None

        # latest MAP estimate of integer copy number (to be initialized and periodically updated by
        #   `DenoisingCallingWorkspace.update_auxiliary_vars)
        self.c_map_st: Optional[types.TensorSharedVariable] = None

        # latest bitmask of CNV-active intervals (to be initialized and periodically updated by
        #   `DenoisingCallingWorkspace.update_auxiliary_vars if q_c_expectation_mode == 'hybrid')
        self.active_class_bitmask_t: Optional[types.TensorSharedVariable] = None

        # copy number emission log posterior
        log_copy_number_emission_stc = np.zeros(
            (self.num_samples, self.num_intervals, calling_config.num_copy_number_states), dtype=types.floatX)
        self.log_copy_number_emission_stc: types.TensorSharedVariable = th.shared(
            log_copy_number_emission_stc, name="log_copy_number_emission_stc", borrow=config.borrow_numpy)

        # class log posterior (to be initialized by `PosteriorInitializer`)
        self.log_q_tau_tk: Optional[types.TensorSharedVariable] = None

        # class emission log posterior
        # (to be initialized by calling `initialize_copy_number_class_inference_vars`)
        self.log_class_emission_tk: Optional[types.TensorSharedVariable] = None

        # class assignment prior probabilities
        # (to be initialized by calling `initialize_copy_number_class_inference_vars`)
        self.class_probs_k: Optional[types.TensorSharedVariable] = None

        # class Markov chain log prior (initialized here and remains constant throughout)
        # (to be initialized by calling `initialize_copy_number_class_inference_vars`)
        self.log_prior_k: Optional[np.ndarray] = None

        # class Markov chain log transition (initialized here and remains constant throughout)
        # (to be initialized by calling `initialize_copy_number_class_inference_vars`)
        self.log_trans_tkk: Optional[np.ndarray] = None

        # GC bias factors
        # (to be initialized by calling `initialize_bias_inference_vars`)
        self.W_gc_tg: Optional[tst.SparseConstant] = None

        # auxiliary data structures for hybrid q_c_expectation_mode calculation
        # (to be initialized by calling `initialize_bias_inference_vars`)
        self.interval_neighbor_index_list: Optional[List[List[int]]] = None

        # denoised copy ratios
        denoised_copy_ratio_st = np.zeros((self.num_samples, self.num_intervals), dtype=types.floatX)
        self.denoised_copy_ratio_st: types.TensorSharedVariable = th.shared(
            denoised_copy_ratio_st, name="denoised_copy_ratio_st", borrow=config.borrow_numpy)

        # initialize posterior
        posterior_initializer.initialize_posterior(denoising_config, calling_config, self)
        self.initialize_bias_inference_vars()
        self.update_auxiliary_vars()

    def initialize_copy_number_class_inference_vars(self):
        """Initializes members required for copy number class inference (must be called in the cohort mode).
        The following members are initialized:
            - `DenoisingCallingWorkspace.log_class_emission_tk`
            - `DenoisingCallingWorkspace.class_probs_k`
            - `DenoisingCallingWorkspace.log_prior_k`
            - `DenoisingCallingWorkspace.log_trans_tkk`
        """
        # class emission log posterior
        log_class_emission_tk = np.zeros(
            (self.num_intervals, self.calling_config.num_copy_number_classes), dtype=types.floatX)
        self.log_class_emission_tk: types.TensorSharedVariable = th.shared(
            log_class_emission_tk, name="log_class_emission_tk", borrow=True)

        # class assignment prior probabilities
        # Note:
        #   The first class is the CNV-silent class (highly biased toward the baseline copy number)
        #   The second class is a CNV-active class (all copy number states are equally probable)
        class_probs_k = np.asarray([1.0 - self.calling_config.p_active, self.calling_config.p_active],
                                   dtype=types.floatX)
        self.class_probs_k: types.TensorSharedVariable = th.shared(
            class_probs_k, name='class_probs_k', borrow=config.borrow_numpy)

        # class Markov chain log prior (initialized here and remains constant throughout)
        self.log_prior_k: np.ndarray = np.log(class_probs_k)

        # class Markov chain log transition (initialized here and remains constant throughout)
        self.log_trans_tkk: np.ndarray = self._get_log_trans_tkk(
            self.dist_t,
            self.calling_config.class_coherence_length,
            self.calling_config.num_copy_number_classes,
            class_probs_k)

    def initialize_bias_inference_vars(self):
        """Initializes `DenoisingCallingWorkspace.W_gc_tg` and `DenoisingCallingWorkspace.interval_neighbor_index_list`
        if required by the model configuration."""
        if self.denoising_config.enable_explicit_gc_bias_modeling:
            self.W_gc_tg = self._create_sparse_gc_bin_tensor_tg(
                self.interval_list, self.denoising_config.num_gc_bins)

        if self.denoising_config.q_c_expectation_mode == 'hybrid':
            self.interval_neighbor_index_list = self._get_interval_neighbor_index_list(
                self.interval_list, self.denoising_config.active_class_padding_hybrid_mode)
        else:
            self.interval_neighbor_index_list = None

    def update_auxiliary_vars(self):
        """Updates `DenoisingCallingWorkspace.c_map_st' and `DenoisingCallingWorkspace.active_class_bitmask_t`."""
        # MAP copy number call
        if self.c_map_st is None:
            c_map_st = np.zeros((self.num_samples, self.num_intervals), dtype=types.small_uint)
            self.c_map_st = th.shared(c_map_st, name="c_map_st", borrow=config.borrow_numpy)
        self.c_map_st.set_value(
            np.argmax(self.log_q_c_stc.get_value(borrow=True), axis=2).astype(types.small_uint),
            borrow=config.borrow_numpy)

        if self.denoising_config.q_c_expectation_mode == 'hybrid':
            _logger.debug("Updating CNV-active class bitmask...")
            if self.active_class_bitmask_t is None:
                active_class_bitmask_t = np.zeros((self.num_intervals,), dtype=bool)
                self.active_class_bitmask_t = th.shared(
                    active_class_bitmask_t, name="active_class_bitmask_t", borrow=config.borrow_numpy)

            # bitmask for intervals of which the probability of being in the silent class is below 0.5
            active_class_bitmask_t: np.ndarray = \
                self.log_q_tau_tk.get_value(borrow=True)[:, 0] < -np.log(2)
            padded_active_class_bitmask_t = np.zeros_like(active_class_bitmask_t)
            for ti, neighbor_index_list in enumerate(self.interval_neighbor_index_list):
                padded_active_class_bitmask_t[ti] = np.any(active_class_bitmask_t[neighbor_index_list])
            self.active_class_bitmask_t.set_value(
                padded_active_class_bitmask_t, borrow=config.borrow_numpy)

    @staticmethod
    def _get_interval_neighbor_index_list(interval_list: List[Interval],
                                          maximum_neighbor_distance: int) -> List[List[int]]:
        """Pads a given interval list, finds the index of overlapping neighbors, and returns a list of indices of
        overlapping neighbors.

        Note:
            It is assumed that the `interval_list` is sorted (this is not asserted).

        Args:
            interval_list: list of intervals
            maximum_neighbor_distance: Maximum distance between intervals to be considered neighbors

        Returns:
            A list of indices of overlapping neighbors with the same length as `interval_list`. Each element
            in a variable-length list, depending on the number of neighbors.
        """
        assert maximum_neighbor_distance >= 0
        num_intervals = len(interval_list)
        padded_interval_list = [interval.get_padded(maximum_neighbor_distance) for interval in interval_list]
        interval_neighbor_index_list = []
        for ti, padded_interval in enumerate(padded_interval_list):
            overlapping_interval_indices = [ti]
            right_ti = ti
            while right_ti < num_intervals - 1:
                right_ti += 1
                if interval_list[right_ti].overlaps_with(padded_interval):
                    overlapping_interval_indices.append(right_ti)
                else:
                    break
            left_ti = ti
            while left_ti > 0:
                left_ti -= 1
                if interval_list[left_ti].overlaps_with(padded_interval):
                    overlapping_interval_indices.append(left_ti)
                else:
                    break
            interval_neighbor_index_list.append(overlapping_interval_indices)
        return interval_neighbor_index_list

    @staticmethod
    def _get_log_trans_tkk(dist_t: np.ndarray,
                           class_coherence_length: float,
                           num_copy_number_classes: int,
                           class_probs_k: np.ndarray) -> np.ndarray:
        """Calculates the log transition probability between copy number classes."""
        class_stay_prob_t = np.exp(-dist_t / class_coherence_length)
        class_not_stay_prob_t = np.ones_like(class_stay_prob_t) - class_stay_prob_t
        delta_kl = np.eye(num_copy_number_classes, dtype=types.floatX)
        trans_tkl = (class_not_stay_prob_t[:, None, None] * class_probs_k[None, None, :]
                     + class_stay_prob_t[:, None, None] * delta_kl[None, :, :])
        return np.log(trans_tkl)

    @staticmethod
    def _create_sparse_gc_bin_tensor_tg(interval_list: List[Interval], num_gc_bins: int) -> tst.SparseConstant:
        """Creates a sparse 2d theano tensor with shape (num_intervals, gc_bin). The sparse
        tensor represents a 1-hot mapping of each interval to its GC bin index. The range [0, 1]
        is uniformly divided into num_gc_bins.
        """
        assert all([GCContentAnnotation.get_key() in interval.annotations.keys() for interval in interval_list]), \
            "Explicit GC bias modeling is enabled, however, some or all intervals lack \"{0}\" annotation".format(
                GCContentAnnotation.get_key())

        def get_gc_bin_idx(gc_content):
            return min(int(gc_content * num_gc_bins), num_gc_bins - 1)

        num_intervals = len(interval_list)
        data = np.ones((num_intervals,))
        indices = [get_gc_bin_idx(interval.get_annotation(GCContentAnnotation.get_key()))
                   for interval in interval_list]
        indptr = np.arange(0, num_intervals + 1)
        scipy_gc_matrix = sp.csr_matrix((data, indices, indptr), shape=(num_intervals, num_gc_bins),
                                        dtype=types.small_uint)
        theano_gc_matrix: tst.SparseConstant = tst.as_sparse(scipy_gc_matrix)
        return theano_gc_matrix

    @staticmethod
    def _get_baseline_copy_number_and_read_depth(sample_metadata_collection: SampleMetadataCollection,
                                                 sample_names: List[str],
                                                 contig_list: List[str]) \
            -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Generates global read depth array, average ploidy array, and baseline copy numbers for all
        samples.

        Args:
            sample_metadata_collection: a instance of `SampleMetadataCollection` containing required metadata
                for all samples in `sample_names`
            sample_names: list of sample names
            contig_list: list of contigs appearing in the modeling interval list

        Returns:
            global read depth, average ploudy, baseline copy number
        """
        assert sample_metadata_collection.all_samples_have_read_depth_metadata(sample_names), \
            "Some samples do not have read depth metadata"
        assert sample_metadata_collection.all_samples_have_ploidy_metadata(sample_names), \
            "Some samples do not have ploidy metadata"
        num_samples = len(sample_names)
        num_contigs = len(contig_list)

        global_read_depth_s = np.zeros((num_samples,), dtype=types.floatX)
        average_ploidy_s = np.zeros((num_samples,), dtype=types.floatX)
        baseline_copy_number_sj = np.zeros((num_samples, num_contigs), dtype=types.small_uint)

        for si, sample_name in enumerate(sample_names):
            sample_read_depth_metadata = sample_metadata_collection.get_sample_read_depth_metadata(sample_name)
            sample_ploidy_metadata = sample_metadata_collection.get_sample_ploidy_metadata(sample_name)

            global_read_depth_s[si] = sample_read_depth_metadata.global_read_depth
            average_ploidy_s[si] = sample_read_depth_metadata.average_ploidy
            sample_baseline_copy_number_j = np.asarray([sample_ploidy_metadata.get_contig_ploidy(contig)
                                                        for contig in contig_list], dtype=types.small_uint)
            baseline_copy_number_sj[si, :] = sample_baseline_copy_number_j[:]

        return global_read_depth_s, average_ploidy_s, baseline_copy_number_sj


class InitialModelParametersSupplier:
    """Base class for suppliers of initial global model parameters"""
    def __init__(self,
                 denoising_model_config: DenoisingModelConfig,
                 calling_config: CopyNumberCallingConfig,
                 shared_workspace: DenoisingCallingWorkspace):
        self.denoising_model_config = denoising_model_config
        self.calling_config = calling_config
        self.shared_workspace = shared_workspace

    @abstractmethod
    def get_init_psi_t(self) -> np.ndarray:
        """Initial interval-specific unexplained variance."""
        raise NotImplementedError

    @abstractmethod
    def get_init_log_mean_bias_t(self) -> np.ndarray:
        """Initial mean bias in log space."""
        raise NotImplementedError

    @abstractmethod
    def get_init_ard_u(self) -> np.ndarray:
        """Initial ARD prior precisions."""
        raise NotImplementedError


class TrivialInitialModelParametersSupplier(InitialModelParametersSupplier):
    """Trivial initial model supplier."""
    def __init__(self,
                 denoising_model_config: DenoisingModelConfig,
                 calling_config: CopyNumberCallingConfig,
                 shared_workspace: DenoisingCallingWorkspace):
        super().__init__(denoising_model_config, calling_config, shared_workspace)

    def get_init_psi_t(self) -> np.ndarray:
        return self.denoising_model_config.psi_t_scale * np.ones(
            (self.shared_workspace.num_intervals,), dtype=types.floatX)

    def get_init_log_mean_bias_t(self) -> np.ndarray:
        return np.zeros((self.shared_workspace.num_intervals,), dtype=types.floatX)

    def get_init_ard_u(self) -> np.ndarray:
        fact = self.denoising_model_config.psi_t_scale * self.denoising_model_config.init_ard_rel_unexplained_variance
        return fact * np.ones((self.denoising_model_config.max_bias_factors,), dtype=types.floatX)


class DenoisingModel(GeneralizedContinuousModel):
    """The gCNV coverage denoising model declaration (continuous RVs only; discrete posteriors are assumed
    to be given)."""
    def __init__(self,
                 denoising_model_config: DenoisingModelConfig,
                 shared_workspace: DenoisingCallingWorkspace,
                 initial_model_parameters_supplier: InitialModelParametersSupplier):
        super().__init__()
        self.shared_workspace = shared_workspace
        register_as_global = self.register_as_global
        register_as_sample_specific = self.register_as_sample_specific

        eps_mapping = denoising_model_config.mapping_error_rate

        # interval-specific unexplained variance
        psi_t = Exponential(name='psi_t', lam=1.0 / denoising_model_config.psi_t_scale,
                            shape=(shared_workspace.num_intervals,),
                            broadcastable=(False,))
        register_as_global(psi_t)

        # sample-specific unexplained variance
        psi_s = Exponential(name='psi_s', lam=1.0 / denoising_model_config.psi_s_scale,
                            shape=(shared_workspace.num_samples,),
                            broadcastable=(False,))
        register_as_sample_specific(psi_s, sample_axis=0)

        # convert "unexplained variance" to negative binomial over-dispersion
        alpha_st = tt.maximum(tt.inv(tt.exp(psi_t.dimshuffle('x', 0) + psi_s.dimshuffle(0, 'x')) - 1.0),
                              _eps)

        # interval-specific mean log bias
        log_mean_bias_t = Normal(name='log_mean_bias_t', mu=0.0, sd=denoising_model_config.log_mean_bias_std,
                                 shape=(shared_workspace.num_intervals,),
                                 broadcastable=(False,),
                                 testval=initial_model_parameters_supplier.get_init_log_mean_bias_t())
        register_as_global(log_mean_bias_t)

        # log-normal read depth centered at the global read depth
        read_depth_mu_s = (np.log(shared_workspace.global_read_depth_s)
                           - 0.5 / denoising_model_config.depth_correction_tau)
        read_depth_s = Lognormal(name='read_depth_s',
                                 mu=read_depth_mu_s,
                                 tau=denoising_model_config.depth_correction_tau,
                                 shape=(shared_workspace.num_samples,),
                                 broadcastable=(False,),
                                 testval=shared_workspace.global_read_depth_s)
        register_as_sample_specific(read_depth_s, sample_axis=0)

        # log bias modelling, starting with the log mean bias
        log_bias_st = tt.tile(log_mean_bias_t, (shared_workspace.num_samples, 1))

        if denoising_model_config.enable_bias_factors:
            # ARD prior precisions
            ard_u = HalfFlat(name='ard_u',
                             shape=(denoising_model_config.max_bias_factors,),
                             broadcastable=(False,),
                             testval=initial_model_parameters_supplier.get_init_ard_u())
            register_as_global(ard_u)

            # bias factors
            W_tu = Normal(name='W_tu', mu=0.0, tau=ard_u.dimshuffle('x', 0),
                          shape=(shared_workspace.num_intervals, denoising_model_config.max_bias_factors),
                          broadcastable=(False, False))
            register_as_global(W_tu)

            # sample-specific bias factor loadings
            z_su = Normal(name='z_su', mu=0.0, sd=1.0,
                          shape=(shared_workspace.num_samples, denoising_model_config.max_bias_factors),
                          broadcastable=(False, False))
            register_as_sample_specific(z_su, sample_axis=0)

            # add contribution to total log bias
            if denoising_model_config.disable_bias_factors_in_active_class:
                prob_silent_class_t = tt.exp(shared_workspace.log_q_tau_tk[:, 0])
                log_bias_st += (prob_silent_class_t.dimshuffle('x', 0) * tt.dot(W_tu, z_su.T).T)
            else:
                log_bias_st += tt.dot(W_tu, z_su.T).T

        # GC bias
        if denoising_model_config.enable_explicit_gc_bias_modeling:
            # sample-specific GC bias factor loadings
            z_sg = Normal(name='z_sg', mu=0.0, sd=denoising_model_config.gc_curve_sd,
                          shape=(shared_workspace.num_samples, denoising_model_config.num_gc_bins),
                          broadcastable=(False, False))
            register_as_sample_specific(z_sg, sample_axis=0)

            # add contribution to total log bias
            log_bias_st += tst.dot(shared_workspace.W_gc_tg, z_sg.T).T

        # useful expressions
        bias_st = tt.exp(log_bias_st)

        # the expected number of erroneously mapped reads
        mean_mapping_error_correction_s = eps_mapping * read_depth_s * shared_workspace.average_ploidy_s

        denoised_copy_ratio_st = shared_workspace.n_st / (read_depth_s.dimshuffle(0, 'x') * bias_st)

        Deterministic(name='denoised_copy_ratio_st', var=denoised_copy_ratio_st)

        mu_stc = ((1.0 - eps_mapping) * read_depth_s.dimshuffle(0, 'x', 'x')
                  * bias_st.dimshuffle(0, 1, 'x')
                  * shared_workspace.copy_number_values_c.dimshuffle('x', 'x', 0)
                  + mean_mapping_error_correction_s.dimshuffle(0, 'x', 'x'))

        Deterministic(name='log_copy_number_emission_stc',
                      var=commons.negative_binomial_logp(
                          mu_stc, alpha_st.dimshuffle(0, 1, 'x'), shared_workspace.n_st.dimshuffle(0, 1, 'x')))

        # n_st (observed)
        if denoising_model_config.q_c_expectation_mode == 'map':
            def _copy_number_emission_logp(_n_st):
                mu_st = ((1.0 - eps_mapping) * read_depth_s.dimshuffle(0, 'x') * bias_st
                         * shared_workspace.c_map_st + mean_mapping_error_correction_s.dimshuffle(0, 'x'))
                log_copy_number_emission_st = commons.negative_binomial_logp(
                    mu_st, alpha_st, _n_st)
                return log_copy_number_emission_st

        elif denoising_model_config.q_c_expectation_mode == 'exact':
            def _copy_number_emission_logp(_n_st):
                _log_copy_number_emission_stc = commons.negative_binomial_logp(
                    mu_stc,
                    alpha_st.dimshuffle(0, 1, 'x'),
                    _n_st.dimshuffle(0, 1, 'x'))
                log_q_c_stc = shared_workspace.log_q_c_stc
                q_c_stc = tt.exp(log_q_c_stc)
                return tt.sum(q_c_stc * (_log_copy_number_emission_stc - log_q_c_stc), axis=2)

        elif denoising_model_config.q_c_expectation_mode == 'hybrid':
            def _copy_number_emission_logp(_n_st):
                active_class_bitmask_t = self.shared_workspace.active_class_bitmask_t
                active_class_indices = active_class_bitmask_t.nonzero()[0]
                silent_class_indices = (1 - active_class_bitmask_t).nonzero()[0]

                # for CNV-active classes, calculate exact posterior expectation
                mu_active_stc = ((1.0 - eps_mapping) * read_depth_s.dimshuffle(0, 'x', 'x')
                                 * bias_st.dimshuffle(0, 1, 'x')[:, active_class_indices, :]
                                 * shared_workspace.copy_number_values_c.dimshuffle('x', 'x', 0)
                                 + mean_mapping_error_correction_s.dimshuffle(0, 'x', 'x'))
                alpha_active_stc = tt.maximum(tt.inv((tt.exp(psi_t.dimshuffle('x', 0)[:, active_class_indices]
                                                             + psi_s.dimshuffle(0, 'x')) - 1.0)).dimshuffle(0, 1, 'x'),
                                              _eps)
                n_active_stc = _n_st.dimshuffle(0, 1, 'x')[:, active_class_indices, :]
                active_class_logp_stc = commons.negative_binomial_logp(mu_active_stc, alpha_active_stc, n_active_stc)
                log_q_c_active_stc = shared_workspace.log_q_c_stc[:, active_class_indices, :]
                q_c_active_stc = tt.exp(log_q_c_active_stc)
                active_class_logp = tt.sum(q_c_active_stc * (active_class_logp_stc - log_q_c_active_stc))

                # for CNV-silent classes, use MAP copy number state
                mu_silent_st = ((1.0 - eps_mapping) * read_depth_s.dimshuffle(0, 'x') * bias_st[:, silent_class_indices]
                                * shared_workspace.c_map_st[:, silent_class_indices]
                                + mean_mapping_error_correction_s.dimshuffle(0, 'x'))
                alpha_silent_st = alpha_st[:, silent_class_indices]
                n_silent_st = _n_st[:, silent_class_indices]
                silent_class_logp = tt.sum(commons.negative_binomial_logp(mu_silent_st, alpha_silent_st, n_silent_st))

                return active_class_logp + silent_class_logp

        elif denoising_model_config.q_c_expectation_mode == 'marginalize':
            def _copy_number_emission_logp(_n_st):
                _log_copy_number_emission_stc = commons.negative_binomial_logp(
                    mu_stc,
                    alpha_st.dimshuffle(0, 1, 'x'),
                    _n_st.dimshuffle(0, 1, 'x'))
                return pm.math.logsumexp(shared_workspace.log_q_c_stc + _log_copy_number_emission_stc, axis=2)

        else:
            raise Exception("Unknown q_c expectation mode; an exception should have been raised earlier")

        DensityDist(name='n_st_obs',
                    logp=_copy_number_emission_logp,
                    observed=shared_workspace.n_st)


class CopyNumberEmissionBasicSampler:
    """Draws posterior samples from log copy number emission probabilities for a given variational
    approximation to the denoising model continuous RVs."""
    def __init__(self,
                 denoising_model_config: DenoisingModelConfig,
                 calling_config: CopyNumberCallingConfig,
                 inference_params: HybridInferenceParameters,
                 shared_workspace: DenoisingCallingWorkspace,
                 denoising_model: DenoisingModel):
        self.model_config = denoising_model_config
        self.calling_config = calling_config
        self.inference_params = inference_params
        self.shared_workspace = shared_workspace
        self.denoising_model = denoising_model
        self._simultaneous_log_copy_number_emission_sampler = None

    def update_approximation(self, approx: pm.approximations.MeanField):
        """Generates a new compiled sampler based on a given approximation.
        Args:
            approx: an instance of PyMC3 mean-field approximation

        Returns:
            None
        """
        self._simultaneous_log_copy_number_emission_sampler = \
            self._get_compiled_simultaneous_log_copy_number_emission_sampler(approx)

    @property
    def is_sampler_initialized(self):
        return self._simultaneous_log_copy_number_emission_sampler is not None

    def draw(self) -> np.ndarray:
        assert self.is_sampler_initialized, "Posterior approximation is not provided yet"
        return self._simultaneous_log_copy_number_emission_sampler()

    @th.configparser.change_flags(compute_test_value="off")
    def _get_compiled_simultaneous_log_copy_number_emission_sampler(self, approx: pm.approximations.MeanField):
        """For a given variational approximation, returns a compiled theano function that draws posterior samples
        from log copy number emission probabilities."""
        log_copy_number_emission_stc = commons.stochastic_node_mean_symbolic(
            approx, self.denoising_model['log_copy_number_emission_stc'],
            size=self.inference_params.log_emission_samples_per_round)
        return th.function(inputs=[], outputs=log_copy_number_emission_stc)


class HHMMClassAndCopyNumberBasicCaller:
    """This class updates copy number and interval class posteriors according to the following hierarchical
    hidden Markov model:

        class_prior_k --> (tau_1) --> (tau_2) --> (tau_3) --> ...
                             |           |           |
                             |           |           |
                             v           v           v
                           (c_s1) -->  (c_s2) -->  (c_s3) --> ...
                             |           |           |
                             |           |           |
                             v           v           v
                            n_s1        n_s2        n_s3

        The posterior probability of `tau` and `c_s`, q(tau) and q(c_s) respectively, are obtained via
        the following variational ansatz:

            \prod_s p(tau, c_s | n) ~ q(tau) \prod_s q(c_s),

        where correlations between intervals are preserved in both chains, however, cross-correlations
        between `tau` and `c` are neglected, including correlations induced between copy numbers of
        different samples. As usual, the posteriors are determined by minimizing the KL divergence w.r.t.
        the true posterior resulting in the following iterative scheme:

        - Given q(tau), the effective copy number prior for the first interval and the effective copy number
          transition probabilities are determined (see _get_update_copy_number_hmm_specs_compiled_function).
          Along with the given emission probabilities to sample read counts, q(c_s) is updated using the
          forward-backward algorithm for each sample (see _update_copy_number_log_posterior)

        - Given q(c_s), the emission probability of each copy number class (tau) is determined
          (see _get_update_log_class_emission_tk_theano_func). The class prior and transition probabilities
          are fixed hyperparameters. Therefore, q(tau) can be updated immediately using a single run
          of forward-backward algorithm (see _update_class_log_posterior).
    """
    CopyNumberForwardBackwardResult = collections.namedtuple(
        'CopyNumberForwardBackwardResult',
        'sample_index, new_log_posterior_tc, copy_number_update_size, log_likelihood')

    def __init__(self,
                 calling_config: CopyNumberCallingConfig,
                 inference_params: HybridInferenceParameters,
                 shared_workspace: DenoisingCallingWorkspace,
                 disable_class_update: bool,
                 temperature: types.TensorSharedVariable):
        self.calling_config = calling_config
        self.inference_params = inference_params
        self.shared_workspace = shared_workspace
        self.disable_class_update = disable_class_update
        self.temperature = temperature

        # generate the 2-class inventory of copy number priors (CNV-silent, CNV-active) for all samples
        # according to their respective germline contig ploidies
        pi_sjkc = np.zeros((shared_workspace.num_samples,
                            shared_workspace.num_contigs,
                            calling_config.num_copy_number_classes,
                            calling_config.num_copy_number_states), dtype=types.floatX)
        for si in range(shared_workspace.num_samples):
            pi_sjkc[si, :, :, :] = self.get_copy_number_prior_for_sample_jkc(
                calling_config.num_copy_number_states,
                calling_config.p_alt,
                shared_workspace.baseline_copy_number_sj[si, :])[:, :, :]
        self.pi_sjkc: types.TensorSharedVariable = th.shared(pi_sjkc, name='pi_sjkc', borrow=config.borrow_numpy)

        # compiled function for forward-backward updates of copy number posterior
        self._hmm_q_copy_number = TheanoForwardBackward(
            log_posterior_probs_output_tc=None,
            resolve_nans=False,
            do_thermalization=True,
            do_admixing=True,
            include_update_size_output=True,
            include_alpha_beta_output=False)

        if not disable_class_update:
            # compiled function for forward-backward update of class posterior
            # Note:
            #   if p_active == 0, we have to deal with inf - inf expressions properly.
            #   setting resolve_nans = True takes care of such ambiguities.
            self._hmm_q_class = TheanoForwardBackward(
                log_posterior_probs_output_tc=shared_workspace.log_q_tau_tk,
                resolve_nans=(calling_config.p_active == 0),
                do_thermalization=True,
                do_admixing=True,
                include_update_size_output=True,
                include_alpha_beta_output=False)

            # compiled function for update of class log emission
            self._update_log_class_emission_tk_theano_func = self._get_update_log_class_emission_tk_theano_func()
        else:
            self._hmm_q_class: Optional[TheanoForwardBackward] = None
            self._update_log_class_emission_tk_theano_func = None

        # compiled function for variational update of copy number HMM specs
        self._get_copy_number_hmm_specs_theano_func = self.get_compiled_copy_number_hmm_specs_theano_func()

    @staticmethod
    def get_copy_number_prior_for_sample_jkc(num_copy_number_states: int,
                                             p_alt: float,
                                             baseline_copy_number_j: np.ndarray) -> np.ndarray:
        """Returns copy-number prior probabilities for each contig (j) and class (k) as a 3d ndarray.

        Args:
            num_copy_number_states: total number of copy-number states
            p_alt: total probability of alt copy-number states
            baseline_copy_number_j: baseline copy-number state for each contig

        Returns:
            a 3d ndarray
        """
        p_baseline = 1.0 - (num_copy_number_states - 1) * p_alt
        pi_jkc = np.zeros((len(baseline_copy_number_j), 2, num_copy_number_states), dtype=types.floatX)
        for j, baseline_state in enumerate(baseline_copy_number_j):
            # the silent class
            pi_jkc[j, 0, :] = p_alt
            pi_jkc[j, 0, baseline_state] = p_baseline
            # the active class
            pi_jkc[j, 1, :] = 1.0 / num_copy_number_states

        return pi_jkc

    def call(self,
             copy_number_update_summary_statistic_reducer,
             class_update_summary_statistic_reducer) -> Tuple[np.ndarray, np.ndarray, float, float]:
        """Perform a round of update of q(tau) and q(c)

        Note:
            This function must be called until q(tau) and q(c) converge to a self-consistent solution.

        Args:
            copy_number_update_summary_statistic_reducer: a function that reduces vectors to scalars and
                is used to compile a summary of copy number posterior updates across intervals for each sample
            class_update_summary_statistic_reducer: a function that reduces vectors to scalars and
                is used to compile a summary of interval class posterior updates across intervals

        Returns:
            copy number update summary (ndarray of size `num_samples`),
            copy number Markov chain log likelihoods (ndarray of size `num_samples`),
            interval class update summary,
            interval class Markov chain log likelihood
        """
        # copy number posterior update
        copy_number_update_s, copy_number_log_likelihoods_s = self._update_copy_number_log_posterior(
            copy_number_update_summary_statistic_reducer)

        if not self.disable_class_update:
            # class posterior update
            self._update_log_class_emission_tk()
            class_update, class_log_likelihood = self._update_class_log_posterior(
                class_update_summary_statistic_reducer)
        else:
            class_update = None
            class_log_likelihood = None

        return copy_number_update_s, copy_number_log_likelihoods_s, class_update, class_log_likelihood

    def _update_copy_number_log_posterior(self, copy_number_update_summary_statistic_reducer) \
            -> Tuple[np.ndarray, np.ndarray]:
        ws = self.shared_workspace
        copy_number_update_s = np.zeros((ws.num_samples,), dtype=types.floatX)
        copy_number_log_likelihoods_s = np.zeros((ws.num_samples,), dtype=types.floatX)
        num_calling_processes = self.calling_config.num_calling_processes

        def _run_single_sample_fb(_sample_index: int):
            # step 1. calculate copy-number HMM log prior and log transition matrix
            pi_jkc = self.pi_sjkc.get_value(borrow=True)[_sample_index, ...]
            cnv_stay_prob_t = self.shared_workspace.cnv_stay_prob_t.get_value(borrow=True)
            log_q_tau_tk = self.shared_workspace.log_q_tau_tk.get_value(borrow=True)
            t_to_j_map = self.shared_workspace.t_to_j_map.get_value(borrow=True)
            hmm_spec = self._get_copy_number_hmm_specs_theano_func(pi_jkc, cnv_stay_prob_t, log_q_tau_tk, t_to_j_map)
            log_prior_c = hmm_spec[0]
            log_trans_tcc = hmm_spec[1]

            prev_log_posterior_tc = ws.log_q_c_stc.get_value(borrow=True)[_sample_index, ...]
            log_copy_number_emission_tc = ws.log_copy_number_emission_stc.get_value(borrow=True)[_sample_index, ...]

            # step 2. run forward-backward and update copy-number posteriors
            _fb_result = self._hmm_q_copy_number.perform_forward_backward(
                log_prior_c, log_trans_tcc, log_copy_number_emission_tc,
                prev_log_posterior_tc=prev_log_posterior_tc,
                admixing_rate=self.inference_params.caller_internal_admixing_rate,
                temperature=self.temperature.get_value()[0])
            new_log_posterior_tc = _fb_result.log_posterior_probs_tc
            copy_number_update_size = copy_number_update_summary_statistic_reducer(_fb_result.update_norm_t)
            log_likelihood = float(_fb_result.log_data_likelihood)

            return self.CopyNumberForwardBackwardResult(
                _sample_index, new_log_posterior_tc, copy_number_update_size, log_likelihood)

        def _update_log_q_c_stc_inplace(log_q_c_stc, _sample_index, new_log_posterior_tc):
            log_q_c_stc[_sample_index, :, :] = new_log_posterior_tc[:, :]
            return log_q_c_stc

        max_chunks = ws.num_samples // num_calling_processes + 1
        for chunk_index in range(max_chunks):
            begin_index = chunk_index * num_calling_processes
            end_index = min((chunk_index + 1) * num_calling_processes, ws.num_samples)
            if begin_index >= ws.num_samples:
                break
            # todo multiprocessing
            # with mp.Pool(processes=num_calling_processes) as pool:
            #     for fb_result in pool.map(_run_single_sample_fb, range(begin_index, end_index)):
            for fb_result in [_run_single_sample_fb(sample_index)
                              for sample_index in range(begin_index, end_index)]:
                # update log posterior in the workspace
                ws.log_q_c_stc.set_value(
                    _update_log_q_c_stc_inplace(
                        ws.log_q_c_stc.get_value(borrow=True),
                        fb_result.sample_index, fb_result.new_log_posterior_tc),
                    borrow=True)
                # update summary stats
                copy_number_update_s[fb_result.sample_index] = fb_result.copy_number_update_size
                copy_number_log_likelihoods_s[fb_result.sample_index] = fb_result.log_likelihood

        return copy_number_update_s, copy_number_log_likelihoods_s

    def _update_log_class_emission_tk(self):
        self._update_log_class_emission_tk_theano_func()

    def _update_class_log_posterior(self, class_update_summary_statistic_reducer) -> Tuple[float, float]:
        fb_result = self._hmm_q_class.perform_forward_backward(
            self.shared_workspace.log_prior_k,
            self.shared_workspace.log_trans_tkk,
            self.shared_workspace.log_class_emission_tk.get_value(borrow=True),
            prev_log_posterior_tc=self.shared_workspace.log_q_tau_tk.get_value(borrow=True),
            admixing_rate=self.inference_params.caller_internal_admixing_rate,
            temperature=self.temperature.get_value()[0])
        class_update_size = class_update_summary_statistic_reducer(fb_result.update_norm_t)
        log_likelihood = float(fb_result.log_data_likelihood)
        return class_update_size, log_likelihood

    def update_auxiliary_vars(self):
        self.shared_workspace.update_auxiliary_vars()

    @staticmethod
    @th.configparser.change_flags(compute_test_value="off")
    def get_compiled_copy_number_hmm_specs_theano_func() -> th.compile.function_module.Function:
        """Returns a compiled function that calculates the interval-class-averaged and probability-sum-normalized
        log copy number transition matrix and log copy number prior for the first interval

        Returned theano function inputs:
            pi_jkc: a 3d tensor containing copy-number priors for each contig (j) and each class (k)
            cnv_stay_prob_t: probability of staying on the same copy-number state at interval `t`
            log_q_tau_tk: log probability of copy-number classes at interval `t`
            t_to_j_map: a mapping from interval indices (t) to contig indices (j); it is used to unpack
                `pi_jkc` to `pi_tkc` (see below)

        Returned theano function outputs:
            log_prior_c_first_interval: log probability of copy-number states for the first interval
            log_trans_tab: log transition probability matrix from interval `t` to interval `t+1`

        Note:
            In the following, we use "a" and "b" subscripts in the variable names to refer to the departure
            and destination states, respectively. Like before, "t" and "k" denote interval and class, and "j"
            refers to contig index.
        """
        # shorthands
        pi_jkc = tt.tensor3(name='pi_jkc')
        cnv_stay_prob_t = tt.vector(name='cnv_stay_prob_t')
        log_q_tau_tk = tt.matrix(name='log_q_tau_tk')
        t_to_j_map = tt.vector(name='t_to_j_map', dtype=tt.scal.uint32)

        # log prior probability for the first interval
        log_prior_c_first_interval = tt.dot(tt.log(pi_jkc[t_to_j_map[0], :, :].T), tt.exp(log_q_tau_tk[0, :]))
        log_prior_c_first_interval -= pm.logsumexp(log_prior_c_first_interval)

        # log transition matrix
        cnv_not_stay_prob_t = tt.ones_like(cnv_stay_prob_t) - cnv_stay_prob_t
        num_copy_number_states = pi_jkc.shape[2]
        delta_ab = tt.eye(num_copy_number_states)

        # map contig to interval and obtain pi_tkc for the rest of the targets
        pi_tkc = pi_jkc[t_to_j_map[1:], :, :]

        # calculate normalized log transition matrix
        # todo use logaddexp
        log_trans_tkab = tt.log(cnv_not_stay_prob_t.dimshuffle(0, 'x', 'x', 'x') * pi_tkc.dimshuffle(0, 1, 'x', 2)
                                + cnv_stay_prob_t.dimshuffle(0, 'x', 'x', 'x') * delta_ab.dimshuffle('x', 'x', 0, 1))
        q_tau_tkab = tt.exp(log_q_tau_tk[1:, :]).dimshuffle(0, 1, 'x', 'x')
        log_trans_tab = tt.sum(q_tau_tkab * log_trans_tkab, axis=1)
        log_trans_tab -= pm.logsumexp(log_trans_tab, axis=2)

        inputs = [pi_jkc, cnv_stay_prob_t, log_q_tau_tk, t_to_j_map]
        outputs = [log_prior_c_first_interval, log_trans_tab]

        return th.function(inputs=inputs, outputs=outputs)

    @th.configparser.change_flags(compute_test_value="off")
    def _get_update_log_class_emission_tk_theano_func(self) -> th.compile.function_module.Function:
        """Returns a compiled function that calculates the log interval class emission probability and
        directly updates `log_class_emission_tk` in the workspace.

        Note:
            In the following,

                xi_tab ~ posterior copy number probability of two subsequent intervals.

            We ignore correlations, i.e. we assume:

              xi_st(a, b) \equiv q_c(c_{s,t} = a, c_{s,t+1} = b)
                          \approx q_c(c_{s,t} = a) q_c(c_{s,t+1} = b)

            If needed, xi can be calculated exactly from the forward-backward tables.
        """
        # shorthands
        cnv_stay_prob_t = self.shared_workspace.cnv_stay_prob_t
        q_c_stc = tt.exp(self.shared_workspace.log_q_c_stc)
        pi_sjkc = self.pi_sjkc
        t_to_j_map = self.shared_workspace.t_to_j_map
        num_copy_number_states = self.calling_config.num_copy_number_states

        # log copy number transition matrix for each class
        cnv_not_stay_prob_t = tt.ones_like(cnv_stay_prob_t) - cnv_stay_prob_t
        delta_ab = tt.eye(num_copy_number_states)

        # calculate log class emission by reducing over samples; see below
        log_class_emission_cum_sum_tk = tt.zeros((self.shared_workspace.num_intervals - 1,
                                                  self.calling_config.num_copy_number_classes),
                                                 dtype=types.floatX)

        def inc_log_class_emission_tk_except_for_first_interval(pi_jkc, q_c_tc, cum_sum_tk):
            """Adds the contribution of a given sample to the log class emission (symbolically).

            Args:
                pi_jkc: copy number prior inventory for the sample
                q_c_tc: copy number posteriors for the sample
                cum_sum_tk: current cumulative sum of log class emission

            Returns:
                Symbolically updated cumulative sum of log class emission
            """
            # map contigs to targets (starting from the second interval)
            pi_tkc = pi_jkc[t_to_j_map[1:], :, :]

            # todo use logaddexp
            log_trans_tkab = tt.log(
                cnv_not_stay_prob_t.dimshuffle(0, 'x', 'x', 'x') * pi_tkc.dimshuffle(0, 1, 'x', 2)
                + cnv_stay_prob_t.dimshuffle(0, 'x', 'x', 'x') * delta_ab.dimshuffle('x', 'x', 0, 1))
            xi_tab = q_c_tc[:-1, :].dimshuffle(0, 1, 'x') * q_c_tc[1:, :].dimshuffle(0, 'x', 1)
            current_log_class_emission_tk = tt.sum(tt.sum(
                xi_tab.dimshuffle(0, 'x', 1, 2) * log_trans_tkab, axis=-1), axis=-1)
            return cum_sum_tk + current_log_class_emission_tk

        reduce_output = th.reduce(inc_log_class_emission_tk_except_for_first_interval,
                                  sequences=[pi_sjkc, q_c_stc],
                                  outputs_info=[log_class_emission_cum_sum_tk])
        log_class_emission_tk_except_for_first_interval = reduce_output[0]

        # the first interval
        pi_skc_first = pi_sjkc[:, t_to_j_map[0], :, :]
        q_skc_first = q_c_stc[:, 0, :].dimshuffle(0, 'x', 1)
        log_class_emission_k_first = tt.sum(tt.sum(tt.log(pi_skc_first) * q_skc_first, axis=0), axis=-1)

        # concatenate first and rest
        log_class_emission_tk = tt.concatenate((log_class_emission_k_first.dimshuffle('x', 0),
                                                log_class_emission_tk_except_for_first_interval))

        return th.function(inputs=[], outputs=[], updates=[
            (self.shared_workspace.log_class_emission_tk, log_class_emission_tk)])
