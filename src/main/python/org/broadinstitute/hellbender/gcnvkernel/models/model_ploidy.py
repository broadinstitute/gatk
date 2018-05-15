import argparse
import inspect
import logging
from typing import List, Dict, Set, Tuple

import numpy as np
import pymc3 as pm
import theano as th
import theano.tensor as tt
from pymc3 import Normal, Deterministic, DensityDist, Bound, Exponential

from . import commons
from .fancy_model import GeneralizedContinuousModel
from .. import config, types
from ..structs.interval import Interval
from ..structs.metadata import IntervalListMetadata, SampleMetadataCollection
from ..tasks.inference_task_base import HybridInferenceParameters

_logger = logging.getLogger(__name__)


class PloidyModelConfig:
    """Germline contig ploidy model hyper-parameters."""
    def __init__(self,
                 contig_ploidy_prior_map: Dict[str, np.ndarray] = None,
                 mean_bias_sd: float = 1e-2,
                 psi_j_scale: float = 1e-3,
                 psi_s_scale: float = 1e-4,
                 mapping_error_rate: float = 1e-2):
        """Initializer.

        Args:
            contig_ploidy_prior_map: map from contigs to prior probabilities of each ploidy state
            mean_bias_sd: standard deviation of mean contig-level coverage bias
            psi_j_scale: typical scale of contig-specific unexplained variance
            psi_s_scale: typical scale of sample-specific unexplained variance
            mapping_error_rate: typical mapping error probability
        """
        assert contig_ploidy_prior_map is not None
        self.mean_bias_sd = mean_bias_sd
        self.psi_j_scale = psi_j_scale
        self.psi_s_scale = psi_s_scale
        self.mapping_error_rate = mapping_error_rate
        self.contig_ploidy_prior_map, self.num_ploidy_states = self._get_validated_contig_ploidy_prior_map(
            contig_ploidy_prior_map)
        self.contig_set = set(contig_ploidy_prior_map.keys())
        self.unordered_contig_list = list(self.contig_set)

    @staticmethod
    def _get_validated_contig_ploidy_prior_map(given_contig_ploidy_prior_map: Dict[str, np.ndarray],
                                               min_prob: float = 0) -> Tuple[Dict[str, np.ndarray], int]:
        given_contigs = set(given_contig_ploidy_prior_map.keys())
        num_ploidy_states: int = 0
        for contig in given_contigs:
            num_ploidy_states = max(num_ploidy_states, given_contig_ploidy_prior_map[contig].size)
        validated_contig_ploidy_prior_map: Dict[str, np.ndarray] = dict()
        for contig in given_contigs:
            padded_validated_prior = np.zeros((num_ploidy_states,), dtype=types.floatX) + min_prob
            given_prior = given_contig_ploidy_prior_map[contig].flatten()
            padded_validated_prior[:given_prior.size] = padded_validated_prior[:given_prior.size] + given_prior
            padded_validated_prior = commons.get_normalized_prob_vector(padded_validated_prior, config.prob_sum_tol)
            validated_contig_ploidy_prior_map[contig] = padded_validated_prior
        return validated_contig_ploidy_prior_map, num_ploidy_states

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

        initializer_params = inspect.signature(PloidyModelConfig.__init__).parameters
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

        process_and_maybe_add("mean_bias_sd",
                              type=float,
                              help="Contig-level mean coverage bias standard deviation",
                              default=initializer_params['mean_bias_sd'].default)

        process_and_maybe_add("mapping_error_rate",
                              type=float,
                              help="Typical mapping error rate",
                              default=initializer_params['mapping_error_rate'].default)

        process_and_maybe_add("psi_j_scale",
                              type=float,
                              help="Typical scale of contig-specific unexplained coverage variance",
                              default=initializer_params['psi_j_scale'].default)

        process_and_maybe_add("psi_s_scale",
                              type=float,
                              help="Typical scale of sample-specific unexplained coverage variance",
                              default=initializer_params['psi_s_scale'].default)

    @staticmethod
    def from_args_dict(args_dict: Dict):
        """Initialize an instance of `PloidyModelConfig` from a dictionary of arguments.

        Args:
            args_dict: a dictionary of arguments; the keys must match argument names in
                `PloidyModelConfig.__init__`

        Returns:
            an instance of `PloidyModelConfig`
        """
        relevant_keys = set(inspect.getfullargspec(PloidyModelConfig.__init__).args)
        relevant_kwargs = {k: v for k, v in args_dict.items() if k in relevant_keys}
        return PloidyModelConfig(**relevant_kwargs)


class PloidyWorkspace:
    """Workspace for storing data structures that are shared between continuous and discrete sectors
    of the germline contig ploidy model."""
    def __init__(self,
                 ploidy_config: PloidyModelConfig,
                 interval_list_metadata: IntervalListMetadata,
                 sample_names: List[str],
                 sample_metadata_collection: SampleMetadataCollection):
        self.interval_list_metadata = interval_list_metadata
        self.sample_metadata_collection = sample_metadata_collection
        self.ploidy_config = ploidy_config
        self.num_contigs = interval_list_metadata.num_contigs
        self.sample_names = sample_names
        self.num_samples: int = len(sample_names)
        self.num_ploidy_states = ploidy_config.num_ploidy_states
        assert all([contig in ploidy_config.contig_set for contig in interval_list_metadata.contig_set]), \
            "Some contigs do not have ploidy priors"
        assert sample_metadata_collection.all_samples_have_coverage_metadata(sample_names), \
            "Some samples do not have coverage metadata"

        # number of intervals per contig as a shared theano tensor
        self.t_j: types.TensorSharedVariable = th.shared(
            interval_list_metadata.t_j.astype(types.floatX), name='t_j', borrow=config.borrow_numpy)

        # count per contig and total count as shared theano tensors
        n_sj = np.zeros((self.num_samples, self.num_contigs), dtype=types.floatX)
        n_s = np.zeros((self.num_samples,), dtype=types.floatX)
        for si, sample_name in enumerate(self.sample_names):
            sample_metadata = sample_metadata_collection.get_sample_coverage_metadata(sample_name)
            n_sj[si, :] = sample_metadata.n_j[:]
            n_s[si] = sample_metadata.n_total
        self.n_sj: types.TensorSharedVariable = th.shared(n_sj, name='n_sj', borrow=config.borrow_numpy)
        self.n_s: types.TensorSharedVariable = th.shared(n_s, name='n_s', borrow=config.borrow_numpy)

        # integer ploidy values
        int_ploidy_values_k = np.arange(0, ploidy_config.num_ploidy_states, dtype=types.small_uint)
        self.int_ploidy_values_k = th.shared(int_ploidy_values_k, name='int_ploidy_values_k',
                                             borrow=config.borrow_numpy)

        # ploidy priors
        p_ploidy_jk = np.zeros((self.num_contigs, self.ploidy_config.num_ploidy_states), dtype=types.floatX)
        for j, contig in enumerate(interval_list_metadata.ordered_contig_list):
            p_ploidy_jk[j, :] = ploidy_config.contig_ploidy_prior_map[contig][:]
        log_p_ploidy_jk = np.log(p_ploidy_jk)
        self.log_p_ploidy_jk: types.TensorSharedVariable = th.shared(log_p_ploidy_jk, name='log_p_ploidy_jk',
                                                                     borrow=config.borrow_numpy)

        # ploidy log posteriors (initial value is immaterial)
        log_q_ploidy_sjk = np.tile(log_p_ploidy_jk, (self.num_samples, 1, 1))
        self.log_q_ploidy_sjk: types.TensorSharedVariable = th.shared(
            log_q_ploidy_sjk, name='log_q_ploidy_sjk', borrow=config.borrow_numpy)

        # ploidy log emission (initial value is immaterial)
        log_ploidy_emission_sjk = np.zeros(
            (self.num_samples, self.num_contigs, ploidy_config.num_ploidy_states), dtype=types.floatX)
        self.log_ploidy_emission_sjk: types.TensorSharedVariable = th.shared(
            log_ploidy_emission_sjk, name="log_ploidy_emission_sjk", borrow=config.borrow_numpy)

        # exclusion mask; mask(j, k) = 1 - delta(j, k)
        contig_exclusion_mask_jj = (np.ones((self.num_contigs, self.num_contigs), dtype=types.small_uint)
                                    - np.eye(self.num_contigs, dtype=types.small_uint))
        self.contig_exclusion_mask_jj = th.shared(contig_exclusion_mask_jj, name='contig_exclusion_mask_jj')

    @staticmethod
    def _get_contig_set_from_interval_list(interval_list: List[Interval]) -> Set[str]:
        return {interval.contig for interval in interval_list}


class PloidyModel(GeneralizedContinuousModel):
    """Declaration of the germline contig ploidy model (continuous variables only; posterior of discrete
    variables are assumed to be known)."""

    PositiveNormal = Bound(Normal, lower=0)  # how cool is this?

    def __init__(self,
                 ploidy_config: PloidyModelConfig,
                 ploidy_workspace: PloidyWorkspace):
        super().__init__()

        # shorthands
        t_j = ploidy_workspace.t_j
        contig_exclusion_mask_jj = ploidy_workspace.contig_exclusion_mask_jj
        n_s = ploidy_workspace.n_s
        n_sj = ploidy_workspace.n_sj
        ploidy_k = ploidy_workspace.int_ploidy_values_k
        q_ploidy_sjk = tt.exp(ploidy_workspace.log_q_ploidy_sjk)
        eps = ploidy_config.mapping_error_rate

        register_as_global = self.register_as_global
        register_as_sample_specific = self.register_as_sample_specific

        # mean per-contig bias
        mean_bias_j = self.PositiveNormal('mean_bias_j',
                                          mu=1.0,
                                          sd=ploidy_config.mean_bias_sd,
                                          shape=(ploidy_workspace.num_contigs,))
        register_as_global(mean_bias_j)

        # contig coverage unexplained variance
        psi_j = Exponential(name='psi_j',
                            lam=1.0 / ploidy_config.psi_j_scale,
                            shape=(ploidy_workspace.num_contigs,))
        register_as_global(psi_j)

        # sample-specific contig unexplained variance
        psi_s = Exponential(name='psi_s',
                            lam=1.0 / ploidy_config.psi_j_scale,
                            shape=(ploidy_workspace.num_samples,))
        register_as_sample_specific(psi_s, sample_axis=0)

        # convert "unexplained variance" to negative binomial over-dispersion
        alpha_sj = tt.inv((tt.exp(psi_j.dimshuffle('x', 0) + psi_s.dimshuffle(0, 'x')) - 1.0))

        # mean ploidy per contig per sample
        mean_ploidy_sj = tt.sum(tt.exp(ploidy_workspace.log_q_ploidy_sjk)
                                * ploidy_workspace.int_ploidy_values_k.dimshuffle('x', 'x', 0), axis=2)

        # mean-field amplification coefficient per contig
        gamma_sj = mean_ploidy_sj * t_j.dimshuffle('x', 0) * mean_bias_j.dimshuffle('x', 0)

        # gamma_rest_sj \equiv sum_{j' \neq j} gamma_sj
        gamma_rest_sj = tt.dot(gamma_sj, contig_exclusion_mask_jj)

        # NB per-contig counts
        mu_num_sjk = (t_j.dimshuffle('x', 0, 'x') * mean_bias_j.dimshuffle('x', 0, 'x')
                      * ploidy_k.dimshuffle('x', 'x', 0))
        mu_den_sjk = gamma_rest_sj.dimshuffle(0, 1, 'x') + mu_num_sjk
        eps_j = eps * t_j / tt.sum(t_j)  # average number of reads erroneously mapped to contig j
        mu_sjk = ((1.0 - eps) * (mu_num_sjk / mu_den_sjk)
                  + eps_j.dimshuffle('x', 0, 'x')) * n_s.dimshuffle(0, 'x', 'x')

        def _get_logp_sjk(_n_sj):
            _logp_sjk = commons.negative_binomial_logp(
                mu_sjk,  # mean
                alpha_sj.dimshuffle(0, 1, 'x'),  # over-dispersion
                _n_sj.dimshuffle(0, 1, 'x'))  # contig counts
            return _logp_sjk

        DensityDist(name='n_sj_obs',
                    logp=lambda _n_sj: tt.sum(q_ploidy_sjk * _get_logp_sjk(_n_sj)),
                    observed=n_sj)

        # for log ploidy emission sampling
        Deterministic(name='logp_sjk', var=_get_logp_sjk(n_sj))


class PloidyEmissionBasicSampler:
    """Draws posterior samples from the ploidy log emission probability for a given variational
    approximation to the ploidy model posterior."""
    def __init__(self, ploidy_model: PloidyModel, samples_per_round: int):
        self.ploidy_model = ploidy_model
        self.samples_per_round = samples_per_round
        self._simultaneous_log_ploidy_emission_sampler = None

    def update_approximation(self, approx: pm.approximations.MeanField):
        """Generates a new compiled sampler based on a given approximation.
        Args:
            approx: an instance of PyMC3 mean-field approximation

        Returns:
            None
        """
        self._simultaneous_log_ploidy_emission_sampler = \
            self._get_compiled_simultaneous_log_ploidy_emission_sampler(approx)

    def is_sampler_initialized(self):
        return self._simultaneous_log_ploidy_emission_sampler is not None

    def draw(self) -> np.ndarray:
        return self._simultaneous_log_ploidy_emission_sampler()

    @th.configparser.change_flags(compute_test_value="off")
    def _get_compiled_simultaneous_log_ploidy_emission_sampler(self, approx: pm.approximations.MeanField):
        """For a given variational approximation, returns a compiled theano function that draws posterior samples
        from the log ploidy emission."""
        log_ploidy_emission_sjk = commons.stochastic_node_mean_symbolic(
            approx, self.ploidy_model['logp_sjk'], size=self.samples_per_round)
        return th.function(inputs=[], outputs=log_ploidy_emission_sjk)


class PloidyBasicCaller:
    """Bayesian update of germline contig ploidy log posteriors."""
    def __init__(self,
                 inference_params: HybridInferenceParameters,
                 ploidy_workspace: PloidyWorkspace):
        self.ploidy_workspace = ploidy_workspace
        self.inference_params = inference_params
        self._update_log_q_ploidy_sjk_theano_func = self._get_update_log_q_ploidy_sjk_theano_func()

    @th.configparser.change_flags(compute_test_value="off")
    def _get_update_log_q_ploidy_sjk_theano_func(self) -> th.compile.function_module.Function:
        new_log_q_ploidy_sjk = (self.ploidy_workspace.log_p_ploidy_jk.dimshuffle('x', 0, 1)
                                + self.ploidy_workspace.log_ploidy_emission_sjk)
        new_log_q_ploidy_sjk -= pm.logsumexp(new_log_q_ploidy_sjk, axis=2)
        old_log_q_ploidy_sjk = self.ploidy_workspace.log_q_ploidy_sjk
        admixed_new_log_q_ploidy_sjk = commons.safe_logaddexp(
            new_log_q_ploidy_sjk + np.log(self.inference_params.caller_external_admixing_rate),
            old_log_q_ploidy_sjk + np.log(1.0 - self.inference_params.caller_external_admixing_rate))
        update_norm_sj = commons.get_hellinger_distance(admixed_new_log_q_ploidy_sjk, old_log_q_ploidy_sjk)
        return th.function(inputs=[],
                           outputs=[update_norm_sj],
                           updates=[(self.ploidy_workspace.log_q_ploidy_sjk, admixed_new_log_q_ploidy_sjk)])

    def call(self) -> np.ndarray:
        return self._update_log_q_ploidy_sjk_theano_func()
