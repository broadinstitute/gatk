import logging
import argparse
import inspect
from typing import List, Dict, Set, Tuple

import numpy as np
import pymc3 as pm
import theano as th
import theano.tensor as tt
from pymc3 import Deterministic, DensityDist, Poisson, Gamma, Uniform, Beta, Normal, Bound, Exponential, HalfNormal
from pymc3.math import logsumexp

from ..tasks.inference_task_base import HybridInferenceParameters
from .fancy_model import GeneralizedContinuousModel
from . import commons
from .. import config, types
from ..structs.interval import Interval
from ..structs.metadata import IntervalListMetadata, SampleMetadataCollection

_logger = logging.getLogger(__name__)


class PloidyModelConfig:
    """Germline contig ploidy model hyper-parameters."""
    def __init__(self,
                 contig_ploidy_prior_map: Dict[str, np.ndarray] = None,
                 mean_bias_sd: float = 1e-2,
                 mapping_error_rate: float = 1e-2):
        """Initializer.

        Args:
            contig_ploidy_prior_map: map from contigs to prior probabilities of each ploidy state
            mean_bias_sd: standard deviation of mean contig-level coverage bias
            mapping_error_rate: typical mapping error probability
        """
        assert contig_ploidy_prior_map is not None
        self.mean_bias_sd = mean_bias_sd
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

    @staticmethod
    def _get_contig_set_from_interval_list(interval_list: List[Interval]) -> Set[str]:
        return {interval.contig for interval in interval_list}


class PloidyModel(GeneralizedContinuousModel):
    """Declaration of the germline contig ploidy model (continuous variables only; posterior of discrete
    variables are assumed to be known)."""

    def __init__(self,
                 ploidy_config: PloidyModelConfig,
                 ploidy_workspace: PloidyWorkspace):
        super().__init__()

        # shorthands
        t_j = ploidy_workspace.t_j
        n_sj = ploidy_workspace.n_sj
        ploidy_k = ploidy_workspace.int_ploidy_values_k
        q_ploidy_sjk = tt.exp(ploidy_workspace.log_q_ploidy_sjk)
        eps = ploidy_config.mapping_error_rate

        register_as_global = self.register_as_global
        register_as_sample_specific = self.register_as_sample_specific

        # per-contig bias
        bias_j = Gamma('bias_j',
                       alpha=100.0,
                       beta=100.0,
                       shape=(ploidy_workspace.num_contigs,))
        register_as_global(bias_j)
        norm_bias_j = bias_j / tt.mean(bias_j)

        # per-sample depth
        depth_s = Uniform('depth_s',
                          lower=0.0,
                          upper=10000.0,
                          shape=(ploidy_workspace.num_samples,))
        register_as_sample_specific(depth_s, sample_axis=0)

        pi_mosaicism_s = Beta(name='pi_mosaicism_s',
                              alpha=1.0,
                              beta=50.0,
                              shape=(ploidy_workspace.num_samples,))
        register_as_sample_specific(pi_mosaicism_s, sample_axis=0)
        f_mosaicism_sj = Beta(name='f_mosaicism_sj',
                              alpha=10.0,
                              beta=1.0,
                              shape=(ploidy_workspace.num_samples, ploidy_workspace.num_contigs,))
        register_as_sample_specific(f_mosaicism_sj, sample_axis=0)
        norm_f_mosaicism_sj = f_mosaicism_sj / tt.max(f_mosaicism_sj, axis=1).dimshuffle(0, 'x')

        # Poisson per-contig counts
        eps_j = HalfNormal('eps_j', sd=0.01, shape=(ploidy_workspace.num_contigs,))
        register_as_global(eps_j)
        mu_sjk =  depth_s.dimshuffle(0, 'x', 'x') * t_j.dimshuffle('x', 0, 'x') * norm_bias_j.dimshuffle('x', 0, 'x') * \
                  (ploidy_workspace.int_ploidy_values_k.dimshuffle('x', 'x', 0) + eps_j.dimshuffle('x', 0, 'x'))
        mu_mosaic_sjk = norm_f_mosaicism_sj.dimshuffle(0, 1, 'x') * mu_sjk

        # unexplained variance
        psi = Uniform(name='psi', upper=10.0)
        register_as_global(psi)

        # convert "unexplained variance" to negative binomial over-dispersion
        alpha = tt.inv((tt.exp(psi) - 1.0))

        def _get_logp_sjk(_n_sj):
            _logp_sjk = logsumexp([tt.log(1 - pi_mosaicism_s.dimshuffle(0, 'x', 'x')) + commons.negative_binomial_logp(mu_sjk, alpha.dimshuffle('x', 'x', 'x'), _n_sj.dimshuffle(0, 1, 'x')),
                                   tt.log(pi_mosaicism_s.dimshuffle(0, 'x', 'x')) + commons.negative_binomial_logp(mu_mosaic_sjk, alpha.dimshuffle('x', 'x', 'x'), _n_sj.dimshuffle(0, 1, 'x'))],
                                  axis=0)[0]
            # _logp_sjk = logsumexp([tt.log(1 - pi_mosaicism_s.dimshuffle(0, 'x', 'x')) + Poisson.dist(mu=mu_sjk).logp(_n_sj.dimshuffle(0, 1, 'x')),
            #                        tt.log(pi_mosaicism_s.dimshuffle(0, 'x', 'x')) + Poisson.dist(mu=mu_mosaic_sjk).logp(_n_sj.dimshuffle(0, 1, 'x'))],
            #                       axis=0)[0]
            return _logp_sjk

        DensityDist(name='n_sj_obs',
                    logp=lambda _n_sj: tt.sum(q_ploidy_sjk * _get_logp_sjk(_n_sj), axis=2),
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
                 ploidy_workspace: PloidyWorkspace,
                 temperature: types.TensorSharedVariable):
        self.ploidy_workspace = ploidy_workspace
        self.inference_params = inference_params
        self.temperature = temperature
        self._update_log_q_ploidy_sjk_theano_func = self._get_update_log_q_ploidy_sjk_theano_func()

    @th.configparser.change_flags(compute_test_value="off")
    def _get_update_log_q_ploidy_sjk_theano_func(self):
        new_log_q_ploidy_sjk = (self.ploidy_workspace.log_p_ploidy_jk.dimshuffle('x', 0, 1)
                                + self.ploidy_workspace.log_ploidy_emission_sjk) / self.temperature.get_value()
        new_log_q_ploidy_sjk -= pm.logsumexp(new_log_q_ploidy_sjk, axis=2)
        old_log_q_ploidy_sjk = self.ploidy_workspace.log_q_ploidy_sjk
        # admixed_new_log_q_ploidy_sjk = commons.safe_logaddexp(
        #     new_log_q_ploidy_sjk + np.log(self.inference_params.caller_admixing_rate),
        #     old_log_q_ploidy_sjk + np.log(1.0 - self.inference_params.caller_admixing_rate))
        admixed_new_log_q_ploidy_sjk = new_log_q_ploidy_sjk
        update_norm_sj = commons.get_hellinger_distance(admixed_new_log_q_ploidy_sjk, old_log_q_ploidy_sjk)
        return th.function(inputs=[],
                           outputs=[update_norm_sj],
                           updates=[(self.ploidy_workspace.log_q_ploidy_sjk, admixed_new_log_q_ploidy_sjk)])

    def call(self) -> np.ndarray:
        return self._update_log_q_ploidy_sjk_theano_func()
