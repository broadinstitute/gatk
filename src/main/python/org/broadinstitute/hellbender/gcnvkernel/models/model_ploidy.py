import argparse
import inspect
import logging

import numpy as np
import pymc3 as pm
import theano as th
import theano.tensor as tt
import pymc3.distributions.dist_math as pm_dist_math
from pymc3 import Deterministic, Dirichlet, Bound, Uniform, Gamma, Potential
from typing import List, Dict, Set, Tuple
from scipy.stats import nbinom
from scipy.misc import logsumexp
import matplotlib.pyplot as plt

from .fancy_model import GeneralizedContinuousModel
from .. import types
from ..structs.interval import Interval
from ..structs.metadata import IntervalListMetadata, SampleMetadataCollection
from ..tasks.inference_task_base import HybridInferenceParameters, HybridInferenceTask

_logger = logging.getLogger(__name__)
np.set_printoptions(threshold=np.inf)


class PloidyModelConfig:
    """Germline ploidy-model hyperparameters."""
    def __init__(self,
                 ploidy_state_priors_map: Dict[List[str], Dict[List[int], float]] = None,
                 ploidy_concentration_scale: float = 0.1,
                 depth_upper_bound: float = 1000.0,
                 error_rate_upper_bound: float = 0.1,
                 contig_bias_lower_bound: float = 0.1,
                 contig_bias_upper_bound: float = 2.0,
                 contig_bias_scale: float = 10.0):
        """Initializer.

        Args:
            ploidy_state_priors_map: Map of the ploidy-state priors. This is a defaultdict(OrderedDict).  The keys
                                     of the defaultdict are the contig tuples.  The keys of the OrderedDict
                                     are the ploidy states, and the values of the OrderedDict are the normalized
                                     prior probabilities.
            ploidy_concentration_scale: Scale factor for the concentration parameters of the per-contig-set
                                        Dirichlet prior on ploidy states
            depth_upper_bound: Upper bound of the uniform prior on the per-sample depth
            error_rate_upper_bound: Upper bound of the uniform prior on the error rate
            contig_bias_lower_bound: Lower bound of the Gamma prior on the per-contig bias
            contig_bias_upper_bound: Upper bound of the Gamma prior on the per-contig bias
            contig_bias_scale: Scale factor for the Gamma prior on the per-contig bias
        """
        assert ploidy_state_priors_map is not None
        self.ploidy_state_priors_map = ploidy_state_priors_map
        self.ploidy_concentration_scale = ploidy_concentration_scale
        self.depth_upper_bound = depth_upper_bound
        self.error_rate_upper_bound = error_rate_upper_bound
        self.contig_bias_lower_bound = contig_bias_lower_bound
        self.contig_bias_upper_bound = contig_bias_upper_bound
        self.contig_bias_scale = contig_bias_scale

    @staticmethod
    def expose_args(args: argparse.ArgumentParser, hide: Set[str] = None):
        """Exposes arguments of `__init__` to a given instance of `ArgumentParser`.

        Args:
            args: an instance of `ArgumentParser`
            hide: a set of arguments not to expose

        Returns:
            None
        """
        group = args.add_argument_group(title="Ploidy-model parameters")
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

        process_and_maybe_add("ploidy_concentration_scale",
                              type=float,
                              help="Scale factor for the concentration parameters of the per-contig-set "
                                   "Dirichlet prior on ploidy states",
                              default=initializer_params['ploidy_concentration_scale'].default)

        process_and_maybe_add("depth_upper_bound",
                              type=float,
                              help="Upper bound of the uniform prior on the per-sample depth",
                              default=initializer_params['depth_upper_bound'].default)

        process_and_maybe_add("error_rate_upper_bound",
                              type=float,
                              help="Upper bound of the uniform prior on the error rate",
                              default=initializer_params['error_rate_upper_bound'].default)

        process_and_maybe_add("contig_bias_lower_bound",
                              type=float,
                              help="Lower bound of the Gamma prior on the per-contig bias",
                              default=initializer_params['contig_bias_lower_bound'].default)

        process_and_maybe_add("contig_bias_upper_bound",
                              type=float,
                              help="Upper bound of the Gamma prior on the per-contig bias",
                              default=initializer_params['contig_bias_upper_bound'].default)

        process_and_maybe_add("contig_bias_scale",
                              type=float,
                              help="Scale factor for the Gamma prior on the per-contig bias",
                              default=initializer_params['contig_bias_scale'].default)

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
    epsilon: float = 1e-10

    """Workspace for storing data structures that are shared between continuous and discrete sectors
    of the germline contig ploidy model."""
    def __init__(self,
                 ploidy_config: PloidyModelConfig,
                 interval_list_metadata: IntervalListMetadata,
                 sample_names: List[str],
                 sample_metadata_collection: SampleMetadataCollection):
        self.ploidy_config = ploidy_config
        self.interval_list_metadata = interval_list_metadata
        self.sample_names = sample_names
        self.sample_metadata_collection = sample_metadata_collection

        assert sample_metadata_collection.all_samples_have_coverage_metadata(sample_names), \
            "Some samples do not have coverage metadata"

        # define useful quantities and shared tensors
        self.eps = self.epsilon

        self.num_samples: int = len(sample_names)
        self.num_contigs = interval_list_metadata.num_contigs
        self.num_counts = sample_metadata_collection.get_sample_coverage_metadata(sample_names[0]).max_count + 1

        # in the below, s = sample index, i = contig-tuple index, j = contig index,
        # k = ploidy-state index, l = ploidy index (equal to ploidy), m = count index

        # process the ploidy-state priors map
        self.contig_tuples: List[Tuple[str]] = list(self.ploidy_config.ploidy_state_priors_map.keys())
        self.num_contig_tuples = len(self.contig_tuples)
        self.ploidy_states_i_k: List[List[Tuple[int]]] = \
            [list(self.ploidy_config.ploidy_state_priors_map[contig_tuple].keys())
             for contig_tuple in self.contig_tuples]
        self.ploidy_state_priors_i_k: List[np.ndarray] = \
            [np.array(list(self.ploidy_config.ploidy_state_priors_map[contig_tuple].values()))
             for contig_tuple in self.contig_tuples]
        self.ploidy_j_k: List[np.ndarray] = []
        self.contigs: List[str] = []
        for i, contig_tuple in enumerate(self.contig_tuples):
            for j, contig in enumerate(contig_tuple):
                self.contigs.append(contig)
                self.ploidy_j_k.append(np.array([ploidy_state[j]
                                                 for ploidy_state in self.ploidy_states_i_k[i]]))

        assert set(self.contigs) == interval_list_metadata.contig_set, \
            "The set of contigs present in the coverage files must match exactly " \
            "the set of contigs present in the ploidy-state-priors file."

        self.contig_to_index_map = {contig: index for index, contig in enumerate(self.contigs)}
        self.num_ploidy_states_j = np.array([len(ploidy_k) for ploidy_k in self.ploidy_j_k])
        self.num_ploidies = np.max([np.max(ploidy_k) for ploidy_k in self.ploidy_j_k]) + 1
        self.is_ploidy_in_ploidy_state_j_kl = [np.zeros((self.num_ploidy_states_j[j], self.num_ploidies))
                                               for j in range(self.num_contigs)]
        self.ploidy_priors_jl = 1E-10 * np.ones((self.num_contigs, self.num_ploidies),
                                                dtype=types.floatX)

        for j in range(self.num_contigs):
            for k in range(self.num_ploidy_states_j[j]):
                ploidy = self.ploidy_j_k[j][k]
                self.is_ploidy_in_ploidy_state_j_kl[j][k, ploidy] = 1

        # count-distribution data
        self.hist_sjm = np.zeros((self.num_samples, self.num_contigs, self.num_counts), dtype=types.med_uint)
        for si, sample_name in enumerate(self.sample_names):
            sample_metadata = sample_metadata_collection.get_sample_coverage_metadata(sample_name)
            for j, contig in enumerate(self.contigs):
                self.hist_sjm[si, j] = sample_metadata.contig_hist_m[contig]
        self.counts_m = np.arange(self.num_counts, dtype=types.med_uint)

        # mask for count bins
        self.hist_mask_sjm = self._construct_mask(self.hist_sjm)

        average_ploidy = 2. # TODO
        self.d_s_testval = np.median(np.sum(self.hist_sjm * self.hist_mask_sjm * self.counts_m, axis=-1) /
                                     np.sum(self.hist_sjm * self.hist_mask_sjm, axis=-1), axis=-1) / average_ploidy

        self.fit_mu_sj = None
        self.fit_mu_sd_sj = None
        self.fit_alpha_sj = None
        self.fit_alpha_sd_sj = None

        self.log_q_ploidy_sjl = None


    @staticmethod
    def _get_contig_set_from_interval_list(interval_list: List[Interval]) -> Set[str]:
        return {interval.contig for interval in interval_list}

    @staticmethod
    def _construct_mask(hist_sjm):
        mask_sjm = np.full(np.shape(hist_sjm), True)
        # mask_sjm[hist_sjm < 10] = False
        # mask_sjm[:, :, :50] = False
        return mask_sjm

    def update_log_q_ploidy_sjl(self, ploidy_model_approx):
        num_samples = 1000
        trace = ploidy_model_approx.sample(num_samples)
        pi_i_sk = [np.mean(trace['pi_%d_sk' % i], axis=0)
                    for i in range(self.num_contig_tuples)]
        mu_j_sk = [np.mean(trace['mu_%d_sk' % j], axis=0)
                    for j in range(self.num_contigs)]
        p_j_sk = [self.fit_alpha_sj[:, j, np.newaxis] /
                  (mu_j_sk[j] + self.fit_alpha_sj[:, j, np.newaxis])
                  for j in range(self.num_contigs)]

        log_ploidy_emission_sjl = \
            np.transpose(
                    [logsumexp(
                        # skl
                        np.log(pi_i_sk[i][:, :, np.newaxis] * self.is_ploidy_in_ploidy_state_j_kl[j][np.newaxis, :, :] + self.eps) +
                                    # skm
                                    np.sum(self.hist_sjm[:, j, np.newaxis, :] *
                                           self.hist_mask_sjm[:, j, np.newaxis, :] *
                                           nbinom.logpmf(k=self.counts_m[np.newaxis, np.newaxis, :],
                                                         n=self.fit_alpha_sj[:, j, np.newaxis, np.newaxis],
                                                         p=p_j_sk[j][:, :, np.newaxis]),
                                           axis=-1)[:, :, np.newaxis],   # sum over m -> skl
                                   axis=1)  # logsumexp over k
                         for i, contig_tuple in enumerate(self.contig_tuples)
                         for j in map(self.contig_to_index_map.get, contig_tuple)], # jsl
                axes=(1, 0, 2)) # sjl

        self.log_q_ploidy_sjl = log_ploidy_emission_sjl - logsumexp(log_ploidy_emission_sjl, axis=2)[:, :, np.newaxis]


class HistogramModel(GeneralizedContinuousModel):
    def __init__(self,
                 ploidy_config: PloidyModelConfig,
                 ploidy_workspace: PloidyWorkspace):
        super().__init__()
        self.ploidy_config = ploidy_config
        self.ploidy_workspace = ploidy_workspace

        # shorthands
        num_samples = ploidy_workspace.num_samples
        num_contigs = ploidy_workspace.num_contigs
        num_counts = ploidy_workspace.num_counts
        counts_m = ploidy_workspace.counts_m
        hist_sjm = ploidy_workspace.hist_sjm
        hist_mask_sjm = ploidy_workspace.hist_mask_sjm
        eps = ploidy_workspace.eps

        register_as_sample_specific = self.register_as_sample_specific

        alpha_max = 1e5
        random_seed = 1
        batch_size_base = 32

        num_occurrences_sj = th.shared(np.sum(hist_sjm * hist_mask_sjm, axis=-1))

        batch_size = batch_size_base
        counts_m_batched = pm.Minibatch(counts_m, batch_size=batch_size, random_seed=random_seed)
        hist_sjm_th = th.shared(hist_sjm)
        hist_mask_sjm_th = th.shared(hist_mask_sjm)

        fit_mu_sj = Uniform('fit_mu_sj',
                            upper=num_counts,
                            shape=(num_samples, num_contigs))
        register_as_sample_specific(fit_mu_sj, sample_axis=0)

        fit_alpha_sj = Uniform('fit_alpha_sj',
                               upper=alpha_max,
                               shape=(num_samples, num_contigs))
        register_as_sample_specific(fit_alpha_sj, sample_axis=0)

        # p_sjm = tt.exp(negative_binomial_logp(mu=fit_mu_sj.dimshuffle(0, 1, 'x') + eps,
        #                                       alpha=fit_alpha_sj.dimshuffle(0, 1, 'x'),
        #                                       value=counts_m_batched[np.newaxis, np.newaxis, :],
        #                                       mask=hist_mask_sjm_th[:, :, counts_m_batched]))
        #
        # pm.Potential(name='logp_hist_sjm',
        #              var=tt.sum(poisson_logp(mu=num_occurrences_sj[:, :, np.newaxis] * p_sjm + eps,
        #                                      value=hist_sjm_th[:, :, counts_m_batched],
        #                                      mask=hist_mask_sjm_th[:, :, counts_m_batched])))

        pm.Potential(name='logp_hist_sjm',
                     var=tt.sum(hist_sjm_th[:, :, counts_m_batched] *
                                negative_binomial_logp(mu=fit_mu_sj[:, :, np.newaxis] + eps,
                                                       alpha=fit_alpha_sj[:, :, np.newaxis],
                                                       value=counts_m_batched[np.newaxis, np.newaxis, :],
                                                       mask=hist_mask_sjm_th[:, :, counts_m_batched])))

class HistogramInferenceTask(HybridInferenceTask):
    def __init__(self,
                 hybrid_inference_params: HybridInferenceParameters,
                 ploidy_config: PloidyModelConfig,
                 ploidy_workspace: PloidyWorkspace):
        _logger.info("Instantiating the histogram model...")
        self.histogram_model = HistogramModel(ploidy_config, ploidy_workspace)

        elbo_normalization_factor = ploidy_workspace.num_samples * ploidy_workspace.num_contigs
        super().__init__(hybrid_inference_params, self.histogram_model, None, None,
                         elbo_normalization_factor=elbo_normalization_factor,
                         advi_task_name="fitting histogram model")

        self.ploidy_config = ploidy_config
        self.ploidy_workspace = ploidy_workspace

    def disengage(self):
        trace = self.continuous_model_approx.sample(1000)

        fit_mu_sj = np.mean(trace['fit_mu_sj'], axis=0)
        fit_mu_sd_sj = np.std(trace['fit_mu_sj'], axis=0)
        fit_alpha_sj = np.mean(trace['fit_alpha_sj'], axis=0)
        fit_alpha_sd_sj = np.std(trace['fit_alpha_sj'], axis=0)

        print(fit_mu_sj)
        print(fit_alpha_sj)

        for s in range(self.ploidy_workspace.num_samples):
            fig, axarr = plt.subplots(2, 1, figsize=(12, 8), gridspec_kw = {'height_ratios':[3, 1]})
            for i, contig_tuple in enumerate(self.ploidy_workspace.contig_tuples):
                for contig in contig_tuple:
                    j = self.ploidy_workspace.contig_to_index_map[contig]
                    hist_mask_m = np.logical_not(self.ploidy_workspace.hist_mask_sjm[s, j])
                    counts_m = self.ploidy_workspace.counts_m
                    hist_norm_m = self.ploidy_workspace.hist_sjm[s, j] / np.sum(self.ploidy_workspace.hist_sjm[s, j] * self.ploidy_workspace.hist_mask_sjm[s, j])
                    axarr[0].semilogy(counts_m, hist_norm_m, c='k', alpha=0.25)
                    axarr[0].semilogy(counts_m, np.ma.array(hist_norm_m, mask=hist_mask_m), c='b', alpha=0.5)
                    mu = fit_mu_sj[s, j]
                    alpha = fit_alpha_sj[s, j]
                    pdf_m = nbinom.pmf(k=counts_m, n=alpha, p=alpha / (mu + alpha))
                    axarr[0].semilogy(counts_m, np.ma.array(pdf_m, mask=hist_mask_m), c='g', lw=2)
                    axarr[0].set_xlim([0, self.ploidy_workspace.num_counts])
            axarr[0].set_ylim([1 / np.max(np.sum(self.ploidy_workspace.hist_sjm[s] * self.ploidy_workspace.hist_mask_sjm[s], axis=-1)), 1E-1])
            axarr[0].set_xlabel('count', size=14)
            axarr[0].set_ylabel('density', size=14)

            j = np.arange(self.ploidy_workspace.num_contigs)

            axarr[1].set_xticks(j)
            axarr[1].set_xticklabels(self.ploidy_workspace.contigs)
            axarr[1].set_xlabel('contig', size=14)
            axarr[1].set_ylabel('ploidy', size=14)

            fig.tight_layout(pad=0.5)
            fig.savefig('/home/slee/working/gatk/test_files/plots/sample_{0}.png'.format(s))

        self.ploidy_workspace.fit_mu_sj = fit_mu_sj
        self.ploidy_workspace.fit_mu_sd_sj = fit_mu_sd_sj
        self.ploidy_workspace.fit_alpha_sj = fit_alpha_sj
        self.ploidy_workspace.fit_alpha_sd_sj = fit_alpha_sd_sj


class PloidyModel(GeneralizedContinuousModel):
    def __init__(self,
                 ploidy_config: PloidyModelConfig,
                 ploidy_workspace: PloidyWorkspace):
        super().__init__()
        self.ploidy_config = ploidy_config
        self.ploidy_workspace = ploidy_workspace

        # shorthands
        ploidy_concentration_scale = ploidy_config.ploidy_concentration_scale
        depth_upper_bound = ploidy_config.depth_upper_bound
        error_rate_upper_bound = ploidy_config.error_rate_upper_bound
        contig_bias_lower_bound = ploidy_config.contig_bias_lower_bound
        contig_bias_upper_bound = ploidy_config.contig_bias_upper_bound
        contig_bias_scale = ploidy_config.contig_bias_scale
        contig_tuples = ploidy_workspace.contig_tuples
        num_samples = ploidy_workspace.num_samples
        num_contigs = ploidy_workspace.num_contigs
        contig_to_index_map = ploidy_workspace.contig_to_index_map
        ploidy_state_priors_i_k = ploidy_workspace.ploidy_state_priors_i_k
        ploidy_j_k = ploidy_workspace.ploidy_j_k
        d_s_testval = ploidy_workspace.d_s_testval
        eps = ploidy_workspace.eps

        register_as_global = self.register_as_global
        register_as_sample_specific = self.register_as_sample_specific

        d_s = Uniform('d_s',
                      upper=depth_upper_bound,
                      shape=num_samples,
                      testval=d_s_testval)
        register_as_sample_specific(d_s, sample_axis=0)

        b_j = Bound(Gamma,
                    lower=contig_bias_lower_bound,
                    upper=contig_bias_upper_bound)('b_j',
                                                   alpha=contig_bias_scale,
                                                   beta=contig_bias_scale,
                                                   shape=num_contigs)
        register_as_global(b_j)
        b_j_norm = Deterministic('b_j_norm', var=b_j / tt.mean(b_j))

        pi_i_sk = []
        for i, contig_tuple in enumerate(contig_tuples):
            if len(ploidy_state_priors_i_k[i]) > 1:
                pi_i_sk.append(Dirichlet('pi_%d_sk' % i,
                                         a=ploidy_concentration_scale * ploidy_state_priors_i_k[i],
                                         shape=(num_samples, len(ploidy_state_priors_i_k[i])),
                                         transform=pm.distributions.transforms.t_stick_breaking(eps),
                                         testval=ploidy_state_priors_i_k[i]))
                register_as_sample_specific(pi_i_sk[i], sample_axis=0)
            else:
                pi_i_sk.append(Deterministic('pi_%d_sk' % i, var=tt.ones((num_samples, 1))))

        error_rate_js = Uniform('error_rate_js',
                                upper=error_rate_upper_bound,
                                shape=(num_contigs, num_samples))
        register_as_sample_specific(error_rate_js, sample_axis=1)

        mu_j_sk = [Deterministic('mu_%d_sk' % j,
                                 var=d_s.dimshuffle(0, 'x') * b_j_norm[j] * \
                                     (tt.maximum(ploidy_j_k[j][np.newaxis, :], error_rate_js[j][:, np.newaxis])) + eps)
                   for j in range(num_contigs)]

        logp_j_sk = [Gamma.dist(mu=self.ploidy_workspace.fit_mu_sj[:, j, np.newaxis],
                                sd=self.ploidy_workspace.fit_mu_sd_sj[:, j, np.newaxis]).logp(mu_j_sk[j])
                     for j in range(num_contigs)]

        Potential('logp', tt.sum([pm.logsumexp(tt.log(pi_i_sk[i] + eps) +
                                               tt.sum([logp_j_sk[contig_to_index_map[contig]]
                                                       for contig in contig_tuple], axis=0),
                                               axis=1)
                                  for i, contig_tuple in enumerate(contig_tuples)]))


def bound(logp, *conditions, **kwargs):
    broadcast_conditions = kwargs.get('broadcast_conditions', True)
    if broadcast_conditions:
        alltrue = pm_dist_math.alltrue_elemwise
    else:
        alltrue = pm_dist_math.alltrue_scalar
    return tt.switch(alltrue(conditions), logp, 0)

def negative_binomial_logp(mu, alpha, value, mask=True):
    return bound(pm_dist_math.binomln(value + alpha - 1, value)
                 + pm_dist_math.logpow(mu / (mu + alpha), value)
                 + pm_dist_math.logpow(alpha / (mu + alpha), alpha),
                 mu > 0, value >= 0, alpha > 0, mask)

def poisson_logp(mu, value, mask=True):
    log_prob = bound(pm_dist_math.logpow(mu, value) - mu, mu > 0, value >= 0, mask)
    # Return zero when mu and value are both zero
    return tt.switch(tt.eq(mu, 0) * tt.eq(value, 0), 0, log_prob)
