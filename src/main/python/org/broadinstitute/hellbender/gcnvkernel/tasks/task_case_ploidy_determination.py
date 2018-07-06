import logging
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import nbinom

from .inference_task_base import HybridInferenceTask, HybridInferenceParameters
from ..inference.fancy_optimizers import FancyAdamax
from ..io.io_ploidy import PloidyModelReader
from ..models.model_ploidy import PloidyModelConfig, PloidyModel, PloidyWorkspace, HistogramInferenceTask

_logger = logging.getLogger(__name__)


class CasePloidyInferenceTask(HybridInferenceTask):
    """Case sample ploidy inference task."""
    def __init__(self,
                 hybrid_inference_params: HybridInferenceParameters,
                 ploidy_config: PloidyModelConfig,
                 ploidy_workspace: PloidyWorkspace,
                 input_model_path: str):
        _logger.info("Fitting histograms...")
        histogram_task = HistogramInferenceTask(hybrid_inference_params, ploidy_config, ploidy_workspace)
        histogram_task.engage()
        histogram_task.disengage()

        _logger.info("Instantiating the germline contig ploidy determination model...")
        ploidy_model = PloidyModel(ploidy_config, ploidy_workspace)

        elbo_normalization_factor = ploidy_workspace.num_samples * ploidy_workspace.num_contigs

        # the optimizer is a custom adamax that only updates sample-specific model variables
        opt = FancyAdamax(learning_rate=hybrid_inference_params.learning_rate,
                          beta1=hybrid_inference_params.adamax_beta1,
                          beta2=hybrid_inference_params.adamax_beta2,
                          sample_specific_only=True)

        super().__init__(hybrid_inference_params, ploidy_model, None, None,
                         elbo_normalization_factor=elbo_normalization_factor,
                         advi_task_name="fitting ploidy model",
                         custom_optimizer=opt)

        self.ploidy_config = ploidy_config
        self.ploidy_workspace = ploidy_workspace

        _logger.info("Loading the model and updating the instantiated model and workspace...")
        PloidyModelReader(self.continuous_model, self.continuous_model_approx, input_model_path)()

    def disengage(self):
        self.ploidy_workspace.update_log_q_ploidy_sjl(self.continuous_model_approx)
        num_samples = 1000
        trace = self.continuous_model_approx.sample(num_samples)
        pi_i_sk = [np.mean(trace['pi_%d_sk' % i], axis=0)
                   for i in range(self.ploidy_workspace.num_contig_tuples)]
        d_s = np.mean(trace['d_s'], axis=0)
        b_j_norm = np.mean(trace['b_j_norm'], axis=0)
        error_rate_js = np.mean(trace['error_rate_js'], axis=0)
        mu_j_sk = [np.mean(trace['mu_%d_sk' % j], axis=0)
                   for j in range(self.ploidy_workspace.num_contigs)]

        fit_mu_sj = self.ploidy_workspace.fit_mu_sj
        fit_mu_sd_sj = self.ploidy_workspace.fit_mu_sd_sj
        fit_alpha_sj = self.ploidy_workspace.fit_alpha_sj

        print("pi_i_sk")
        print(pi_i_sk)
        print("d_s")
        print(d_s)
        print("b_j_norm")
        print(b_j_norm)
        q_ploidy_sjl = np.exp(self.ploidy_workspace.log_q_ploidy_sjl)
        for s, q_ploidy_jl in enumerate(q_ploidy_sjl):
            print('sample_{0}:'.format(s), np.argmax(q_ploidy_jl, axis=1))
        for s in range(self.ploidy_workspace.num_samples):
            l_j = np.argmax(q_ploidy_sjl[s], axis=1)
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

            k_j = [np.argmax(pi_i_sk[i][s])
                   for i, contig_tuple in enumerate(self.ploidy_workspace.contig_tuples)
                   for j in range(len(contig_tuple))]
            mu_j = np.array([mu_j_sk[j][s, k_j[j]] for j in range(self.ploidy_workspace.num_contigs)])

            j = np.arange(self.ploidy_workspace.num_contigs)

            axarr[1].errorbar(j, fit_mu_sj[s] / (d_s[s] * b_j_norm),
                              yerr=fit_mu_sd_sj[s] / (d_s[s] * b_j_norm),
                              c='g', fmt='o', elinewidth=2, alpha=0.5)
            axarr[1].scatter(j, l_j, c='r')
            axarr[1].set_xticks(j)
            axarr[1].set_xticklabels(self.ploidy_workspace.contigs)
            axarr[1].set_xlabel('contig', size=14)
            axarr[1].set_ylabel('ploidy', size=14)
            axarr[1].set_ylim([0, np.shape(q_ploidy_sjl)[2]])

            fig.tight_layout(pad=0.5)
            fig.savefig('/home/slee/working/gatk/test_files/plots/sample_{0}.png'.format(s))
