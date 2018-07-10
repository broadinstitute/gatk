import logging
import os
from collections import defaultdict, OrderedDict
from typing import List, Dict

import numpy as np
import pandas as pd
import pymc3 as pm
import matplotlib.pyplot as plt
from scipy.stats import nbinom
from scipy.special import erf

from . import io_commons
from . import io_consts
from .. import types
from ..models import commons as model_commons
from ..models.model_ploidy import PloidyModelConfig
from ..models.model_ploidy import PloidyWorkspace, PloidyModel
from ..structs.metadata import SampleReadDepthMetadata, SamplePloidyMetadata

_logger = logging.getLogger(__name__)


class PloidyModelWriter:
    """Writes global ploidy model parameters to disk."""
    def __init__(self,
                 ploidy_config: PloidyModelConfig,
                 ploidy_workspace: PloidyWorkspace,
                 ploidy_model: PloidyModel,
                 ploidy_model_approx: pm.MeanField,
                 output_path: str):
        io_commons.assert_output_path_writable(output_path)
        self.ploidy_config = ploidy_config
        self.ploidy_workspace = ploidy_workspace
        self.ploidy_model = ploidy_model
        self.ploidy_model_approx = ploidy_model_approx
        self.output_path = output_path
        (self._approx_var_set, self._approx_mu_map,
         self._approx_std_map) = io_commons.extract_mean_field_posterior_parameters(self.ploidy_model_approx)

    def __call__(self):
        # write gcnvkernel version
        io_commons.write_gcnvkernel_version(self.output_path)

        # write ploidy config
        io_commons.write_dict_to_json_file(
            os.path.join(self.output_path, io_consts.default_ploidy_config_json_filename),
            self.ploidy_config.__dict__,
            {'ploidy_state_priors_map'})

        # write global variables in the posterior
        for var_name in self.ploidy_model.global_var_registry:
            assert var_name in self._approx_var_set, \
                "a variable named {0} does not exist in the approximation".format(var_name)
            _logger.info("Writing {0}...".format(var_name))
            var_mu = self._approx_mu_map[var_name]
            var_std = self._approx_std_map[var_name]
            var_mu_out_path = os.path.join(self.output_path, 'mu_' + var_name + '.tsv')
            io_commons.write_ndarray_to_tsv(var_mu_out_path, var_mu)
            var_std_out_path = os.path.join(self.output_path, 'std_' + var_name + '.tsv')
            io_commons.write_ndarray_to_tsv(var_std_out_path, var_std)


class SamplePloidyWriter:
    """Writes sample-specific ploidy model parameters and associated workspace variables to disk."""
    def __init__(self,
                 ploidy_config: PloidyModelConfig,
                 ploidy_workspace: PloidyWorkspace,
                 ploidy_model: PloidyModel,
                 ploidy_model_approx: pm.MeanField,
                 output_path: str):
        self.ploidy_config = ploidy_config
        self.ploidy_workspace = ploidy_workspace
        self.ploidy_model = ploidy_model
        self.ploidy_model_approx = ploidy_model_approx
        self.output_path = output_path
        (self._approx_var_set, self._approx_mu_map,
         self._approx_std_map) = io_commons.extract_mean_field_posterior_parameters(self.ploidy_model_approx)

    @staticmethod
    def _write_sample_contig_ploidy(sample_posterior_path: str,
                                    sample_ploidy_metadata: SamplePloidyMetadata,
                                    extra_comment_lines: List[str] = None,
                                    comment=io_consts.default_comment_char,
                                    delimiter=io_consts.default_delimiter_char):
        with open(os.path.join(sample_posterior_path, io_consts.default_sample_contig_ploidy_tsv_filename), 'w') as f:
            if extra_comment_lines is not None:
                for comment_line in extra_comment_lines:
                    f.write(comment + comment_line + '\n')
            header = delimiter.join([io_consts.contig_column_name,
                                     io_consts.ploidy_column_name,
                                     io_consts.ploidy_gq_column_name])
            f.write(header + '\n')
            for j, contig in enumerate(sample_ploidy_metadata.contig_list):
                f.write(delimiter.join([contig,
                                        repr(sample_ploidy_metadata.ploidy_j[j]),
                                        repr(sample_ploidy_metadata.ploidy_genotyping_quality_j[j])]) + '\n')

    @staticmethod
    def _write_sample_read_depth(sample_posterior_path: str,
                                 sample_read_depth_metadata: SampleReadDepthMetadata,
                                 extra_comment_lines: List[str] = None,
                                 comment=io_consts.default_comment_char,
                                 delimiter=io_consts.default_delimiter_char):
        with open(os.path.join(sample_posterior_path, io_consts.default_sample_read_depth_tsv_filename), 'w') as f:
            if extra_comment_lines is not None:
                for comment_line in extra_comment_lines:
                    f.write(comment + comment_line + '\n')
            header = delimiter.join([io_consts.global_read_depth_column_name,
                                     io_consts.average_ploidy_column_name])
            f.write(header + '\n')
            f.write(delimiter.join([repr(sample_read_depth_metadata.global_read_depth),
                                    repr(sample_read_depth_metadata.average_ploidy)]) + '\n')

    @staticmethod
    def _write_sample_ploidy_plot(sample_posterior_path: str,
                                  ploidy_workspace: PloidyWorkspace,
                                  sample_ploidy_metadata: SamplePloidyMetadata,
                                  sample_index: int):
        fig, axarr = plt.subplots(2, 1, figsize=(12, 8), gridspec_kw = {'height_ratios':[3, 1]})
        si = sample_index
        num_counts = ploidy_workspace.num_counts
        counts_m = ploidy_workspace.counts_m
        for i, contig_tuple in enumerate(ploidy_workspace.contig_tuples):
            for contig in contig_tuple:
                j = ploidy_workspace.contig_to_index_map[contig]
                hist_mask_m = np.logical_not(ploidy_workspace.hist_mask_sjm[si, j])

                hist_norm_m = ploidy_workspace.hist_sjm[si, j] / \
                              np.sum(ploidy_workspace.hist_sjm[si, j] * ploidy_workspace.hist_mask_sjm[si, j])
                axarr[0].semilogy(counts_m, hist_norm_m, c='k', alpha=0.25,
                                  label='masked data' if j == 0 else None)
                axarr[0].semilogy(counts_m, np.ma.array(hist_norm_m, mask=hist_mask_m), c='b', alpha=0.5,
                                  label='data' if j == 0 else None)
                mu = ploidy_workspace.fit_mu_sj[si, j]
                alpha = ploidy_workspace.fit_alpha_sj[si, j]
                tau = alpha / (mu * (alpha + mu))
                hist_norm = 0.5 * (1 + erf((num_counts - 0.5 - mu) * np.sqrt(tau / 2.)))
                pdf_m = nbinom.pmf(k=counts_m, n=alpha, p=alpha / (mu + alpha)) / hist_norm
                axarr[0].semilogy(counts_m, np.ma.array(pdf_m, mask=hist_mask_m), c='g', lw=2,
                                  label='histogram model' if j == 0 else None)
                axarr[0].set_xlim([0, ploidy_workspace.num_counts])
        axarr[0].set_ylim([1 / np.max(np.sum(ploidy_workspace.hist_sjm[si] * ploidy_workspace.hist_mask_sjm[si],
                                             axis=-1)),
                           1])
        axarr[0].set_title(sample_ploidy_metadata.sample_name, fontsize=16)
        axarr[0].set_xlabel('count', size=14)
        axarr[0].set_ylabel('density', size=14)
        axarr[0].legend(loc='upper right')

        d = np.mean(ploidy_workspace.ploidy_model_approx_trace['d_s'][:, si])
        b_j_norm = np.mean(ploidy_workspace.ploidy_model_approx_trace['b_j_norm'], axis=0)
        j = np.arange(ploidy_workspace.num_contigs)

        axarr[1].errorbar(j, ploidy_workspace.fit_mu_sj[si] / (d * b_j_norm),
                          yerr=ploidy_workspace.fit_mu_sd_sj[si] / (d * b_j_norm),
                          c='g', fmt='o', elinewidth=2, alpha=0.5, label='histogram model')
        axarr[1].scatter(j, sample_ploidy_metadata.ploidy_j, c='r', label='ploidy model')
        axarr[1].set_xticks(j)
        axarr[1].set_xticklabels(ploidy_workspace.contigs)
        axarr[1].set_xlabel('contig', size=14)
        axarr[1].set_ylabel('ploidy', size=14)
        axarr[1].set_ylim([0, ploidy_workspace.num_ploidies])
        axarr[1].legend(loc='lower left')

        fig.tight_layout(pad=0.5)
        fig.savefig(os.path.join(sample_posterior_path, io_consts.default_sample_ploidy_plot_svg_filename))
        plt.close()


    def __call__(self):
        for si, sample_name in enumerate(self.ploidy_workspace.sample_names):
            sample_name_comment_line = [io_consts.sample_name_sam_header_prefix + sample_name]
            sample_posterior_path = os.path.join(self.output_path, io_consts.sample_folder_prefix + repr(si))
            io_commons.assert_output_path_writable(sample_posterior_path, try_creating_output_path=True)
            _logger.info("Saving posteriors for sample \"{0}\" in \"{1}\"...".format(
                sample_name, sample_posterior_path))

            # find best contig ploidy calls and calculate ploidy genotyping quality
            ploidy_j = np.zeros((self.ploidy_workspace.num_contigs,), dtype=types.small_uint)
            ploidy_genotyping_quality_j = np.zeros((self.ploidy_workspace.num_contigs,), dtype=types.floatX)
            log_q_ploidy_jl = self.ploidy_workspace.log_q_ploidy_sjl[si, :, :]
            for j in range(self.ploidy_workspace.num_contigs):
                ploidy_j[j], ploidy_genotyping_quality_j[j] = model_commons.perform_genotyping(log_q_ploidy_jl[j, :])

            # generate sample ploidy metadata
            sample_ploidy_metadata = SamplePloidyMetadata(
                sample_name, ploidy_j, ploidy_genotyping_quality_j,
                self.ploidy_workspace.interval_list_metadata.ordered_contig_list)

            # generate sample read depth metadata
            sample_read_depth_metadata = SampleReadDepthMetadata.generate_sample_read_depth_metadata(
                self.ploidy_workspace.sample_metadata_collection.get_sample_coverage_metadata(sample_name),
                sample_ploidy_metadata,
                self.ploidy_workspace.interval_list_metadata)

            # write contig ploidy
            self._write_sample_contig_ploidy(
                sample_posterior_path, sample_ploidy_metadata, extra_comment_lines=sample_name_comment_line)

            # write read depth
            self._write_sample_read_depth(
                sample_posterior_path, sample_read_depth_metadata, extra_comment_lines=sample_name_comment_line)

            # write sample name
            io_commons.write_sample_name_to_txt_file(sample_posterior_path, sample_name)

            # write sample-specific posteriors in the approximation
            io_commons.write_mean_field_sample_specific_params(
                si, sample_posterior_path, self._approx_var_set, self._approx_mu_map, self._approx_std_map,
                self.ploidy_model, sample_name_comment_line)

            # write ploidy plot
            self._write_sample_ploidy_plot(
                sample_posterior_path, self.ploidy_workspace, sample_ploidy_metadata, si)


class PloidyModelReader:
    """Reads ploidy model parameters from disk and updates the provided approximation accordingly.

    Note:
        It is assumed that the provided model instance and approximation are compatible with the model
        parameters to be read. This has to be asserted beforehand by the CLI tool.
    """
    def __init__(self,
                 ploidy_model: PloidyModel,
                 ploidy_model_approx: pm.MeanField,
                 input_path: str):
        self.ploidy_model = ploidy_model
        self.ploidy_model_approx = ploidy_model_approx
        self.input_path = input_path

    def __call__(self):
        # check if the model is created with the same gcnvkernel version
        io_commons.check_gcnvkernel_version_from_path(self.input_path)

        # read model params
        io_commons.read_mean_field_global_params(self.input_path, self.ploidy_model_approx, self.ploidy_model)


def get_ploidy_state_priors_map_from_tsv_file(input_path: str,
                                              comment=io_consts.default_comment_char,
                                              delimiter=io_consts.default_delimiter_char) \
        -> Dict[List[str], Dict[List[int], float]]:
    """Reads the ploidy-state priors from a file.

    Returns:
        A map of the ploidy-state priors. This is a defaultdict(OrderedDict).
        The keys of the defaultdict are the contig tuples.
        The keys of the OrderedDict are the ploidy states, and
        the values of the OrderedDict are the normalized prior probabilities.
    """
    ploidy_state_priors_pd = pd.read_csv(input_path, delimiter=delimiter, comment=comment,
                                         dtype={'CONTIG_TUPLE': str,
                                                'PLOIDY_STATE': str,
                                                'RELATIVE_PROBABILITY': float})
    columns = [str(x) for x in ploidy_state_priors_pd.columns.values]
    assert columns == ['CONTIG_TUPLE', 'PLOIDY_STATE', 'RELATIVE_PROBABILITY']

    # read in the relative (unnormalized) probabilities
    raw_ploidy_state_priors_map = defaultdict(dict)

    for _, row in ploidy_state_priors_pd.iterrows():
        contig_tuple = tuple(contig for contig in row['CONTIG_TUPLE'][1:-1].split(','))
        contig_tuple_set = set(contig_tuple)
        assert len(contig_tuple_set) == len(contig_tuple), \
            "Contig tuples cannot contain duplicate contigs."
        ploidy_state = tuple(int(x) for x in row['PLOIDY_STATE'][1:-1].split(','))
        assert len(contig_tuple) == len(ploidy_state)
        relative_prob = row['RELATIVE_PROBABILITY']
        assert relative_prob > 0, \
            "Relative probabilities must be positive.  " \
            "Ploidy states with zero probability should not be included in the priors file."
        assert ploidy_state not in raw_ploidy_state_priors_map[contig_tuple], \
            "Relative probability should be specified only once for each contig-tuple--ploidy-state combination."
        raw_ploidy_state_priors_map[contig_tuple][ploidy_state] = relative_prob

    contig_set = set()
    for i, contig_tuple in enumerate(raw_ploidy_state_priors_map.keys()):
        for j, contig in enumerate(contig_tuple):
            assert contig not in contig_set, "Contig tuples must be disjoint."
            contig_set.add(contig)

    # normalize the probabilities
    ploidy_state_priors_map = defaultdict(OrderedDict)
    for contig_tuple in raw_ploidy_state_priors_map:
        normalizing_factor = sum(raw_ploidy_state_priors_map[contig_tuple].values())
        for ploidy_state, relative_prob in raw_ploidy_state_priors_map[contig_tuple].items():
            ploidy_state_priors_map[contig_tuple][ploidy_state] = \
                raw_ploidy_state_priors_map[contig_tuple][ploidy_state] / normalizing_factor

    return ploidy_state_priors_map