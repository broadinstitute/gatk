import logging
import os
from typing import List, Dict

import numpy as np
import pandas as pd
import pymc3 as pm

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
            {'contig_ploidy_prior_map', 'contig_set', 'num_ploidy_states'})

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
            log_q_ploidy_jk = self.ploidy_workspace.log_q_ploidy_sjk.get_value(borrow=True)[si, :, :]
            for j in range(self.ploidy_workspace.num_contigs):
                ploidy_j[j], ploidy_genotyping_quality_j[j] = model_commons.perform_genotyping(log_q_ploidy_jk[j, :])

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


def get_contig_ploidy_prior_map_from_tsv_file(input_path: str,
                                              comment=io_consts.default_comment_char,
                                              delimiter=io_consts.default_delimiter_char) -> Dict[str, np.ndarray]:
    contig_ploidy_prior_pd = pd.read_csv(input_path, delimiter=delimiter, comment=comment)
    columns = [str(x) for x in contig_ploidy_prior_pd.columns.values]
    assert len(columns) > 1
    assert columns[0] == io_consts.ploidy_prior_contig_name_column
    contig_list = [str(x) for x in contig_ploidy_prior_pd['CONTIG_NAME'].values]
    assert all([len(column) > len(io_consts.ploidy_prior_prefix)
                and column[:len(io_consts.ploidy_prior_prefix)] == io_consts.ploidy_prior_prefix
                for column in columns[1:]])
    ploidy_values = [int(column[len(io_consts.ploidy_prior_prefix):]) for column in columns[1:]]
    num_ploidy_states = np.max(ploidy_values) + 1
    contig_ploidy_prior_map: Dict[str, np.ndarray] = dict()
    for contig in contig_list:
        contig_ploidy_prior_map[contig] = np.zeros((num_ploidy_states,), dtype=types.floatX)
    for ploidy in range(num_ploidy_states):
        column_name = io_consts.ploidy_prior_prefix + str(ploidy)
        if column_name in columns:
            values = [float(x) for x in contig_ploidy_prior_pd[column_name].values]
            for j, contig in enumerate(contig_list):
                contig_ploidy_prior_map[contig][ploidy] = values[j]
    return contig_ploidy_prior_map
