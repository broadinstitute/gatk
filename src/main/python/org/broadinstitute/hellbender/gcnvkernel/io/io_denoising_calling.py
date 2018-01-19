import logging

import numpy as np
import pymc3 as pm
import pandas as pd
import os
from typing import List, Optional

from .. import config
from ..models.model_denoising_calling import DenoisingCallingWorkspace, DenoisingModel
from ..models.model_denoising_calling import CopyNumberCallingConfig, DenoisingModelConfig
from . import io_commons
from . import io_consts
from . import io_intervals_and_counts

_logger = logging.getLogger(__name__)


class DenoisingModelExporter:
    """Writes global denoising model parameters to disk."""
    def __init__(self,
                 denoising_config: DenoisingModelConfig,
                 calling_config: CopyNumberCallingConfig,
                 denoising_calling_workspace: DenoisingCallingWorkspace,
                 denoising_model: DenoisingModel,
                 denoising_model_approx: pm.MeanField,
                 output_path: str):
        io_commons.assert_output_path_writable(output_path)
        self.denoising_config = denoising_config
        self.calling_config = calling_config
        self.denoising_calling_workspace = denoising_calling_workspace
        self.denoising_model = denoising_model
        self.denoising_model_approx = denoising_model_approx
        self.output_path = output_path

    @staticmethod
    def _export_class_log_posterior(output_path, log_q_tau_tk):
        io_commons.write_ndarray_to_tsv(
            os.path.join(output_path, io_consts.default_class_log_posterior_tsv_filename), log_q_tau_tk)

    def __call__(self):
        # export gcnvkernel version
        io_commons.export_gcnvkernel_version(self.output_path)

        # export denoising config
        io_commons.export_dict_to_json_file(
            os.path.join(self.output_path, io_consts.default_denoising_config_json_filename),
            self.denoising_config.__dict__, set())

        # export calling config
        io_commons.export_dict_to_json_file(
            os.path.join(self.output_path, io_consts.default_calling_config_json_filename),
            self.calling_config.__dict__, set())

        # export global variables in the workspace
        self._export_class_log_posterior(
            self.output_path, self.denoising_calling_workspace.log_q_tau_tk.get_value(borrow=True))

        # export global variables in the posterior
        io_commons.export_meanfield_global_params(
            self.output_path, self.denoising_model_approx, self.denoising_model)


class DenoisingModelImporter:
    """Reads global denoising model parameters from disk."""
    def __init__(self,
                 denoising_config: DenoisingModelConfig,
                 calling_config: CopyNumberCallingConfig,
                 denoising_calling_workspace: DenoisingCallingWorkspace,
                 denoising_model: DenoisingModel,
                 denoising_model_approx: pm.MeanField,
                 input_path: str):
        self.denoising_config = denoising_config
        self.calling_config = calling_config
        self.denoising_calling_workspace = denoising_calling_workspace
        self.denoising_model = denoising_model
        self.denoising_model_approx = denoising_model_approx
        self.input_path = input_path

    def __call__(self):
        # check if the model is created with the same gcnvkernel version
        io_commons.check_gcnvkernel_version_from_path(self.input_path)

        # import global workspace variables
        self.denoising_calling_workspace.log_q_tau_tk.set_value(
            io_commons.read_ndarray_from_tsv(
                os.path.join(self.input_path, io_consts.default_class_log_posterior_tsv_filename)),
            borrow=config.borrow_numpy)

        # import global posterior parameters
        io_commons.import_meanfield_global_params(
            self.input_path, self.denoising_model_approx, self.denoising_model)


def get_sample_posterior_path(calls_path: str, sample_index: int):
    return os.path.join(calls_path, io_consts.sample_folder_prefix + repr(sample_index))


class SampleDenoisingAndCallingPosteriorsExporter:
    """Exports sample-specific model parameters and associated workspace variables to disk."""
    def __init__(self,
                 denoising_config: DenoisingModelConfig,
                 calling_config: CopyNumberCallingConfig,
                 denoising_calling_workspace: DenoisingCallingWorkspace,
                 denoising_model: DenoisingModel,
                 denoising_model_approx: pm.MeanField,
                 output_path: str):
        io_commons.assert_output_path_writable(output_path)
        self.denoising_config = denoising_config
        self.calling_config = calling_config
        self.denoising_calling_workspace = denoising_calling_workspace
        self.denoising_model = denoising_model
        self.denoising_model_approx = denoising_model_approx
        self.output_path = output_path

    @staticmethod
    def export_ndarray_tc_with_copy_number_header(sample_posterior_path: str,
                                                  ndarray_tc: np.ndarray,
                                                  output_file_name: str,
                                                  delimiter='\t',
                                                  comment='@',
                                                  extra_comment_lines: Optional[List[str]] = None):
        assert isinstance(ndarray_tc, np.ndarray)
        assert ndarray_tc.ndim == 2
        num_copy_number_states = ndarray_tc.shape[1]
        copy_number_header_columns = [io_consts.copy_number_column_prefix + str(cn)
                                      for cn in range(num_copy_number_states)]
        with open(os.path.join(sample_posterior_path, output_file_name), 'w') as f:
            if extra_comment_lines is not None:
                for comment_line in extra_comment_lines:
                    f.write(comment + comment_line + '\n')
            f.write(delimiter.join(copy_number_header_columns) + '\n')
            for ti in range(ndarray_tc.shape[0]):
                f.write(delimiter.join([repr(x) for x in ndarray_tc[ti, :]]) + '\n')

    def __call__(self):
        # export gcnvkernel version
        io_commons.export_gcnvkernel_version(self.output_path)

        # export denoising config
        io_commons.export_dict_to_json_file(
            os.path.join(self.output_path, io_consts.default_denoising_config_json_filename),
            self.denoising_config.__dict__, set())

        # export calling config
        io_commons.export_dict_to_json_file(
            os.path.join(self.output_path, io_consts.default_calling_config_json_filename),
            self.calling_config.__dict__, set())

        # extract meanfield parameters
        approx_var_set, approx_mu_map, approx_std_map = io_commons.extract_meanfield_posterior_parameters(
            self.denoising_model_approx)

        for si, sample_name in enumerate(self.denoising_calling_workspace.sample_names):
            sample_name_comment_line = [io_consts.sample_name_sam_header_prefix + sample_name]
            sample_posterior_path = get_sample_posterior_path(self.output_path, si)
            _logger.info("Saving posteriors for sample \"{0}\" in \"{1}\"...".format(
                sample_name, sample_posterior_path))
            io_commons.assert_output_path_writable(sample_posterior_path, try_creating_output_path=True)

            # export sample-specific posteriors in the approximation
            io_commons.export_meanfield_sample_specific_params(
                si, sample_posterior_path, approx_var_set, approx_mu_map, approx_std_map,
                self.denoising_model, sample_name_comment_line)

            # export sample name
            io_commons.write_sample_name_to_txt_file(sample_posterior_path, sample_name)

            # export copy number log posterior
            self.export_ndarray_tc_with_copy_number_header(
                sample_posterior_path,
                self.denoising_calling_workspace.log_q_c_stc.get_value(borrow=True)[si, ...],
                io_consts.default_copy_number_log_posterior_tsv_filename,
                extra_comment_lines=sample_name_comment_line)

            # export copy number log emission
            self.export_ndarray_tc_with_copy_number_header(
                sample_posterior_path,
                self.denoising_calling_workspace.log_copy_number_emission_stc.get_value(borrow=True)[si, ...],
                io_consts.default_copy_number_log_emission_tsv_filename,
                extra_comment_lines=sample_name_comment_line)

            # export baseline copy numbers
            baseline_copy_number_t = self.denoising_calling_workspace.baseline_copy_number_sj[
                si, self.denoising_calling_workspace.t_to_j_map.get_value(borrow=True)]
            io_commons.write_ndarray_to_tsv(
                os.path.join(sample_posterior_path, io_consts.default_baseline_copy_number_tsv_filename),
                baseline_copy_number_t)


class SampleDenoisingAndCallingPosteriorsImporter:
    """Imports sample-specific model parameters and associated workspace variables from disk."""
    def __init__(self,
                 denoising_calling_workspace: DenoisingCallingWorkspace,
                 denoising_model: DenoisingModel,
                 denoising_model_approx: pm.MeanField,
                 input_calls_path: str):
        self.denoising_calling_workspace = denoising_calling_workspace
        self.denoising_model = denoising_model
        self.denoising_model_approx = denoising_model_approx
        self.input_calls_path = input_calls_path

    @staticmethod
    def import_ndarray_tc_with_copy_number_header(sample_posterior_path: str,
                                                  input_file_name: str,
                                                  delimiter='\t',
                                                  comment='@') -> np.ndarray:
        ndarray_tc_tsv_file = os.path.join(sample_posterior_path, input_file_name)
        ndarray_tc_pd = pd.read_csv(ndarray_tc_tsv_file, delimiter=delimiter, comment=comment)
        imported_columns = [str(column_name) for column_name in ndarray_tc_pd.columns.values]
        num_imported_columns = len(imported_columns)
        expected_copy_number_header_columns =\
            [io_consts.copy_number_column_prefix + str(cn) for cn in range(num_imported_columns)]
        assert imported_columns == expected_copy_number_header_columns
        imported_ndarray_tc = ndarray_tc_pd.values
        assert imported_ndarray_tc.ndim == 2
        return imported_ndarray_tc

    def _import_sample_copy_number_log_posterior(self,
                                                 sample_posterior_path: str,
                                                 delimiter='\t',
                                                 comment='@') -> np.ndarray:
        imported_log_q_c_tc = self.import_ndarray_tc_with_copy_number_header(
            sample_posterior_path,
            io_consts.default_copy_number_log_posterior_tsv_filename,
            delimiter=delimiter,
            comment=comment)
        assert imported_log_q_c_tc.shape == (self.denoising_calling_workspace.num_intervals,
                                             self.denoising_calling_workspace.calling_config.num_copy_number_states)
        return imported_log_q_c_tc

    def _import_sample_copy_number_log_emission(self,
                                                sample_posterior_path: str,
                                                delimiter='\t',
                                                comment='@') -> np.ndarray:
        imported_log_emission_tc = self.import_ndarray_tc_with_copy_number_header(
            sample_posterior_path,
            io_consts.default_copy_number_log_emission_tsv_filename,
            delimiter=delimiter,
            comment=comment)
        assert imported_log_emission_tc.shape ==\
               (self.denoising_calling_workspace.num_intervals,
                self.denoising_calling_workspace.calling_config.num_copy_number_states)
        return imported_log_emission_tc

    def __call__(self):
        # assert that the interval list is the same
        interval_list_tsv_file = os.path.join(self.input_calls_path, io_consts.default_interval_list_filename)
        assert os.path.exists(interval_list_tsv_file)
        imported_interval_list = io_intervals_and_counts.load_interval_list_tsv_file(interval_list_tsv_file)
        assert imported_interval_list == self.denoising_calling_workspace.interval_list

        for si in range(self.denoising_calling_workspace.num_samples):
            sample_posterior_path = get_sample_posterior_path(self.input_calls_path, si)
            assert os.path.exists(sample_posterior_path)

            # import sample-specific posteriors and update approximation
            io_commons.import_meanfield_sample_specific_params(
                sample_posterior_path, si, self.denoising_calling_workspace.sample_names[si],
                self.denoising_model_approx, self.denoising_model)

            # import copy number posterior and emission and update workspace
            log_q_c_tc = self._import_sample_copy_number_log_posterior(sample_posterior_path)
            log_copy_number_emission_tc = self._import_sample_copy_number_log_emission(sample_posterior_path)

            def update_log_q_c_stc_for_sample(log_q_c_stc):
                log_q_c_stc[si, ...] = log_q_c_tc[...]
                return log_q_c_stc

            def update_log_copy_number_emission_stc_for_sample(log_copy_number_emission_stc):
                log_copy_number_emission_stc[si, ...] = log_copy_number_emission_tc[...]
                return log_copy_number_emission_stc

            self.denoising_calling_workspace.log_q_c_stc.set_value(
                update_log_q_c_stc_for_sample(
                    self.denoising_calling_workspace.log_q_c_stc.get_value(borrow=True)),
                borrow=True)

            self.denoising_calling_workspace.log_copy_number_emission_stc.set_value(
                update_log_copy_number_emission_stc_for_sample(
                    self.denoising_calling_workspace.log_copy_number_emission_stc.get_value(borrow=True)),
                borrow=True)

        # update auxiliary workspace variables
        self.denoising_calling_workspace.update_auxiliary_vars()
