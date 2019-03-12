import logging
import os
from typing import List, Optional

import numpy as np
import pandas as pd
import pymc3 as pm

from . import io_commons
from . import io_consts
from . import io_intervals_and_counts
from .. import config
from ..models import commons
from ..models.model_denoising_calling import CopyNumberCallingConfig, DenoisingModelConfig
from ..models.model_denoising_calling import DenoisingCallingWorkspace, DenoisingModel
from ..utils import math

_logger = logging.getLogger(__name__)


class DenoisingModelWriter:
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
    def _write_class_log_posterior(output_path, log_q_tau_tk):
        io_commons.write_ndarray_to_tsv(
            os.path.join(output_path, io_consts.default_class_log_posterior_tsv_filename), log_q_tau_tk)

    def __call__(self):
        # write gcnvkernel version
        io_commons.write_gcnvkernel_version(self.output_path)

        # write denoising config
        io_commons.write_dict_to_json_file(
            os.path.join(self.output_path, io_consts.default_denoising_config_json_filename),
            self.denoising_config.__dict__, set())

        # write calling config
        io_commons.write_dict_to_json_file(
            os.path.join(self.output_path, io_consts.default_calling_config_json_filename),
            self.calling_config.__dict__, set())

        # write global variables in the workspace
        self._write_class_log_posterior(
            self.output_path, self.denoising_calling_workspace.log_q_tau_tk.get_value(borrow=True))

        # write global variables in the posterior
        io_commons.write_mean_field_global_params(
            self.output_path, self.denoising_model_approx, self.denoising_model)


class DenoisingModelReader:
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

        # read global workspace variables
        self.denoising_calling_workspace.log_q_tau_tk.set_value(
            io_commons.read_ndarray_from_tsv(
                os.path.join(self.input_path, io_consts.default_class_log_posterior_tsv_filename)),
            borrow=config.borrow_numpy)

        # read global posterior parameters
        io_commons.read_mean_field_global_params(
            self.input_path, self.denoising_model_approx, self.denoising_model)


def get_sample_posterior_path(calls_path: str, sample_index: int):
    return os.path.join(calls_path, io_consts.sample_folder_prefix + repr(sample_index))


class SampleDenoisingAndCallingPosteriorsWriter:
    """Writes sample-specific model parameters and associated workspace variables to disk."""
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
    def write_ndarray_tc_with_copy_number_header(sample_posterior_path: str,
                                                 ndarray_tc: np.ndarray,
                                                 output_file_name: str,
                                                 comment=io_consts.default_comment_char,
                                                 delimiter=io_consts.default_delimiter_char,
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
        # write gcnvkernel version
        io_commons.write_gcnvkernel_version(self.output_path)

        # write denoising config
        io_commons.write_dict_to_json_file(
            os.path.join(self.output_path, io_consts.default_denoising_config_json_filename),
            self.denoising_config.__dict__, set())

        # write calling config
        io_commons.write_dict_to_json_file(
            os.path.join(self.output_path, io_consts.default_calling_config_json_filename),
            self.calling_config.__dict__, set())

        # extract mean-field parameters
        approx_var_set, approx_mu_map, approx_std_map = io_commons.extract_mean_field_posterior_parameters(
            self.denoising_model_approx)

        # compute approximate denoised copy ratios
        denoising_copy_ratios_approx_generator = commons.get_sampling_generator_for_model_approximation(
            model_approx=self.denoising_model_approx, model_var_name='denoised_copy_ratios')
        denoised_copy_ratios_mean, denoised_copy_ratios_variance =\
            math.calculate_mean_and_variance_online(denoising_copy_ratios_approx_generator)
        denoised_copy_ratios_mean = np.transpose(denoised_copy_ratios_mean)
        denoised_copy_ratios_std = np.transpose(np.sqrt(denoised_copy_ratios_variance))

        for si, sample_name in enumerate(self.denoising_calling_workspace.sample_names):
            sample_name_comment_line = [io_consts.sample_name_sam_header_prefix + sample_name]
            sample_posterior_path = get_sample_posterior_path(self.output_path, si)
            _logger.info("Saving posteriors for sample \"{0}\" in \"{1}\"...".format(
                sample_name, sample_posterior_path))
            io_commons.assert_output_path_writable(sample_posterior_path, try_creating_output_path=True)

            # write sample-specific posteriors in the approximation
            io_commons.write_mean_field_sample_specific_params(
                si, sample_posterior_path, approx_var_set, approx_mu_map, approx_std_map,
                self.denoising_model, sample_name_comment_line)

            # write sample name
            io_commons.write_sample_name_to_txt_file(sample_posterior_path, sample_name)

            # write copy number log posterior
            self.write_ndarray_tc_with_copy_number_header(
                sample_posterior_path,
                self.denoising_calling_workspace.log_q_c_stc.get_value(borrow=True)[si, ...],
                io_consts.default_copy_number_log_posterior_tsv_filename,
                extra_comment_lines=sample_name_comment_line)

            # write copy number log emission
            self.write_ndarray_tc_with_copy_number_header(
                sample_posterior_path,
                self.denoising_calling_workspace.log_copy_number_emission_stc.get_value(borrow=True)[si, ...],
                io_consts.default_copy_number_log_emission_tsv_filename,
                extra_comment_lines=sample_name_comment_line)

            # write baseline copy numbers
            baseline_copy_number_t = self.denoising_calling_workspace.baseline_copy_number_sj[
                si, self.denoising_calling_workspace.t_to_j_map.get_value(borrow=True)]
            io_commons.write_ndarray_to_tsv(
                os.path.join(sample_posterior_path, io_consts.default_baseline_copy_number_tsv_filename),
                baseline_copy_number_t,
                extra_comment_lines=sample_name_comment_line,
                header=io_consts.baseline_copy_number_column_name,
                write_shape_info=False)

            # write denoised copy ratio means
            denoised_copy_ratio_mu_s = denoised_copy_ratios_mean[:, si]
            io_commons.write_ndarray_to_tsv(
                os.path.join(sample_posterior_path, io_consts.default_denoised_copy_ratios_mean_tsv_filename),
                denoised_copy_ratio_mu_s,
                extra_comment_lines=sample_name_comment_line,
                header=io_consts.denoised_copy_ratio_mean_column_name,
                write_shape_info=False
            )

            # write denoised copy ratio standard deviations
            denoised_copy_ratio_std_s = denoised_copy_ratios_std[:, si]
            io_commons.write_ndarray_to_tsv(
                os.path.join(sample_posterior_path, io_consts.default_denoised_copy_ratios_std_tsv_filename),
                denoised_copy_ratio_std_s,
                extra_comment_lines=sample_name_comment_line,
                header=io_consts.denoised_copy_ratio_std_column_name,
                write_shape_info=False
            )


class SampleDenoisingAndCallingPosteriorsReader:
    """Reads sample-specific model parameters and associated workspace variables from disk."""
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
    def read_ndarray_tc_with_copy_number_header(sample_posterior_path: str,
                                                input_file_name: str,
                                                comment=io_consts.default_comment_char,
                                                delimiter=io_consts.default_delimiter_char) -> np.ndarray:
        """Reads a TSV-formatted dim-2 (intervals x copy-number) ndarray from a sample posterior path."""
        ndarray_tc_tsv_file = os.path.join(sample_posterior_path, input_file_name)
        ndarray_tc_pd = pd.read_csv(ndarray_tc_tsv_file, delimiter=delimiter, comment=comment)
        read_columns = [str(column_name) for column_name in ndarray_tc_pd.columns.values]
        num_read_columns = len(read_columns)
        expected_copy_number_header_columns =\
            [io_consts.copy_number_column_prefix + str(cn) for cn in range(num_read_columns)]
        assert read_columns == expected_copy_number_header_columns
        read_ndarray_tc = ndarray_tc_pd.values
        assert read_ndarray_tc.ndim == 2
        return read_ndarray_tc

    def _read_sample_copy_number_log_posterior(self,
                                               sample_posterior_path: str,
                                               comment=io_consts.default_comment_char,
                                               delimiter=io_consts.default_delimiter_char) -> np.ndarray:
        read_log_q_c_tc = self.read_ndarray_tc_with_copy_number_header(
            sample_posterior_path,
            io_consts.default_copy_number_log_posterior_tsv_filename,
            delimiter=delimiter,
            comment=comment)
        assert read_log_q_c_tc.shape == (self.denoising_calling_workspace.num_intervals,
                                             self.denoising_calling_workspace.calling_config.num_copy_number_states)
        return read_log_q_c_tc

    def _read_sample_copy_number_log_emission(self,
                                              sample_posterior_path: str,
                                              comment=io_consts.default_comment_char,
                                              delimiter=io_consts.default_delimiter_char) -> np.ndarray:
        read_log_emission_tc = self.read_ndarray_tc_with_copy_number_header(
            sample_posterior_path,
            io_consts.default_copy_number_log_emission_tsv_filename,
            delimiter=delimiter,
            comment=comment)
        assert read_log_emission_tc.shape ==\
               (self.denoising_calling_workspace.num_intervals,
                self.denoising_calling_workspace.calling_config.num_copy_number_states)
        return read_log_emission_tc

    def __call__(self):
        # assert that the interval list is the same
        interval_list_tsv_file = os.path.join(self.input_calls_path, io_consts.default_interval_list_filename)
        assert os.path.exists(interval_list_tsv_file)
        read_interval_list = io_intervals_and_counts.load_interval_list_tsv_file(interval_list_tsv_file)
        assert read_interval_list == self.denoising_calling_workspace.interval_list

        for si in range(self.denoising_calling_workspace.num_samples):
            sample_posterior_path = get_sample_posterior_path(self.input_calls_path, si)
            assert os.path.exists(sample_posterior_path)

            # read sample-specific posteriors and update approximation
            io_commons.read_mean_field_sample_specific_params(
                sample_posterior_path, si, self.denoising_calling_workspace.sample_names[si],
                self.denoising_model_approx, self.denoising_model)

            # read copy number posterior and emission and update workspace
            log_q_c_tc = self._read_sample_copy_number_log_posterior(sample_posterior_path)
            log_copy_number_emission_tc = self._read_sample_copy_number_log_emission(sample_posterior_path)

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
