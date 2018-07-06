import logging
from typing import Callable, Optional

import numpy as np
import pymc3 as pm
import theano as th

from .inference_task_base import Sampler, Caller, CallerUpdateSummary, HybridInferenceTask, HybridInferenceParameters
from .. import config, types
from ..models.model_denoising_calling import DenoisingModel, DenoisingModelConfig, \
    CopyNumberEmissionBasicSampler, InitialModelParametersSupplier, \
    DenoisingCallingWorkspace, CopyNumberCallingConfig, HHMMClassAndCopyNumberBasicCaller

_logger = logging.getLogger(__name__)


class HHMMClassAndCopyNumberCaller(Caller):
    """This class is a wrapper around `HHMMClassAndCopyNumberBasicCaller` to be used in a cohort denoising and
    calling task."""
    def __init__(self,
                 calling_config: CopyNumberCallingConfig,
                 hybrid_inference_params: HybridInferenceParameters,
                 shared_workspace: DenoisingCallingWorkspace,
                 temperature: types.TensorSharedVariable):
        self.hybrid_inference_params = hybrid_inference_params
        self.copy_number_basic_caller = HHMMClassAndCopyNumberBasicCaller(
            calling_config, hybrid_inference_params, shared_workspace, False, temperature)
        self.shared_workspace = shared_workspace
        self.hybrid_inference_params = hybrid_inference_params
        self.log_q_c_stc_snapshot: np.ndarray = None
        self.log_q_tau_tk_snapshot: np.ndarray = None

    def snapshot(self):
        self.log_q_c_stc_snapshot = self.shared_workspace.log_q_c_stc.get_value()
        self.log_q_tau_tk_snapshot = self.shared_workspace.log_q_tau_tk.get_value()

    def call(self) -> 'HHMMClassAndCopyNumberCallerUpdateSummary':
        (copy_number_update_s, copy_number_log_likelihoods_s,
         class_update, class_log_likelihood) = self.copy_number_basic_caller.call(
            self.hybrid_inference_params.caller_summary_statistics_reducer,
            self.hybrid_inference_params.caller_summary_statistics_reducer)
        return HHMMClassAndCopyNumberCallerUpdateSummary(
            copy_number_update_s, copy_number_log_likelihoods_s,
            class_update, class_log_likelihood,
            self.hybrid_inference_params.caller_summary_statistics_reducer)

    def finalize(self):
        assert self.log_q_c_stc_snapshot is not None, "Snapshot is not taken -- forgot calling snapshot()?"
        assert self.log_q_tau_tk_snapshot is not None, "Snapshot is not taken -- forgot calling snapshot()?"
        log_q_c_stc_latest = self.shared_workspace.log_q_c_stc.get_value(borrow=True)
        log_q_tau_tk_latest = self.shared_workspace.log_q_tau_tk.get_value(borrow=True)

        # admix q_c_stc with the snapshot
        self.shared_workspace.log_q_c_stc.set_value(
            np.logaddexp(
                log_q_c_stc_latest + np.log(self.hybrid_inference_params.caller_external_admixing_rate),
                self.log_q_c_stc_snapshot + np.log(1 - self.hybrid_inference_params.caller_external_admixing_rate)),
            borrow=True)

        # admix q_tau_tk with the snapshot
        self.shared_workspace.log_q_tau_tk.set_value(
            np.logaddexp(
                log_q_tau_tk_latest + np.log(self.hybrid_inference_params.caller_external_admixing_rate),
                self.log_q_tau_tk_snapshot + np.log(1 - self.hybrid_inference_params.caller_external_admixing_rate)),
            borrow=True)

    def update_auxiliary_vars(self):
        self.copy_number_basic_caller.update_auxiliary_vars()


class HHMMClassAndCopyNumberCallerUpdateSummary(CallerUpdateSummary):
    """Copy number and interval class posterior update summary. All vector-, matrix-, and tensor- valued
    updates are reduced to scalar with a single reducer function."""
    def __init__(self,
                 copy_number_update_s: np.ndarray,
                 copy_number_log_likelihoods_s: np.ndarray,
                 class_update: float,
                 class_log_likelihood: float,
                 reducer: Callable[[np.ndarray], float]):
        self.copy_number_update_s = copy_number_update_s
        self.copy_number_log_likelihoods_s = copy_number_log_likelihoods_s
        self.class_update = class_update
        self.class_log_likelihood = class_log_likelihood
        self.copy_number_update_reduced = reducer(copy_number_update_s)

    def __repr__(self):
        return "d_q_tau: {0:2.6f}, d_q_c: {1:2.6f}".format(
            self.class_update, self.copy_number_update_reduced)

    def reduce_to_scalar(self) -> float:
        """Returns the largest of interval class update and copy number update. This scalar value
        is used for checking convergence (i.e. self-consistency between class and copy number posteriors)."""
        return np.max([self.class_update, self.copy_number_update_reduced])


class CopyNumberEmissionSampler(Sampler):
    """This class is a wrapper around `CopyNumberEmissionBasicSampler` to be used in a `HybridInferenceTask`."""
    def __init__(self,
                 hybrid_inference_params: HybridInferenceParameters,
                 denoising_config: DenoisingModelConfig,
                 calling_config: CopyNumberCallingConfig,
                 shared_workspace: DenoisingCallingWorkspace,
                 denoising_model: DenoisingModel):
        super().__init__(hybrid_inference_params)
        self.shared_workspace = shared_workspace
        self.calling_config = calling_config
        self.copy_number_emission_basic_sampler = CopyNumberEmissionBasicSampler(
            denoising_config, calling_config, hybrid_inference_params, shared_workspace, denoising_model)

    def update_approximation(self, approx: pm.approximations.MeanField):
        self.copy_number_emission_basic_sampler.update_approximation(approx)

    def draw(self) -> np.ndarray:
        return self.copy_number_emission_basic_sampler.draw()

    def reset(self):
        self.shared_workspace.log_copy_number_emission_stc.set_value(
            np.zeros((self.shared_workspace.num_samples,
                      self.shared_workspace.num_intervals,
                      self.calling_config.num_copy_number_states),
                     dtype=types.floatX), borrow=config.borrow_numpy)

    def increment(self, update):
        self.shared_workspace.log_copy_number_emission_stc.set_value(
            self.shared_workspace.log_copy_number_emission_stc.get_value(borrow=True) + update,
            borrow=True)

    def get_latest_log_emission_posterior_mean_estimator(self) -> np.ndarray:
        return self.shared_workspace.log_copy_number_emission_stc.get_value(borrow=True)


class CohortDenoisingAndCallingWarmUpTask(HybridInferenceTask):
    """Cohort denoising and calling warm-up task (no sampling/calling -- just DA-ADVI)."""
    def __init__(self,
                 denoising_config: DenoisingModelConfig,
                 hybrid_inference_params: HybridInferenceParameters,
                 shared_workspace: DenoisingCallingWorkspace,
                 initial_param_supplier: InitialModelParametersSupplier):
        _logger.info("Instantiating the denoising model (warm-up)...")
        denoising_model = DenoisingModel(denoising_config, shared_workspace, initial_param_supplier)

        elbo_normalization_factor = shared_workspace.num_samples * shared_workspace.num_intervals
        super().__init__(hybrid_inference_params, denoising_model, None, None,
                         elbo_normalization_factor=elbo_normalization_factor,
                         advi_task_name="denoising (warm-up)")

    def disengage(self):
        pass


class CohortDenoisingAndCallingMainTask(HybridInferenceTask):
    """Cohort denoising and calling inference task (w/ sampling and calling)."""
    def __init__(self,
                 denoising_config: DenoisingModelConfig,
                 calling_config: CopyNumberCallingConfig,
                 hybrid_inference_params: HybridInferenceParameters,
                 shared_workspace: DenoisingCallingWorkspace,
                 initial_param_supplier: InitialModelParametersSupplier,
                 warm_up_task: Optional[CohortDenoisingAndCallingWarmUpTask]):
        _logger.info("Instantiating the denoising model (main)...")
        denoising_model = DenoisingModel(
            denoising_config, shared_workspace, initial_param_supplier)

        if hybrid_inference_params.disable_sampler:
            copy_number_emission_sampler = None
        else:
            _logger.info("Instantiating the sampler...")
            copy_number_emission_sampler = CopyNumberEmissionSampler(
                hybrid_inference_params, denoising_config, calling_config, shared_workspace, denoising_model)

        if hybrid_inference_params.disable_caller:
            copy_number_caller = None
        else:
            _logger.info("Instantiating the copy number caller...")
            initial_temperature = hybrid_inference_params.initial_temperature
            self.temperature: types.TensorSharedVariable = th.shared(
                np.asarray([initial_temperature], dtype=types.floatX))
            copy_number_caller = HHMMClassAndCopyNumberCaller(
                calling_config, hybrid_inference_params, shared_workspace, self.temperature)

        # initialize hybrid ADVI
        elbo_normalization_factor = shared_workspace.num_samples * shared_workspace.num_intervals
        super().__init__(hybrid_inference_params, denoising_model, copy_number_emission_sampler, copy_number_caller,
                         elbo_normalization_factor=elbo_normalization_factor,
                         advi_task_name="denoising (main)",
                         calling_task_name="CNV calling")

        if warm_up_task is not None:
            _logger.info("A warm-up task was provided -- copying mean-field parameter values, temperature, "
                         "and optimizer moments from the warm-up task...")

            # temperature
            self.temperature.set_value(warm_up_task.temperature.get_value())

            # mean-field parameters
            for main_param, warm_up_param in zip(self.continuous_model_approx.params,
                                                 warm_up_task.continuous_model_approx.params):
                main_param.set_value(warm_up_param.get_value())

            # optimizer state
            self.fancy_opt.initialize_state_from_instance(warm_up_task.fancy_opt)

    def disengage(self):
        pass
