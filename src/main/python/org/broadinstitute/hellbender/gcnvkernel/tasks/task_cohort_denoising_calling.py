import numpy as np
import pymc3 as pm
import logging
from typing import Callable
import theano as th

from .inference_task_base import Sampler, Caller, CallerUpdateSummary, HybridInferenceTask, HybridInferenceParameters
from .. import config, types
from ..models.model_denoising_calling import DenoisingModel, DenoisingModelConfig,\
    CopyNumberEmissionBasicSampler, InitialModelParametersSupplier,\
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

    def call(self) -> 'HHMMClassAndCopyNumberCallerUpdateSummary':
        (copy_number_update_s, copy_number_log_likelihoods_s,
         class_update, class_log_likelihood) = self.copy_number_basic_caller.call(
            self.hybrid_inference_params.caller_summary_statistics_reducer,
            self.hybrid_inference_params.caller_summary_statistics_reducer)
        return HHMMClassAndCopyNumberCallerUpdateSummary(
            copy_number_update_s, copy_number_log_likelihoods_s,
            class_update, class_log_likelihood,
            self.hybrid_inference_params.caller_summary_statistics_reducer)

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


class CohortDenoisingAndCallingTask(HybridInferenceTask):
    """Cohort denoising and calling inference task."""
    def __init__(self,
                 denoising_config: DenoisingModelConfig,
                 calling_config: CopyNumberCallingConfig,
                 hybrid_inference_params: HybridInferenceParameters,
                 shared_workspace: DenoisingCallingWorkspace,
                 initial_param_supplier: InitialModelParametersSupplier):
        _logger.info("Instantiating the denoising model...")
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

        elbo_normalization_factor = shared_workspace.num_samples * shared_workspace.num_intervals
        super().__init__(hybrid_inference_params, denoising_model, copy_number_emission_sampler, copy_number_caller,
                         elbo_normalization_factor=elbo_normalization_factor,
                         advi_task_name="denoising",
                         calling_task_name="CNV calling")

    def disengage(self):
        pass
