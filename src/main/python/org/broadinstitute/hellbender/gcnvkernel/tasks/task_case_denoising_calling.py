import logging
from typing import Callable

import numpy as np
import theano as th

from .inference_task_base import Caller, CallerUpdateSummary, HybridInferenceTask, HybridInferenceParameters
from .. import types
from ..inference.fancy_optimizers import FancyAdamax
from ..io.io_denoising_calling import DenoisingModelReader
from ..models.model_denoising_calling import DenoisingModel, DenoisingModelConfig, \
    InitialModelParametersSupplier, DenoisingCallingWorkspace, CopyNumberCallingConfig, \
    HHMMClassAndCopyNumberBasicCaller

_logger = logging.getLogger(__name__)


class HMMCopyNumberCaller(Caller):
    """This class is a wrapper around `HHMMClassAndCopyNumberBasicCaller` to be used in a case denoising and
    calling task.

    Note:
        Update of copy number class is disabled since it is considered a global part of the model and is
        set by the imported model.
    """
    def __init__(self,
                 calling_config: CopyNumberCallingConfig,
                 hybrid_inference_params: HybridInferenceParameters,
                 shared_workspace: DenoisingCallingWorkspace,
                 temperature: types.TensorSharedVariable):
        self.hybrid_inference_params = hybrid_inference_params
        self.copy_number_basic_caller = HHMMClassAndCopyNumberBasicCaller(
            calling_config, hybrid_inference_params, shared_workspace, True, temperature)
        self.shared_workspace = shared_workspace
        self.hybrid_inference_params = hybrid_inference_params
        self.log_q_c_stc_snapshot: np.ndarray = None

    def snapshot(self):
        self.log_q_c_stc_snapshot = self.shared_workspace.log_q_c_stc.get_value()

    def call(self) -> 'HMMCopyNumberCallerUpdateSummary':
        copy_number_update_s, copy_number_log_likelihoods_s, _, _ = self.copy_number_basic_caller.call(
            self.hybrid_inference_params.caller_summary_statistics_reducer,
            self.hybrid_inference_params.caller_summary_statistics_reducer)
        return HMMCopyNumberCallerUpdateSummary(
            copy_number_update_s, copy_number_log_likelihoods_s,
            self.hybrid_inference_params.caller_summary_statistics_reducer)

    def finalize(self):
        assert self.log_q_c_stc_snapshot is not None, "Snapshot is not taken -- forgot calling snapshot()?"
        log_q_c_stc_latest = self.shared_workspace.log_q_c_stc.get_value(borrow=True)

        # admix q_c_stc with the snapshot
        self.shared_workspace.log_q_c_stc.set_value(
            np.logaddexp(
                log_q_c_stc_latest + np.log(self.hybrid_inference_params.caller_external_admixing_rate),
                self.log_q_c_stc_snapshot + np.log(1 - self.hybrid_inference_params.caller_external_admixing_rate)),
            borrow=True)

    def update_auxiliary_vars(self):
        self.copy_number_basic_caller.update_auxiliary_vars()


class HMMCopyNumberCallerUpdateSummary(CallerUpdateSummary):
    def __init__(self,
                 copy_number_update_s: np.ndarray,
                 copy_number_log_likelihoods_s: np.ndarray,
                 reducer: Callable[[np.ndarray], float]):
        self.copy_number_update_s = copy_number_update_s
        self.copy_number_log_likelihoods_s = copy_number_log_likelihoods_s
        self.copy_number_update_reduced = reducer(copy_number_update_s)

    def __repr__(self):
        return "d_q_c: {0:2.6f}".format(self.copy_number_update_reduced)

    def reduce_to_scalar(self) -> float:
        return self.copy_number_update_reduced


class CaseDenoisingCallingTask(HybridInferenceTask):
    def __init__(self,
                 denoising_config: DenoisingModelConfig,
                 calling_config: CopyNumberCallingConfig,
                 hybrid_inference_params: HybridInferenceParameters,
                 shared_workspace: DenoisingCallingWorkspace,
                 initial_param_supplier: InitialModelParametersSupplier,
                 input_model_path: str):
        # the copy number emission sampler is the same as the cohort task
        from .task_cohort_denoising_calling import CopyNumberEmissionSampler

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
            copy_number_caller = HMMCopyNumberCaller(
                calling_config, hybrid_inference_params, shared_workspace, self.temperature)

        elbo_normalization_factor = shared_workspace.num_samples * shared_workspace.num_intervals

        # the optimizer is a custom adamax that only updates sample-specific model variables
        opt = FancyAdamax(learning_rate=hybrid_inference_params.learning_rate,
                          beta1=hybrid_inference_params.adamax_beta1,
                          beta2=hybrid_inference_params.adamax_beta2,
                          sample_specific_only=True)

        super().__init__(hybrid_inference_params, denoising_model, copy_number_emission_sampler, copy_number_caller,
                         elbo_normalization_factor=elbo_normalization_factor,
                         advi_task_name="denoising",
                         calling_task_name="CNV calling",
                         custom_optimizer=opt)

        _logger.info("Loading the model and updating the instantiated model and workspace...")
        DenoisingModelReader(denoising_config, calling_config, shared_workspace, denoising_model,
                             self.continuous_model_approx, input_model_path)()

    def disengage(self):
        pass

