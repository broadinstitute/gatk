import logging

from .inference_task_base import HybridInferenceTask, HybridInferenceParameters
from ..inference.fancy_optimizers import FancyAdamax
from ..io.io_ploidy import PloidyModelReader
from ..models.model_ploidy import PloidyModelConfig, PloidyModel, PloidyWorkspace

_logger = logging.getLogger(__name__)


class CasePloidyInferenceTask(HybridInferenceTask):
    """Case sample ploidy inference task."""
    def __init__(self,
                 hybrid_inference_params: HybridInferenceParameters,
                 ploidy_config: PloidyModelConfig,
                 ploidy_workspace: PloidyWorkspace,
                 input_model_path: str):
        # the caller and sampler are the same as the cohort tool
        from .task_cohort_ploidy_determination import PloidyCaller, PloidyEmissionSampler

        _logger.info("Instantiating the germline contig ploidy determination model...")
        ploidy_model = PloidyModel(ploidy_config, ploidy_workspace)

        _logger.info("Instantiating the ploidy emission sampler...")
        ploidy_emission_sampler = PloidyEmissionSampler(hybrid_inference_params, ploidy_model, ploidy_workspace)

        _logger.info("Instantiating the ploidy caller...")
        ploidy_caller = PloidyCaller(hybrid_inference_params, ploidy_workspace)

        elbo_normalization_factor = ploidy_workspace.num_samples * ploidy_workspace.num_contigs

        # the optimizer is a custom adamax that only updates sample-specific model variables
        opt = FancyAdamax(learning_rate=hybrid_inference_params.learning_rate,
                          beta1=hybrid_inference_params.adamax_beta1,
                          beta2=hybrid_inference_params.adamax_beta2,
                          sample_specific_only=True)

        super().__init__(hybrid_inference_params, ploidy_model, ploidy_emission_sampler, ploidy_caller,
                         elbo_normalization_factor=elbo_normalization_factor,
                         advi_task_name="denoising",
                         calling_task_name="ploidy calling",
                         custom_optimizer=opt)

        self.ploidy_config = ploidy_config
        self.ploidy_workspace = ploidy_workspace

        _logger.info("Loading the model and updating the instantiated model and workspace...")
        PloidyModelReader(self.continuous_model, self.continuous_model_approx, input_model_path)()

    def disengage(self):
        pass
