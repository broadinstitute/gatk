import logging

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

        self.ploidy_workspace = ploidy_workspace

        _logger.info("Loading the model and updating the instantiated model and workspace...")
        PloidyModelReader(self.continuous_model, self.continuous_model_approx, input_model_path)()

    def disengage(self):
        self.ploidy_workspace.update_ploidy_model_approx_trace(
            self.continuous_model_approx, self.hybrid_inference_params.log_emission_samples_per_round)
        self.ploidy_workspace.update_log_q_ploidy_sjl()