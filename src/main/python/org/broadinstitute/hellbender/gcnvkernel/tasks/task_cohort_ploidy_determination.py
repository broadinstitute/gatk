import logging

from .inference_task_base import HybridInferenceTask, HybridInferenceParameters
from ..models.model_ploidy import PloidyModelConfig, PloidyModel, PloidyWorkspace, HistogramInferenceTask

_logger = logging.getLogger(__name__)


class CohortPloidyInferenceTask(HybridInferenceTask):
    """Cohort germline contig ploidy determination task."""
    def __init__(self,
                 hybrid_inference_params: HybridInferenceParameters,
                 ploidy_config: PloidyModelConfig,
                 ploidy_workspace: PloidyWorkspace):
        _logger.info("Fitting histograms...")
        histogram_task = HistogramInferenceTask(hybrid_inference_params, ploidy_config, ploidy_workspace)
        histogram_task.engage()
        histogram_task.disengage()

        _logger.info("Instantiating the germline contig ploidy determination model...")
        self.ploidy_model = PloidyModel(ploidy_config, ploidy_workspace)

        elbo_normalization_factor = ploidy_workspace.num_samples * ploidy_workspace.num_contigs
        super().__init__(hybrid_inference_params, self.ploidy_model, None, None,
                         elbo_normalization_factor=elbo_normalization_factor,
                         advi_task_name="fitting ploidy model")

        self.ploidy_workspace = ploidy_workspace

    def disengage(self):
        self.ploidy_workspace.update_ploidy_model_approx_trace(
            self.continuous_model_approx, self.hybrid_inference_params.log_emission_samples_per_round)
        self.ploidy_workspace.update_log_q_ploidy_sjl()
        self.ploidy_workspace.update_read_depth_s()
