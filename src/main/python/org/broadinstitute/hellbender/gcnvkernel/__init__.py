from pymc3 import __version__ as pymc3_version

from ._version import __version__
from .io import io_commons, io_consts, io_ploidy, io_denoising_calling, \
    io_intervals_and_counts, io_metadata, io_adamax, io_vcf_parsing, test_io_vcf_parsing
# model configs and workspaces
from .models.model_denoising_calling import CopyNumberCallingConfig, DenoisingModelConfig, DenoisingCallingWorkspace
from .models.model_denoising_calling import TrivialInitialModelParametersSupplier as DefaultDenoisingModelInitializer
from .models.model_ploidy import PloidyModelConfig, PloidyWorkspace
# post-processing
from .postprocess import test_viterbiSegmentationEngine
from .postprocess.viterbi_segmentation import ViterbiSegmentationEngine
# pre-processing and io
from .preprocess.interval_list_mask import IntervalListMask
# structs
from .structs.interval import Interval
# metadata
from .structs.metadata import IntervalListMetadata, SampleMetadataCollection, \
    SampleCoverageMetadata, SamplePloidyMetadata
# inference tasks
from .tasks.inference_task_base import HybridInferenceParameters
from .tasks.task_case_denoising_calling import CaseDenoisingCallingTask
from .tasks.task_case_ploidy_determination import CasePloidyInferenceTask
from .tasks.task_cohort_denoising_calling import CohortDenoisingAndCallingMainTask, CohortDenoisingAndCallingWarmUpTask
from .tasks.task_cohort_ploidy_determination import CohortPloidyInferenceTask
from .tasks.inference_task_base import ConvergenceError
from .utils import cli_commons, math

assert pymc3_version == "3.1", "gcnvkernel currently only supports PyMC3 3.1; version found: {0}; " \
                               "please upgrade or downgrade the PyMC3 module in your python environment " \
                               "accordingly.".format(pymc3_version)
