from pymc3 import __version__ as pymc3_version
from ._version import __version__

# inference tasks
from .tasks.inference_task_base import HybridInferenceParameters
from .tasks.task_cohort_ploidy_determination import CohortPloidyInferenceTask
from .tasks.task_cohort_denoising_calling import CohortDenoisingAndCallingTask
from .tasks.task_case_denoising_calling import CaseDenoisingCallingTask
from .tasks.task_case_ploidy_determination import CasePloidyInferenceTask

# model configs and workspaces
from .models.model_denoising_calling import CopyNumberCallingConfig, DenoisingModelConfig, DenoisingCallingWorkspace
from .models.model_denoising_calling import TrivialInitialModelParametersSupplier as DefaultDenoisingModelInitializer
from .models.model_ploidy import PloidyModelConfig, PloidyWorkspace

# metadata
from .structs.metadata import IntervalListMetadata, SampleMetadataCollection, \
    SampleCoverageMetadata, SamplePloidyMetadata

# pre-processing and io
from .preprocess.interval_list_mask import IntervalListMask
from .io import io_commons, io_consts, io_ploidy, io_denoising_calling, \
    io_intervals_and_counts, io_metadata, io_adamax
from .utils import cli_commons

# post-processing
from .postprocess.viterbi_segmentation import ViterbiSegmentationEngine

# structs
from .structs.interval import Interval

assert pymc3_version == "3.1", "gcnvkernel currently only supports PyMC3 3.1; version found: {0}; " \
                               "please upgrade or downgrade the PyMC3 module in your python environment " \
                               "accordingly.".format(pymc3_version)
