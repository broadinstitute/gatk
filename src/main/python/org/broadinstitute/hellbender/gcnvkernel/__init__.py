import warnings

from pymc import __version__ as pymc_version

# gcnvkernel uses PyMC's variational inference internals, so it is sensitive to PyMC releases.
# MINIMUM_PYMC_VERSION is enforced; TESTED_PYMC_VERSION is the release this code is developed
# and validated against, and anything else only warns, so that environments which cannot pin
# PyMC exactly (Linux distributions, conda-forge, Homebrew, ...) are still able to run.
MINIMUM_PYMC_VERSION = "5.10.1"
TESTED_PYMC_VERSION = "5.10.1"


def _version_tuple(version: str):
    return tuple(int(part) for part in version.split(".")[:3] if part.isdigit())


if _version_tuple(pymc_version) < _version_tuple(MINIMUM_PYMC_VERSION):
    raise ImportError("gcnvkernel requires PyMC {0} or later; version found: {1}; "
                      "please upgrade the PyMC module in your python environment "
                      "accordingly.".format(MINIMUM_PYMC_VERSION, pymc_version))

if pymc_version != TESTED_PYMC_VERSION:
    warnings.warn("gcnvkernel is tested against PyMC {0}; version found: {1}. "
                  "Results are expected to be correct but have not been validated "
                  "against this version.".format(TESTED_PYMC_VERSION, pymc_version))

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
