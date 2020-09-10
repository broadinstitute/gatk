"""Runtime-specific data location constants and associated helper methods.

TODO(everyone):
* Add more datasets and data types as needed.
* Update logic to determine runtime as needed.
* Update Cloud Storage and persistent disk paths as needed.
"""

import collections
from enum import Enum, auto
import os


class Runtime(Enum):
  TERRA = auto()
  ML4H_VM = auto()
  DEFAULT = auto()  # Such as an AI Platform Notebooks VM.
  # TODO(sam): PARTNERS = 4


class Dataset(Enum):
  FAKE = auto()
  UKB = auto()


class DataType(Enum):
  RESTING_ECG_HD5 = auto()
  RESTING_ECG_SVG = auto()
  EXERCISE_ECG_HD5 = auto()
  BRAIN_MRI_DICOM = auto()
  CARDIAC_MRI_DICOM = auto()
  MRI_HD5 = auto()

# Three-level dictionary of data locations.
FOLDERS = collections.defaultdict(lambda: collections.defaultdict(dict))

# Terra does not currently support attaching persistent disks, so instead we
# read the subset of the data available from Cloud Storage. Alternately, at a
# later date we might prefer to make use of the workspace-associated bucket for
# fake data so that it could be make publicly available.
FOLDERS[Runtime.TERRA] = {
    Dataset.UKB: {
        DataType.EXERCISE_ECG_HD5: 'gs://ml4cvd/deflaux/ukbb_tensors/',
        DataType.RESTING_ECG_HD5: 'gs://ml4cvd/deflaux/ukbb_tensors/',
        DataType.RESTING_ECG_SVG: 'gs://ml4cvd/ecg_views_11_04_2019_svg/',
        DataType.MRI_HD5: 'gs://ml4cvd/deflaux/ukbb_tensors/',
        DataType.BRAIN_MRI_DICOM: 'gs://bulkml4cvd/brainmri/t1_structural_07_26_2019/zipped_t1_dicoms/',
        DataType.CARDIAC_MRI_DICOM: 'gs://ml4cvd/data/mris/cardiac/',
    },
    Dataset.FAKE: {
        # If fake data is not available, put in the path to the real data.
        # Dependent code must gracefully handle 'not found' conditions.
        DataType.EXERCISE_ECG_HD5: 'gs://ml4cvd/projects/fake_hd5s/',
        DataType.RESTING_ECG_HD5: 'gs://ml4cvd/projects/fake_hd5s/',
        DataType.RESTING_ECG_SVG: 'gs://ml4cvd/ecg_views_fake/',
        DataType.MRI_HD5: 'gs://ml4cvd/projects/fake_hd5s/',
        DataType.BRAIN_MRI_DICOM: 'gs://bulkml4cvd/brainmri/t1_structural_07_26_2019/zipped_t1_dicoms/',
        DataType.CARDIAC_MRI_DICOM: 'gs://ml4cvd/data/mris/cardiac/',
    },
}

# The full data is available on an attached persistent disk, read from there.
# [Reading from Cloud Storage would also work in this environment, but its not
# the preferred source in this runtime environment.]
FOLDERS[Runtime.ML4H_VM] = {
    Dataset.UKB: {
        DataType.EXERCISE_ECG_HD5: '/mnt/ml4cvd/deflaux/ukbb_tensors/',
        DataType.RESTING_ECG_HD5: '/mnt/ml4cvd/deflaux/ukbb_tensors/',
        DataType.RESTING_ECG_SVG: '/mnt/ml4cvd/ecg_views_11_04_2019_svg/',
        DataType.MRI_HD5: '/mnt/ml4cvd/deflaux/ukbb_tensors/',
        DataType.BRAIN_MRI_DICOM: 'gs://bulkml4cvd/brainmri/t1_structural_07_26_2019/zipped_t1_dicoms/',
        DataType.CARDIAC_MRI_DICOM: '/mnt/ml4cvd/data/mris/cardiac/',
    },
    Dataset.FAKE: {
        # If fake data is not available, put in the path to the real data.
        # Dependent code must gracefully handle 'not found' conditions.
        DataType.EXERCISE_ECG_HD5: '/mnt/ml4cvd/projects/fake_hd5s/',
        DataType.RESTING_ECG_HD5: '/mnt/ml4cvd/projects/fake_hd5s/',
        DataType.RESTING_ECG_SVG: '/mnt/ml4cvd/ecg_views_fake/',
        DataType.MRI_HD5: '/mnt/ml4cvd/projects/fake_hd5s/',
        DataType.BRAIN_MRI_DICOM: 'gs://bulkml4cvd/brainmri/t1_structural_07_26_2019/zipped_t1_dicoms/',
        DataType.CARDIAC_MRI_DICOM: '/ml4cvd/data/mris/cardiac/',
    },
}

# Configure other environements such as AI Platform Notebooks to
# also read from Cloud Storage.
FOLDERS[Runtime.DEFAULT] = FOLDERS[Runtime.TERRA]


def determine_runtime():
  """Infer the runtime environment."""
  if 'GOOGLE_PROJECT' in os.environ:
    return Runtime.TERRA
  elif os.path.isdir('/mnt/disks/'):
    return Runtime.ML4H_VM
  else:
    return Runtime.DEFAULT


def determine_dataset(sample_id):
  """Given a sample id, determine to which dataset it belongs."""
  if 'fake' in str(sample_id):
    return Dataset.FAKE
  return Dataset.UKB


def get_resting_ecg_hd5_folder(sample_id):
  return FOLDERS[determine_runtime()][
      determine_dataset(sample_id)
  ][DataType.RESTING_ECG_HD5]


def get_resting_ecg_svg_folder(sample_id):
  return FOLDERS[determine_runtime()][
      determine_dataset(sample_id)
  ][DataType.RESTING_ECG_SVG]


def get_exercise_ecg_hd5_folder(sample_id):
  return FOLDERS[determine_runtime()][
      determine_dataset(sample_id)
  ][DataType.EXERCISE_ECG_HD5]


def get_brain_mri_folder(sample_id):
  return FOLDERS[determine_runtime()][
      determine_dataset(sample_id)
  ][DataType.BRAIN_MRI_DICOM]


def get_cardiac_mri_folder(sample_id):
  return FOLDERS[determine_runtime()][
      determine_dataset(sample_id)
  ][DataType.CARDIAC_MRI_DICOM]


def get_mri_hd5_folder(sample_id):
  return FOLDERS[determine_runtime()][
      determine_dataset(sample_id)
  ][DataType.MRI_HD5]


def get_mri_folders(sample_id):
  return [
      FOLDERS[determine_runtime()][
          determine_dataset(sample_id)
      ][DataType.BRAIN_MRI_DICOM],
      FOLDERS[determine_runtime()][
          determine_dataset(sample_id)
      ][DataType.CARDIAC_MRI_DICOM],
  ]
