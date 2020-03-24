"""Runtime-specific data location constants and associated helper methods.

TODO(everyone):
* Add more datasets and data types as needed.
* Update logic to determine runtime as needed.
* Update Cloud Storage and persistent disk paths as needed.
"""

import collections
import os
from enum import Enum


class Runtime(Enum):
  TERRA = 1
  ML4CVD_VM = 2
  DEFAULT = 3  # Such as an AI Platform Notebooks VM.
  # TODO(sam): PARTNERS = 4


class Dataset(Enum):
  FAKE = 1
  UKB = 2


class DataType(Enum):
  RESTING_ECG_HD5 = 1
  RESTING_ECG_SVG = 2
  EXERCISE_ECG_HD5 = 3
  BRAIN_MRI = 4
  CARDIAC_MRI = 5

# Three-level dictionary of data locations.
FOLDERS = collections.defaultdict(lambda: collections.defaultdict(dict))

# Terra does not currently support attaching persistent disks, so instead we
# read the subset of the data available from Cloud Storage. Alternately, at a
# later date we might prefer to make use of the workspace-associated bucket for
# fake data so that it could be make publicly available.
FOLDERS[Runtime.TERRA] = {
    Dataset.UKB: {
        DataType.EXERCISE_ECG_HD5: 'gs://ml4cvd/ecg_bike/',
        DataType.RESTING_ECG_HD5: 'gs://ml4cvd/rest-ecg-hd5s/2019-11-19/',
        DataType.RESTING_ECG_SVG: 'gs://ml4cvd/ecg_views_11_04_2019_svg/',
        DataType.BRAIN_MRI: 'gs://bulkml4cvd/brainmri/t1_structural_07_26_2019/zipped_t1_dicoms/',
        DataType.CARDIAC_MRI: 'gs://ml4cvd/data/mris/cardiac/',
    },
    Dataset.FAKE: {
        # If fake data is not available, put in the path to the real data.
        # Dependent code must gracefully handle 'not found' conditions.
        DataType.EXERCISE_ECG_HD5: 'gs://ml4cvd/ecg_bike/',
        DataType.RESTING_ECG_HD5: 'gs://ml4cvd/projects/fake_ecgs/',
        DataType.RESTING_ECG_SVG: 'gs://ml4cvd/ecg_views_fake/',
        DataType.BRAIN_MRI: 'gs://bulkml4cvd/brainmri/t1_structural_07_26_2019/zipped_t1_dicoms/',
        DataType.CARDIAC_MRI: 'gs://ml4cvd/projects/fake_mris/',
    },
}

# The full data is available on an attached persistent disk, read from there.
# [Reading from Cloud Storage would also work in this environment, but its not
# the preferred source in this runtime environment.]
FOLDERS[Runtime.ML4CVD_VM] = {
    Dataset.UKB: {
        DataType.EXERCISE_ECG_HD5: '/mnt/disks/ecg-bike-tensors/2019-10-10/',
        DataType.RESTING_ECG_HD5: '/mnt/disks/data/ml4cvd/rest-ecg-hd5s/2019-11-19/',
        DataType.RESTING_ECG_SVG: '/mnt/disks/data/ml4cvd/ecg_views_11_04_2019_svg/',
        DataType.BRAIN_MRI: '/mnt/disks/data/bulkml4cvd/brainmri/t1_structural_07_26_2019/zipped_t1_dicoms/',
        DataType.CARDIAC_MRI: '/mnt/disks/data/ml4cvd/data/mris/cardiac/',
    },
    Dataset.FAKE: {
        # If fake data is not available, put in the path to the real data.
        # Dependent code must gracefully handle 'not found' conditions.
        DataType.EXERCISE_ECG_HD5: '/mnt/disks/ecg-bike-tensors/2019-10-10/',
        DataType.RESTING_ECG_HD5: '/mnt/disks/data/ml4cvd/projects/fake_ecgs/',
        DataType.RESTING_ECG_SVG: '/mnt/disks/data/ml4cvd/ecg_views_fake/',
        DataType.BRAIN_MRI: '/mnt/disks/data/bulkml4cvd/brainmri/t1_structural_07_26_2019/zipped_t1_dicoms/',
        DataType.CARDIAC_MRI: '/mnt/disks/data/ml4cvd/projects/fake_mris/',
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
    return Runtime.ML4CVD_VM
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
  ][DataType.BRAIN_MRI]


def get_cardiac_mri_folder(sample_id):
  return FOLDERS[determine_runtime()][
      determine_dataset(sample_id)
  ][DataType.CARDIAC_MRI]


def get_mri_folders(sample_id):
  return [
      FOLDERS[determine_runtime()][
          determine_dataset(sample_id)
      ][DataType.BRAIN_MRI],
      FOLDERS[determine_runtime()][
          determine_dataset(sample_id)
      ][DataType.CARDIAC_MRI],
  ]
