import logging

GCS_BUCKET = 'ml4cvd'
GCS_BLOB_PATH = 'dataflow_experiment/output/tensors'

REQUIREMENTS_FILE = '/Users/kyuksel/github/ml/env/requirements_ml4cvd_py35.txt'
SETUP_FILE = '/Users/kyuksel/github/ml/tensorize/setup.py'

RUN_NAME = 'ky-test'

FS_OUTPUT_FILE = '/Users/kyuksel/ml4cvd/dataflow/output/%s' % RUN_NAME
GCS_OUTPUT_FILE = 'gs://ml4cvd/dataflow_experiment/output/%s' % RUN_NAME

DIRECT_RUNNER = 'DirectRunner'
DATAFLOW_RUNNER = 'DataflowRunner'

LOG_LEVEL = logging.DEBUG

TENSOR_EXT = 'hd5'

JOIN_CHAR = '_'
CONCAT_CHAR = '-'
HD5_GROUP_CHAR = '/'
