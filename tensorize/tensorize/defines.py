import logging

# ------------- USER VARIABLES --------------
# ID of the dataflow pipeline run
RUN_NAME = 'winter-has-come'
# Type of the data to be tensorized
FIELD_TYPE = 'categorical'
# BigQuery dataset where the data will be drawn from - set to 'ukbb_dev' for smaller dataset
DATASET = 'ukbb7089_r10data'
# Runner type - set to 'DirectRunner' for local development
RUNNER = 'DataflowRunner'
# Repo root
REPO_ROOT = '/Users/JonSnow/github/ml'
# Logging level - set to DEBUG for more verbosity
LOG_LEVEL = logging.INFO
# --------------------------------------------

# Setup file and requirements derived from the Python environment via 'pip freeze > requirements.txt'
# They are used by Dataflow to pip install the necessary packages
REQUIREMENTS_FILE = f"{REPO_ROOT}/env/requirements_ml4cvd_dataflow.txt"
SETUP_FILE = f"{REPO_ROOT}/tensorize/setup.py"

OUTPUT_FOLDER = f"tensors_{DATASET}_{FIELD_TYPE}"
GCS_BUCKET = 'ml4cvd'
GCS_BLOB_PATH = f"data/tensors/{OUTPUT_FOLDER}"

# Needed to be able to submit a pipeline to Dataflow
DATAFLOW_STAGING_LOCATION = f"gs://{GCS_BUCKET}/dataflow/staging"
DATAFLOW_TEMP_LOCATION = f"gs://{GCS_BUCKET}/dataflow/temp"
PROJECT = 'broad-ml4cvd'
REGION = 'us-east1'

TENSOR_EXT = 'hd5'

JOIN_CHAR = '_'
CONCAT_CHAR = '-'
HD5_GROUP_CHAR = '/'
