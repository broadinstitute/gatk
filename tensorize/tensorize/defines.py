import argparse
import datetime


def parse_args():
    now_string = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')

    parser = argparse.ArgumentParser()

    parser.add_argument('--id', default=f"run_{now_string}",
                        help='User-defined identifier for this pipeline run. '
                             'Per Google: the name must consist of only the characters [-a-z0-9], '
                             'starting with a letter and ending with a letter or number.')
    parser.add_argument('--tensor_type', default="categorical",
                        help='Type of data to be tensorized',
                        choices=['categorical', 'continuous'])
    parser.add_argument('--bigquery_dataset', default='ukbb_dev',
                        help='BigQuery dataset where the data will be drawn from')
    parser.add_argument('--beam_runner', default='DirectRunner',
                        help='Apache Beam runner that will execute the pipeline',
                        choices=['DirectRunner', 'DataflowRunner'])
    parser.add_argument('--repo_root',
                        help='Root directory of the cloned ml repo')
    parser.add_argument('--gcp_project', default='broad-ml4cvd',
                        help='Name of the Google Cloud Platform project')
    parser.add_argument('--gcp_region', default='us-east1',
                        help='Google Cloud Platform region')
    parser.add_argument('--gcs_bucket', default='ml4cvd',
                        help='Name of the Google Cloud Storage bucket where tensors will be written to')
    parser.add_argument('--gcs_output_path',
                        help='gs:// folder path excluding the bucket name where tensors will be written to '
                             'e.g. specifying /path/to/folder will write to gs://<gcs_bucket>/path/to/folder')
    parser.add_argument("--logging_level", default='INFO', help="Logging level",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])

    return parser.parse_args()


args = parse_args()

# ------------- COMMAND LINE ARGS --------------
RUN_NAME = args.id
TENSOR_TYPE = args.tensor_type
BIGQUERY_DATASET = args.bigquery_dataset
RUNNER = args.beam_runner
GCS_BUCKET = args.gcs_bucket
GCS_OUTPUT_PATH = args.gcs_output_path
LOG_LEVEL = args.logging_level
PROJECT = args.gcp_project
REGION = args.gcp_region
# -----------------------------------------------

# ------------- DERIVED ARGS --------------------
# Setup file and requirements derived from the Python environment via 'pip freeze > requirements.txt'
# They are used by Dataflow to pip install the necessary packages
REQUIREMENTS_FILE = f"{args.repo_root}/env/requirements_ml4cvd_dataflow.txt"
SETUP_FILE = f"{args.repo_root}/tensorize/setup.py"

# Needed to be able to submit a pipeline to Dataflow
DATAFLOW_STAGING_LOCATION = f"gs://{args.gcs_bucket}/dataflow/staging"
DATAFLOW_TEMP_LOCATION = f"gs://{args.gcs_bucket}/dataflow/temp"
# ------------------------------------------------

# ------------- CONSTANTS ------------------------
TENSOR_EXT = 'hd5'
JOIN_CHAR = '_'
CONCAT_CHAR = '-'
HD5_GROUP_CHAR = '/'
# ------------------------------------------------
