import argparse
import datetime
import logging

import apache_beam as beam
from apache_beam.options.pipeline_options import PipelineOptions, GoogleCloudOptions, StandardOptions

from ml4h.defines import GCS_BUCKET
from ml4h.tensorize.database import tensorize


def parse_args():
    now_string = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--id', default=f"run_{now_string}",
        help='User-defined identifier for this pipeline run. '
             'Per Google: the name must consist of only the characters [-a-z0-9], '
             'starting with a letter and ending with a letter or number.',
    )
    parser.add_argument(
        '--tensor_type', default="categorical",
        help='Type of data to be tensorized',
        choices=['categorical', 'continuous', 'icd', 'disease', 'death', 'phecode_disease'],
    )
    parser.add_argument(
        '--bigquery_dataset', default='ukbb_dev',
        help='BigQuery dataset where the data will be drawn from',
    )
    parser.add_argument(
        '--beam_runner', default='DirectRunner',
        help='Apache Beam runner that will execute the pipeline',
        choices=['DirectRunner', 'DataflowRunner'],
    )
    parser.add_argument(
        '--repo_root',
        help='Root directory of the cloned ml repo',
    )
    parser.add_argument(
        '--gcp_project', default='broad-ml4cvd',
        help='Name of the Google Cloud Platform project',
    )
    parser.add_argument(
        '--gcp_region', default='us-east1',
        help='Google Cloud Platform region',
    )
    # parser.add_argument('--gcs_bucket', default='ml4h',
    #                     help='Name of the Google Cloud Storage bucket where tensors will be written to')
    parser.add_argument(
        '--gcs_output_path',
        help='gs:// folder path excluding the bucket name where tensors will be written to '
             'e.g. specifying /path/to/folder will write to gs://<gcs_bucket>/path/to/folder',
    )
    parser.add_argument(
        "--logging_level", default='INFO', help="Logging level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    logging.getLogger().setLevel(args.logging_level)

    packaging_args = [
        '--requirements_file={}'.format(f"{args.repo_root}/env/requirements_ml4h_dataflow.txt"),
        '--setup_file={}'.format(f"{args.repo_root}/setup.py"),
    ]

    pipeline_opts = PipelineOptions(flags=packaging_args)
    google_cloud_options = pipeline_opts.view_as(GoogleCloudOptions)
    google_cloud_options.region = args.gcp_region
    google_cloud_options.project = args.gcp_project
    google_cloud_options.job_name = args.id
    google_cloud_options.staging_location = f"gs://{GCS_BUCKET}/dataflow/staging"
    google_cloud_options.temp_location = f"gs://{GCS_BUCKET}/dataflow/temp"
    pipeline_opts.view_as(StandardOptions).runner = args.beam_runner

    pipeline = beam.Pipeline(options=pipeline_opts)

    tensorize.tensorize_sql_fields(
        pipeline,
        args.gcs_output_path,
        args.bigquery_dataset,
        args.tensor_type,
    )
