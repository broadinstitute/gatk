import logging

import apache_beam as beam
from apache_beam.options.pipeline_options import PipelineOptions, GoogleCloudOptions, StandardOptions

from tensorize import tensorize
from tensorize.defines import REQUIREMENTS_FILE, SETUP_FILE, LOG_LEVEL, GCS_BLOB_PATH, RUNNER


def get_pipeline_options(argv, runner: str) -> PipelineOptions:
    pipeline_opts = PipelineOptions(flags=argv)

    google_cloud_options = pipeline_opts.view_as(GoogleCloudOptions)
    google_cloud_options.region = 'us-east1'
    google_cloud_options.project = 'broad-ml4cvd'
    google_cloud_options.job_name = 'ky-dataflow-test'
    google_cloud_options.staging_location = 'gs://ml4cvd/dataflow_experiment/staging_location'
    google_cloud_options.temp_location = 'gs://ml4cvd/dataflow_experiment/temp_location'

    pipeline_opts.view_as(StandardOptions).runner = runner

    return pipeline_opts


if __name__ == "__main__":
    logging.getLogger().setLevel(LOG_LEVEL)

    packaging_args = [
        '--requirements_file={}'.format(REQUIREMENTS_FILE),
        '--setup_file={}'.format(SETUP_FILE)
    ]

    pipeline_options = get_pipeline_options(packaging_args, RUNNER)
    pipeline = beam.Pipeline(options=pipeline_options)

    tensorize.tensorize_categorical_continuous_fields(pipeline, GCS_BLOB_PATH)
