"""Phenotype Dataset -> Tensors on disk

Currently: * reads from hardcoded ukbb_dev.phenotypes
           * groups by sample_id
           * counts rows
           * dumps to local files
Requires directrunner in python 3.5 environment
"""

import logging
import time
import os

import apache_beam as beam
from apache_beam.options.pipeline_options import PipelineOptions
from apache_beam.options.pipeline_options import GoogleCloudOptions
from apache_beam.options.pipeline_options import StandardOptions

from apache_beam.io.gcp.internal.clients import bigquery

LIMIT=500
DATASET = 'ukbb_dev'
FIELDIDS = 'shared_data.tensorization_fieldids'
OUTPUT_DIRECTORY = os.path.join(os.getcwd(), 'outputs')
if not os.path.exists(OUTPUT_DIRECTORY):
    os.makedirs(OUTPUT_DIRECTORY)
timestr = time.strftime("%Y%m%d_%H%M%S")
OUTPUT_FILE = os.path.join(OUTPUT_DIRECTORY, 'test_%s' % timestr)


def run(argv=None):
    pipeline_options = PipelineOptions(flags=argv)

    google_cloud_options = pipeline_options.view_as(GoogleCloudOptions)
    google_cloud_options.region = 'us-east1'
    google_cloud_options.project = 'broad-ml4cvd'
    google_cloud_options.job_name = 'pb-dataflow-test'
    google_cloud_options.staging_location = 'gs://ml4cvd/projects/pbatra/staging'
    google_cloud_options.temp_location = 'gs://ml4cvd/projects/pbatra/temp'
    #pipeline_options.view_as(StandardOptions).runner = 'DataflowRunner'
    pipeline_options.view_as(StandardOptions).runner = 'DirectRunner'

    p = beam.Pipeline(options=pipeline_options)


    def process_entries(sample_list):
        #values that comes out is not a list of dicts, as expected, but something called a Unwindowed value. Convert to list first, then looks like a list of dicts. So basically, you get a (key, [{dict1}, {dict2}, {dict3}])
        (word, values) = sample_list
        list_dicts = list(values)
        return (word, list_dicts)

    # Query table in BQ
    table_data = (
            p
            | 'QueryTable' >> beam.io.Read(
        beam.io.BigQuerySource(
            query="""select a.* from
            `%s.phenotype` a
            inner join
            `%s` b
            on a.fieldid = b.fieldid
            limit %s"""
            % (DATASET, FIELDIDS, LIMIT),
            use_standard_sql=True
        )
    ) |
            # Each row is a dictionary where the keys are the BigQuery columns
            'create key' >> beam.Map(lambda row: (row['sample_id'], row))
      |
            # group by key
            'group by key' >> beam.GroupByKey()
      |
            # count entries per key
            'process samples' >> beam.Map(process_entries)
    )

    table_data | 'Writing to file %s' % OUTPUT_FILE >> beam.io.WriteToText(OUTPUT_FILE)

    result = p.run()
    result.wait_until_finish()


if __name__ == "__main__":
    logging.getLogger().setLevel(logging.INFO)
    run()
