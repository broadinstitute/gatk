import apache_beam as beam
from apache_beam.options.pipeline_options import PipelineOptions
from apache_beam.options.pipeline_options import GoogleCloudOptions
from apache_beam.options.pipeline_options import StandardOptions

from apache_beam.io.gcp.internal.clients import bigquery

from .defines import DATAFLOW_RUNNER, GCS_OUTPUT_FILE, DIRECT_RUNNER, FS_OUTPUT_FILE
from .utils import count_ones


RUNNER = DATAFLOW_RUNNER
OUTPUT_FILE = GCS_OUTPUT_FILE


def run(argv=None):

    pipeline_options = PipelineOptions(flags=argv)

    google_cloud_options = pipeline_options.view_as(GoogleCloudOptions)
    google_cloud_options.region = 'us-east1'
    google_cloud_options.project = 'broad-ml4cvd'
    google_cloud_options.job_name = 'ky-dataflow-test'
    google_cloud_options.staging_location = 'gs://ml4cvd/dataflow_experiment/staging_location'
    google_cloud_options.temp_location = 'gs://ml4cvd/dataflow_experiment/temp_location'

    pipeline_options.view_as(StandardOptions).runner = RUNNER

    p = beam.Pipeline(options=pipeline_options)

    # Query table in BQ
    table_data = (
            p
            | 'QueryTable' >> beam.io.Read(
        beam.io.BigQuerySource(
            query='select * from `lubitz.coding`',
            use_standard_sql=True
        )
    ) |
            # Each row is a dictionary where the keys are the BigQuery columns
            'create key' >> beam.Map(lambda elem: (elem['coding_file_id'], 1))
            |
            # group by key
            'group by key' >> beam.GroupByKey()
            |
            # count entries per key
            'count' >> beam.Map(count_ones)
    )

    #should be replaced by the schema.json

    table_schema = bigquery.TableSchema()
    coding_file_id_field = bigquery.TableFieldSchema()
    coding_file_id_field.name = 'coding_file_id'
    coding_file_id_field.type = 'INTEGER'
    table_schema.fields.append(coding_file_id_field)

    coding_field = bigquery.TableFieldSchema()
    coding_field.name = 'coding'
    coding_field.type = 'STRING'
    table_schema.fields.append(coding_field)

    meaning_field = bigquery.TableFieldSchema()
    meaning_field.name = 'meaning'
    meaning_field.type = 'STRING'
    table_schema.fields.append(meaning_field)

    # write_table_spec = bigquery.TableReference(
    #     projectId='broad-ml4cvd',
    #     datasetId='pb_working',
    #     tableId=RUN_NAME
    # )
    #
    # write transformed to new table
    # table_data | 'Writing to BigQuery' >> beam.io.WriteToBigQuery(
    #    write_table_spec,
    #    schema=table_schema,
    #    write_disposition=beam.io.BigQueryDisposition.WRITE_TRUNCATE,
    #    create_disposition=beam.io.BigQueryDisposition.CREATE_IF_NEEDED
    # )

    #write to file

    table_data | 'Writing to file %s' % OUTPUT_FILE >> beam.io.WriteToText(OUTPUT_FILE)

    result = p.run()
    result.wait_until_finish()
