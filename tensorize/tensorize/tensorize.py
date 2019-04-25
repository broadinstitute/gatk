import logging
import os
import tempfile
import time

import apache_beam as beam
import h5py
from apache_beam.options.pipeline_options import PipelineOptions

from .utils import count_ones, process_entries


def run(pipeline_options: PipelineOptions, output_file: str):
    p = beam.Pipeline(options=pipeline_options)

    bigquery_source = beam.io.BigQuerySource(
        query='select * from `lubitz.coding`',
        use_standard_sql=True
    )

    # Query table in BQ
    table_data = (
        p
        | 'QueryTable' >> beam.io.Read(bigquery_source)
        # Each row is a dictionary where the keys are the BigQuery columns

        | 'CreateKey' >> beam.Map(lambda elem: (elem['coding_file_id'], 1))
        # group by key

        | 'GroupByKey' >> beam.GroupByKey()
        # count entries per key

        | 'Count' >> beam.Map(count_ones)
    )

    # Should be replaced by the schema.json
    # table_schema = bigquery.TableSchema()
    # coding_file_id_field = bigquery.TableFieldSchema()
    # coding_file_id_field.name = 'coding_file_id'
    # coding_file_id_field.type = 'INTEGER'
    # table_schema.fields.append(coding_file_id_field)
    #
    # coding_field = bigquery.TableFieldSchema()
    # coding_field.name = 'coding'
    # coding_field.type = 'STRING'
    # table_schema.fields.append(coding_field)
    #
    # meaning_field = bigquery.TableFieldSchema()
    # meaning_field.name = 'meaning'
    # meaning_field.type = 'STRING'
    # table_schema.fields.append(meaning_field)

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

    # Write to file
    table_data | 'WriteToFile' >> beam.io.WriteToText(output_file)

    result = p.run()
    result.wait_until_finish()


def run2(pipeline_options: PipelineOptions, output_file: str):
    limit = 500

    p = beam.Pipeline(options=pipeline_options)

    bigquery_source = beam.io.BigQuerySource(
        # query='select * from `ukbb7089_r10data.phenotype` limit %s' % limit,
        query='select * from `ukbb7089_r10data.phenotype`',
        use_standard_sql=True
    )

    # Query table in BQ
    table_data = (
            p
            | 'QueryTable' >> beam.io.Read(bigquery_source)

            # Each row is a dictionary where the keys are the BigQuery columns
            | 'CreateKey' >> beam.Map(lambda row: (row['sample_id'], row))

            # Group by key
            | 'GroupByKey' >> beam.GroupByKey()

            # Create dict of sample_id -> list(row dicts)
            | 'ProcessSamples' >> beam.Map(process_entries)
    )

    # Write to file
    table_data | 'WriteToFile' >> beam.io.WriteToText(output_file)

    result = p.run()
    result.wait_until_finish()


def write_tensor_from_sql(sampleid_to_rows):
    sample_id, rows = sampleid_to_rows

    try:
        with tempfile.TemporaryDirectory() as temp_dir:
            tensor_file = '{}.hd5'.format(sample_id)
            tensor_path = '{}/{}'.format(temp_dir, tensor_file)
            logging.info("Writing tensor {} ...".format(tensor_file))
            with h5py.File(tensor_path, 'w') as hd5:
                for row in rows:
                    field_id = row['FieldID']
                    field = row['Field']
                    instance = row['instance']
                    array_idx = row['array_idx']
                    value = row['value']
                    meaning = row['meaning']

                    dataset_name = _dataset_name_from_meaning('categorical', [field, meaning, instance, array_idx])
                    float_category = _to_float_or_false(value)
                    if float_category is not False:
                        hd5.create_dataset(dataset_name, data=[float_category])
                    else:
                        logging.warning('Cannot cast float from: {} categorical field: {} means: {} sample id: {}'.format(
                            value, field_id, meaning, sample_id))
    except:
        logging.exception("problem with processing sample id '{}'".format(sample_id))


def _dataset_name_from_meaning(group: str, fields: List[str]) -> str:
    clean_fields = []
    for f in fields:
        clean_fields.append(''.join(e for e in f if e.isalnum() or e == ' '))
    joined = JOIN_CHAR.join(clean_fields).replace('  ', CONCAT_CHAR).replace(' ', CONCAT_CHAR)
    return group + HD5_GROUP_CHAR + joined


def _to_float_or_false(s):
    try:
        return float(s)
    except ValueError:
        return False
