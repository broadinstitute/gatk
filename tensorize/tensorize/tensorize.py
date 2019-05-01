import logging
import tempfile

import apache_beam as beam
import h5py
from apache_beam import Pipeline
from apache_beam.options.pipeline_options import PipelineOptions
from tensorize.defines import TENSOR_EXT, GCS_BUCKET, FIELD_TYPE, DATASET

from .utils import count_ones, _dataset_name_from_meaning, _to_float_or_false, get_gcs_bucket


def example(pipeline_options: PipelineOptions, output_file: str):
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


def tensorize_categorical_continuous_fields(pipeline: Pipeline, output_path: str):
    categorical_query = """
        SELECT c.meaning, p_f_d.sample_id, p_f_d.fieldid, p_f_d.field, p_f_d.instance, p_f_d.array_idx, p_f_d.value
        FROM
        (
            SELECT f_d.field, p.sample_id, p.fieldid, p.instance, p.array_idx, p.value, p.coding_file_id
            FROM `ukbb_dev.phenotype` p
            INNER JOIN 
            (
              SELECT d.fieldid, d.field
              FROM `shared_data.tensorization_fieldids` f
              INNER JOIN `ukbb_dev.dictionary` d
              ON f.fieldid = d.fieldid
              WHERE d.valuetype IN ('Categorical single', 'Categorical multiple')
            ) AS f_d
            ON f_d.fieldid = p.fieldid
        ) AS p_f_d
        INNER JOIN `ukbb_dev.coding` c
        ON c.coding=p_f_d.value AND c.coding_file_id=p_f_d.coding_file_id"""

    continuous_query = """
        SELECT f_d.field, p.sample_id, p.fieldid, p.instance, p.array_idx, p.value
        FROM `ukbb_dev.phenotype` p
        INNER JOIN 
        (
          SELECT d.fieldid, d.field
          FROM `shared_data.tensorization_fieldids` f
          INNER JOIN `ukbb_dev.dictionary` d
          ON f.fieldid = d.fieldid
          WHERE d.valuetype IN ('Continuous', 'Integer')
        ) AS f_d
        ON f_d.fieldid = p.fieldid"""

    query = None
    if FIELD_TYPE == 'categorical':
        query = categorical_query
    elif FIELD_TYPE == 'continuous':
        query = continuous_query
    else:
        raise ValueError("Can tensorize only categorical or continuous fields, got ", FIELD_TYPE)

    bigquery_source = beam.io.BigQuerySource(query=query, use_standard_sql=True)

    # Query table in BQ
    steps = (
            pipeline
            | 'QueryTables' >> beam.io.Read(bigquery_source)

            # Each row is a dictionary where the keys are the BigQuery columns
            | 'CreateKey' >> beam.Map(lambda row: (row['sample_id'], row))

            # Group by key
            | 'GroupByKey' >> beam.GroupByKey()

            # Format into hd5 files and upload to GCS
            | 'CreateHd5sAndUploadToGCS' >> beam.Map(write_tensor_from_sql, output_path)
    )

    result = pipeline.run()
    result.wait_until_finish()


# Defining this in global scope because passing it explicitly into a method used by beam.Map()
# gives a 'client not picklable` error.
output_bucket = get_gcs_bucket(GCS_BUCKET)


def write_tensor_from_sql(sampleid_to_rows, output_path):
    # GroupByKey output is not a list of dicts, as expected, but something called an 'UnwindowedValue'.
    # We convert it to list first to arrive at (key, [{dict1}, {dict2}, {dict3}...])
    (sample_id, values) = sampleid_to_rows
    rows = list(values)

    try:
        with tempfile.TemporaryDirectory() as temp_dir:
            tensor_file = '{}.{}'.format(sample_id, TENSOR_EXT)
            tensor_path = '{}/{}'.format(temp_dir, tensor_file)
            gcs_blob = output_bucket.blob('{}/{}'.format(output_path, tensor_file))
            logging.info("Writing tensor {} to {} ...".format(tensor_file, gcs_blob.public_url))
            with h5py.File(tensor_path, 'w') as hd5:
                for row in rows:
                    field_id = row['fieldid']
                    field = row['field']
                    instance = row['instance']
                    array_idx = row['array_idx']
                    value = row['value']

                    dataset_name = None
                    if FIELD_TYPE == 'categorical':
                        meaning = row['meaning']
                        dataset_name = _dataset_name_from_meaning('categorical', [field, meaning, str(instance), str(array_idx)])
                    elif FIELD_TYPE == 'continuous':
                        dataset_name = _dataset_name_from_meaning('continuous', [str(field_id), field, str(instance), str(array_idx)])
                    else:
                        continue

                    float_value = _to_float_or_false(value)
                    if float_value is not False:
                            hd5.create_dataset(dataset_name, data=[float_value])
                    else:
                        logging.warning("Cannot cast to float from '{}' for field id '{}' and sample id '{}'".format(value, field_id, sample_id))
            gcs_blob.upload_from_filename(tensor_path)

    except:
        logging.exception("Problem with processing sample id '{}'".format(sample_id))
