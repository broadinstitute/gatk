import logging
import tempfile

import apache_beam as beam
import h5py
from apache_beam import Pipeline
from apache_beam.options.pipeline_options import PipelineOptions
from tensorize.defines import TENSOR_EXT, GCS_BUCKET

from .utils import count_ones, _dataset_name_from_meaning, _to_float_or_false, get_gcs_bucket, get_field_type


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
    # TODO: Don't hardcode LIMIT, etc. in the query
    # limit = 450
    query = """
        SELECT d.field, d.valuetype, p_sub_f_c.meaning, p_sub_f_c.sample_id, p_sub_f_c.fieldid, p_sub_f_c.instance, p_sub_f_c.array_idx, p_sub_f_c.value FROM
        (
            SELECT c.meaning, p_sub.sample_id, p_sub.fieldid, p_sub.instance, p_sub.array_idx, p_sub.value FROM
                (
                    SELECT p.sample_id, p.fieldid, p.instance, p.array_idx, p.value, p.coding_file_id FROM `ukbb_dev.phenotype` p
                    INNER JOIN (SELECT * FROM `shared_data.tensorization_fieldids`) f
                    ON p.fieldid = f.fieldid
                ) AS p_sub
            JOIN `ukbb_dev.coding` c
            ON c.coding=p_sub.value AND c.coding_file_id=p_sub.coding_file_id
        ) AS p_sub_f_c
        JOIN `ukbb_dev.dictionary` d
        ON p_sub_f_c.fieldid = d.fieldid"""
    
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
                    value_type = row['valuetype']
                    meaning = row['meaning']

                    field_type = get_field_type(value_type)

                    dataset_name = None
                    if field_type == 'categorical':
                        dataset_name = _dataset_name_from_meaning('categorical', [field, meaning, str(instance), str(array_idx)])
                    elif field_type == 'continuous':
                        dataset_name = _dataset_name_from_meaning('continuous', [str(field_id), field, str(instance), str(array_idx)])

                    if dataset_name is not None:
                        float_value = _to_float_or_false(value)
                        if float_value is not False:
                            hd5.create_dataset(dataset_name, data=[float_value])
                    else:
                        logging.warning("Cannot cast to float from '{}' for field id '{}' and sample id '{}'".format(value, field_id, sample_id))
            gcs_blob.upload_from_filename(tensor_path)

    except:
        logging.exception("Problem with processing sample id '{}'".format(sample_id))
