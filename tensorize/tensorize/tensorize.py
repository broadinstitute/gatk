import logging
import tempfile

import apache_beam as beam
import h5py
from apache_beam import Pipeline
from tensorize.defines import TENSOR_EXT, GCS_BUCKET, FIELD_TYPE, DATASET

from .utils import _dataset_name_from_meaning, to_float_or_false, get_gcs_bucket


def tensorize_categorical_continuous_fields(pipeline: Pipeline, output_path: str):
    categorical_query = """
        SELECT c.meaning, p_f_d.sample_id, p_f_d.fieldid, p_f_d.field, p_f_d.instance, p_f_d.array_idx, p_f_d.value
        FROM
        (
            SELECT f_d.field, p.sample_id, p.fieldid, p.instance, p.array_idx, p.value, p.coding_file_id
            FROM `{}.phenotype` p
            INNER JOIN 
            (
              SELECT d.fieldid, d.field
              FROM `shared_data.tensorization_fieldids` f
              INNER JOIN `{}.dictionary` d
              ON f.fieldid = d.fieldid
              WHERE d.valuetype IN ('Categorical single', 'Categorical multiple')
            ) AS f_d
            ON f_d.fieldid = p.fieldid
        ) AS p_f_d
        INNER JOIN `{}.coding` c
        ON c.coding=p_f_d.value AND c.coding_file_id=p_f_d.coding_file_id
    """.format(DATASET, DATASET, DATASET)

    continuous_query = """
        SELECT f_d.field, p.sample_id, p.fieldid, p.instance, p.array_idx, p.value
        FROM `{}.phenotype` p
        INNER JOIN 
        (
          SELECT d.fieldid, d.field
          FROM `shared_data.tensorization_fieldids` f
          INNER JOIN `{}.dictionary` d
          ON f.fieldid = d.fieldid
          WHERE d.valuetype IN ('Continuous', 'Integer')
        ) AS f_d
        ON f_d.fieldid = p.fieldid
    """.format(DATASET, DATASET)

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

                    hd5_dataset_name = None
                    if FIELD_TYPE == 'categorical':
                        meaning = row['meaning']
                        hd5_dataset_name = _dataset_name_from_meaning('categorical', [field, meaning, str(instance), str(array_idx)])
                    elif FIELD_TYPE == 'continuous':
                        hd5_dataset_name = _dataset_name_from_meaning('continuous', [str(field_id), field, str(instance), str(array_idx)])
                    else:
                        continue

                    float_value = to_float_or_false(value)
                    if float_value is not False:
                            hd5.create_dataset(hd5_dataset_name, data=[float_value])
                    else:
                        logging.warning("Cannot cast to float from '{}' for field id '{}' and sample id '{}'".format(value, field_id, sample_id))
            gcs_blob.upload_from_filename(tensor_path)

    except:
        logging.exception("Problem with processing sample id '{}'".format(sample_id))
