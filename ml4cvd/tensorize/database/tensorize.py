import logging
import tempfile

import h5py
import apache_beam as beam
from apache_beam import Pipeline
from google.cloud import storage

from ml4cvd.defines import TENSOR_EXT, GCS_BUCKET
from ml4cvd.tensorize.utils import dataset_name_from_meaning, to_float_or_false


def tensorize_categorical_continuous_fields(pipeline: Pipeline,
                                            output_path: str,
                                            bigquery_dataset: str,
                                            tensor_type: str):
    categorical_query = f"""
        SELECT c.meaning, p_f_d.sample_id, p_f_d.fieldid, p_f_d.field, p_f_d.instance, p_f_d.array_idx, p_f_d.value
        FROM
        (
            SELECT f_d.field, p.sample_id, p.fieldid, p.instance, p.array_idx, p.value, p.coding_file_id
            FROM `{bigquery_dataset}.phenotype` p
            INNER JOIN 
            (
              SELECT d.fieldid, d.field
              FROM `shared_data.tensorization_fieldids` f
              INNER JOIN `{bigquery_dataset}.dictionary` d
              ON f.fieldid = d.fieldid
              WHERE d.valuetype IN ('Categorical single', 'Categorical multiple')
            ) AS f_d
            ON f_d.fieldid = p.fieldid
        ) AS p_f_d
        INNER JOIN `{bigquery_dataset}.coding` c
        ON c.coding=p_f_d.value AND c.coding_file_id=p_f_d.coding_file_id
    """

    continuous_query = f"""
        SELECT f_d.field, p.sample_id, p.fieldid, p.instance, p.array_idx, p.value
        FROM `{bigquery_dataset}.phenotype` p
        INNER JOIN 
        (
          SELECT d.fieldid, d.field
          FROM `shared_data.tensorization_fieldids` f
          INNER JOIN `{bigquery_dataset}.dictionary` d
          ON f.fieldid = d.fieldid
          WHERE d.valuetype IN ('Continuous', 'Integer')
        ) AS f_d
        ON f_d.fieldid = p.fieldid
    """

    query = None
    if tensor_type == 'categorical':
        query = categorical_query
    elif tensor_type == 'continuous':
        query = continuous_query
    else:
        raise ValueError("Can tensorize only categorical or continuous fields, got ", tensor_type)

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
            | 'CreateHd5sAndUploadToGCS' >> beam.Map(write_tensor_from_sql, output_path, tensor_type)
    )

    result = pipeline.run()
    result.wait_until_finish()


# We are instantiating the GCS client in global scope because passing it explicitly into a method used by beam.Map()
# gives a 'client not picklable` error. We're also enclosing it in a 'try' block so when this module is imported and if
# the client instantiation fails (when training models in Docker for instance), the whole run doesn't error out.
try:
    gcs_client = storage.Client()
    output_bucket = gcs_client.get_bucket(GCS_BUCKET)
except OSError:
    output_bucket ='nope'
    logging.warning(f"no GCS storage client")


def write_tensor_from_sql(sampleid_to_rows, output_path, tensor_type):
    # GroupByKey output is not a list of dicts, as expected, but something called an 'UnwindowedValue'.
    # We convert it to list first to arrive at (key, [{dict1}, {dict2}, {dict3}...])
    (sample_id, values) = sampleid_to_rows
    rows = list(values)

    try:
        with tempfile.TemporaryDirectory() as temp_dir:
            tensor_file = f"{sample_id}{TENSOR_EXT}"
            tensor_path = f"{temp_dir}/{tensor_file}"
            gcs_blob = output_bucket.blob(f"{output_path}/{tensor_file}")
            logging.info(f"Writing tensor {tensor_file} to {gcs_blob.public_url} ...")
            with h5py.File(tensor_path, 'w') as hd5:
                for row in rows:
                    field_id = row['fieldid']
                    field = row['field']
                    instance = row['instance']
                    array_idx = row['array_idx']
                    value = row['value']

                    hd5_dataset_name = None
                    if tensor_type == 'categorical':
                        meaning = row['meaning']
                        hd5_dataset_name = dataset_name_from_meaning('categorical', [field, meaning, str(instance), str(array_idx)])
                    elif tensor_type == 'continuous':
                        hd5_dataset_name = dataset_name_from_meaning('continuous', [str(field_id), field, str(instance), str(array_idx)])
                    else:
                        continue

                    float_value = to_float_or_false(value)
                    if float_value is not False:
                        hd5.create_dataset(hd5_dataset_name, data=[float_value])
                    else:
                        logging.warning(f"Cannot cast to float from '{value}' for field id '{field_id}' and sample id '{sample_id}'")
            gcs_blob.upload_from_filename(tensor_path)

    except:
        logging.exception(f"Problem with processing sample id '{sample_id}'")