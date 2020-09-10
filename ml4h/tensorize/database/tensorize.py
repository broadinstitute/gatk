import logging
import tempfile
from typing import Dict, List

import h5py
import apache_beam as beam
from apache_beam import Pipeline
from google.cloud import storage

from ml4h.defines import TENSOR_EXT, GCS_BUCKET, JOIN_CHAR, CONCAT_CHAR, HD5_GROUP_CHAR, dataset_name_from_meaning


def tensorize_sql_fields(pipeline: Pipeline, output_path: str, sql_dataset: str, tensor_type: str):

    if tensor_type == 'categorical':
        query = _get_categorical_query(sql_dataset)
    elif tensor_type == 'continuous':
        query = _get_continuous_query(sql_dataset)
    elif tensor_type == 'icd':
        query = _get_icd_query(sql_dataset)
    elif tensor_type == 'disease':
        query = _get_disease_query(sql_dataset)
    elif tensor_type == 'phecode_disease':
        query = _get_phecode_query(sql_dataset)
    elif tensor_type == 'death':
        query = _get_death_and_censor_query(sql_dataset)
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
# except OSError:
except:
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
                if tensor_type == 'icd':
                    icds = sorted(list(set([row['value'] for row in rows])))
                    hd5.create_dataset('icd', (1,), data=JOIN_CHAR.join(icds), dtype=h5py.special_dtype(vlen=str))
                elif tensor_type == 'categorical':
                    for row in rows:
                        hd5_dataset_name = dataset_name_from_meaning('categorical', [row['field'], row['meaning'], str(row['instance']), str(row['array_idx'])])
                        _write_float_or_warn(sample_id, row, hd5_dataset_name, hd5)
                elif tensor_type == 'continuous':
                    for row in rows:
                        hd5_dataset_name = dataset_name_from_meaning('continuous', [str(row['fieldid']), row['field'], str(row['instance']), str(row['array_idx'])])
                        _write_float_or_warn(sample_id, row, hd5_dataset_name, hd5)
                elif tensor_type in ['disease', 'phecode_disease']:
                    for row in rows:
                        hd5.create_dataset('categorical' + HD5_GROUP_CHAR + row['disease'].lower(), data=[float(row['has_disease'])])
                        hd5_date = 'dates' + HD5_GROUP_CHAR + row['disease'].lower() + '_date'
                        hd5.create_dataset(hd5_date, (1,), data=str(row['censor_date']), dtype=h5py.special_dtype(vlen=str))
                elif tensor_type == 'death':
                    for row in rows:
                        hd5.create_dataset('categorical' + HD5_GROUP_CHAR + 'death', data=[float(row['has_died'])])
                        d = 'dates' + HD5_GROUP_CHAR
                        hd5.create_dataset(d+'enroll_date', (1,), data=str(row['enroll_date']), dtype=h5py.special_dtype(vlen=str))
                        hd5.create_dataset(d+'death_censor', (1,), data=str(row['death_censor_date']), dtype=h5py.special_dtype(vlen=str))
                        hd5.create_dataset(d+'phenotype_censor', (1,), data=str(row['phenotype_censor_date']), dtype=h5py.special_dtype(vlen=str))
            gcs_blob.upload_from_filename(tensor_path)
    except:
        logging.exception(f"Problem with processing sample id '{sample_id}'")


def _write_float_or_warn(sample_id, row, hd5_dataset_name, hd5):
    try:
        float_value = float(row['value'])
        hd5.create_dataset(hd5_dataset_name, data=[float_value])
    except ValueError:
        logging.warning(f"Cannot cast to float from '{row['value']}' for field id '{row['fieldid']}' and sample id '{sample_id}'")


def _get_categorical_query(dataset):
    return f"""
        SELECT c.meaning, p_f_d.sample_id, p_f_d.fieldid, p_f_d.field, p_f_d.instance, p_f_d.array_idx, p_f_d.value
        FROM
        (
            SELECT f_d.field, p.sample_id, p.fieldid, p.instance, p.array_idx, p.value, p.coding_file_id
            FROM `{dataset}.phenotype` p
            INNER JOIN
            (
              SELECT d.fieldid, d.field
              FROM `shared_data.tensorization_fieldids` f
              INNER JOIN `{dataset}.dictionary` d
              ON f.fieldid = d.fieldid
              WHERE d.valuetype IN ('Categorical single', 'Categorical multiple')
            ) AS f_d
            ON f_d.fieldid = p.fieldid
        ) AS p_f_d
        INNER JOIN `{dataset}.coding` c
        ON c.coding=p_f_d.value AND c.coding_file_id=p_f_d.coding_file_id
    """


def _get_continuous_query(dataset):
    return f"""
        SELECT f_d.field, p.sample_id, p.fieldid, p.instance, p.array_idx, p.value
        FROM `{dataset}.phenotype` p
        INNER JOIN
        (
          SELECT d.fieldid, d.field
          FROM `shared_data.tensorization_fieldids` f
          INNER JOIN `{dataset}.dictionary` d
          ON f.fieldid = d.fieldid
          WHERE d.valuetype IN ('Continuous', 'Integer')
        ) AS f_d
        ON f_d.fieldid = p.fieldid
    """


def _get_icd_query(dataset):
    return f"""
        SELECT sample_id, value FROM `{dataset}.phenotype` WHERE fieldid IN (41202, 41204, 40001, 40002, 40006);
    """


def _get_disease_query(dataset):
    return f"""
        SELECT sample_id, disease, has_disease, censor_date FROM `{dataset}.disease` WHERE has_disease=1;
    """


def _get_death_and_censor_query(dataset):
    return f"""
        SELECT distinct(d.sample_id), d.has_died, d.death_censor_date, c.phenotype_censor_date, d.enroll_date
        FROM `{dataset}.disease` d INNER JOIN `{dataset}.censor` c
         ON c.sample_id = d.sample_id ORDER BY d.sample_id;
    """


def _get_phecode_query(dataset):
    return f"""
        SELECT sample_id, disease, has_disease, censor_date FROM `{dataset}.phecodes_nonzero` WHERE has_disease=1;
    """
