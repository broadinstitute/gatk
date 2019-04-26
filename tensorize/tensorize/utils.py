import logging

import h5py
import tempfile
from typing import List

from google.cloud import storage
from tensorize.defines import TENSOR_EXT, JOIN_CHAR, CONCAT_CHAR, HD5_GROUP_CHAR, GCS_BUCKET, GCS_BLOB_PATH


def write_tensor_from_sql(sampleid_to_rows):
    sample_id, rows = sampleid_to_rows

    # TODO: Pass in gcs_client as argument
    gcs_client = storage.Client()
    # TODO: Remove hardcoded gcs project/bucket/blob path
    gcs_bucket = gcs_client.get_bucket(GCS_BUCKET)
    blob_path = GCS_BLOB_PATH

    try:
        with tempfile.TemporaryDirectory() as temp_dir:
            tensor_file = '{}.{}'.format(sample_id, TENSOR_EXT)
            tensor_path = '{}/{}'.format(temp_dir, tensor_file)
            gcs_blob = gcs_bucket.blob('{}/{}'.format(blob_path, tensor_file))
            logging.info("Writing tensor {} to {} ...".format(tensor_file, gcs_blob.public_url))
            with h5py.File(tensor_path, 'w') as hd5:
                for row in rows:
                    field_id = row['fieldid']
                    field = row['field']
                    instance = row['instance']
                    array_idx = row['array_idx']
                    value = row['value']
                    meaning = row['meaning']

                    dataset_name = _dataset_name_from_meaning('categorical', [field, meaning, str(instance), str(array_idx)])
                    float_category = _to_float_or_false(value)
                    if float_category is not False:
                        hd5.create_dataset(dataset_name, data=[float_category])
                    else:
                        logging.warning('Cannot cast float from: {} categorical field: {} means: {} sample id: {}'.format(
                            value, field_id, meaning, sample_id))
            gcs_blob.upload_from_filename(tensor_path)

    except:
        logging.exception("Problem with processing sample id '{}'".format(sample_id))


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


def count_ones(word_ones):
    (word, ones) = word_ones
    return word, sum(ones)


def process_entries(sample_list):
    """ Values that come out is not a list of dicts, as expected, but something called a Unwindowed value.
    Convert to list first, then looks like a list of dicts. So basically, you get a (key, [{dict1}, {dict2}, {dict3}])
    """
    (word, values) = sample_list
    list_dicts = list(values)

    return word, list_dicts
