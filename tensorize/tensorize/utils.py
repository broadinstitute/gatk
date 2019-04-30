from typing import List

from google.cloud import storage
from google.cloud.storage import Bucket
from tensorize.defines import JOIN_CHAR, CONCAT_CHAR, HD5_GROUP_CHAR


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


def get_gcs_bucket(bucket_name: str) -> Bucket:
    gcs_client = storage.Client()
    return gcs_client.get_bucket(bucket_name)


def get_field_type(value_type: str) -> str:
    if value_type == 'Categorical single' or value_type == 'Categorical multiple':
        return 'categorical'
    elif value_type == 'Continuous' or value_type == 'Integer':
        return 'continuous'
    else:
        return value_type
