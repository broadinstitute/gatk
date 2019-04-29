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


def process_entries(sample_list):
    """ Values that come out is not a list of dicts, as expected, but something called a Unwindowed value.
    Convert to list first, then looks like a list of dicts. So basically, you get a (key, [{dict1}, {dict2}, {dict3}])
    """
    (word, values) = sample_list
    list_dicts = list(values)
    return word, list_dicts


def get_gcs_bucket(bucket_name: str) -> Bucket:
    gcs_client = storage.Client()
    return gcs_client.get_bucket(bucket_name)
