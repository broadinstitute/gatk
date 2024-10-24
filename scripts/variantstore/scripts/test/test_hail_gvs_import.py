import unittest

from hail_gvs_import import gcs_generate_avro_args
from unittest.mock import MagicMock


class ShamBlob:
    def __init__(self, name):
        self.name = name


class ShamBucket:
    def __init__(self, name):
        self.name = name

    def list_blobs(self):
        pass


class TestGenerateAvroArguments(unittest.TestCase):

    def test_generate_avro_arguments(self):

        bucket_name = 'fc-workspace'
        blob_prefix = f'submission/workflow/workflow-id/call-Foo/avro'

        # non superpartitioned
        blobs = [
            ShamBlob(f'{blob_prefix}/sample_mapping/sample_mapping_000000000000.avro')
        ]
        bucket = ShamBucket(bucket_name)
        bucket.list_blobs = MagicMock(return_value=blobs)

        actual = gcs_generate_avro_args(
            bucket=bucket,
            blob_prefix=blob_prefix,
            key='sample_mapping'
        )
        bucket.list_blobs.assert_called_with(prefix=f'{blob_prefix}/sample_mapping/')
        expected = [f'gs://{bucket.name}/{b.name}' for b in blobs]
        self.assertEqual(actual, expected)

        # superpartitioned
        blobs = [
            ShamBlob(f'{blob_prefix}/vets/vet_001/vet_001_000000000000.avro'),
            ShamBlob(f'{blob_prefix}/vets/vet_001/vet_001_000000000001.avro'),
            ShamBlob(f'{blob_prefix}/vets/vet_001/vet_001_000000000002.avro'),
            ShamBlob(f'{blob_prefix}/vets/vet_002/vet_002_000000000000.avro'),
            ShamBlob(f'{blob_prefix}/vets/vet_002/vet_002_000000000001.avro'),
            ShamBlob(f'{blob_prefix}/vets/vet_002/vet_002_000000000002.avro'),
            ShamBlob(f'{blob_prefix}/vets/vet_003/vet_003_000000000000.avro'),
            ShamBlob(f'{blob_prefix}/vets/vet_003/vet_003_000000000001.avro'),
            ShamBlob(f'{blob_prefix}/vets/vet_003/vet_003_000000000002.avro'),
        ]

        bucket.list_blobs = MagicMock(return_value=blobs)
        actual = gcs_generate_avro_args(
            bucket=bucket,
            blob_prefix=blob_prefix,
            key='vets'
        )
        bucket.list_blobs.assert_called_with(prefix=f'{blob_prefix}/vets/')

        def vet_avros_for_prefix(prefix):
            return [f'gs://{bucket.name}/{blob_prefix}/vets/{prefix}/{prefix}_00000000000{i}.avro' for i in range(3)]

        expected = [
            vet_avros_for_prefix(f'vet_00{i+1}') for i in range(3)
        ]
        self.assertEqual(actual, expected)
