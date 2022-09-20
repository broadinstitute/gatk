import unittest

from hail_gvs_import import generate_avro_args
from unittest.mock import MagicMock


class ShamBlob:
    def __init__(self, name):
        self.name = name


class ShamBucket:
    def __init__(self, name):
        self.name = name

    def list_blobs(self):
        pass


class TestGenerateAvroDict(unittest.TestCase):

    def test_generate_avro_args(self):

        bucket_name = 'fc-workspace'
        object_prefix = f'gs://{bucket_name}/submission/workflow/workflow-id/call-Foo/avro'

        # non superpartitioned
        blobs = [
            ShamBlob(f'{object_prefix}/sample_mapping/sample_mapping_000000000000.avro')
        ]
        bucket = ShamBucket(bucket_name)
        bucket.list_blobs = MagicMock(return_value=blobs)

        actual = generate_avro_args(
            bucket=bucket,
            object_prefix=object_prefix,
            key='sample_mapping'
        )
        bucket.list_blobs.assert_called_with(prefix=f'{object_prefix}/sample_mapping/')
        expected = [b.name for b in blobs]
        self.assertEqual(actual, expected)

        # superpartitioned
        blobs = [
            ShamBlob(f'{object_prefix}/vets/vet_001/vet_001_000000000000.avro'),
            ShamBlob(f'{object_prefix}/vets/vet_001/vet_001_000000000001.avro'),
            ShamBlob(f'{object_prefix}/vets/vet_001/vet_001_000000000002.avro'),
            ShamBlob(f'{object_prefix}/vets/vet_002/vet_002_000000000000.avro'),
            ShamBlob(f'{object_prefix}/vets/vet_002/vet_002_000000000001.avro'),
            ShamBlob(f'{object_prefix}/vets/vet_002/vet_002_000000000002.avro'),
            ShamBlob(f'{object_prefix}/vets/vet_003/vet_003_000000000000.avro'),
            ShamBlob(f'{object_prefix}/vets/vet_003/vet_003_000000000001.avro'),
            ShamBlob(f'{object_prefix}/vets/vet_003/vet_003_000000000002.avro'),
        ]

        bucket.list_blobs = MagicMock(return_value=blobs)
        actual = generate_avro_args(
            bucket=bucket,
            object_prefix=object_prefix,
            key='vets'
        )
        bucket.list_blobs.assert_called_with(prefix=f'{object_prefix}/vets/')

        def vet_avros_for_prefix(prefix):
            return [f'{object_prefix}/vets/{prefix}/{prefix}_00000000000{i}.avro' for i in range(3)]

        expected = [
            vet_avros_for_prefix(f'vet_00{i+1}') for i in range(3)
        ]
        self.assertEqual(actual, expected)
