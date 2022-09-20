import unittest

from hail_gvs_import import generate_avro_args


class TestGenerateAvroDict(unittest.TestCase):
    class MockBucket:
        def __init__(self, name, blobs):
            self.name = name
            self.blobs = blobs

        def list_blobs(self):
            return self.blobs

    def test_generate_avro_args(self):

        # test non superpartitioned
        bucket = MockBucket("")

        # test superpartitioned
        avro_prefix = "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro"
        avro_dict = generate_avro_args(
            avro_prefix=avro_prefix,
            gcs_listing=listing
        )
        actual = generate_avro_args(avro_dict)
        self.assertEqual(actual.strip(), expected.strip())
