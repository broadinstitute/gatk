import json
import unittest

from generate_hail_gvs_import import generate_avro_dict, generate_avro_args


class TestGenerateAvroDict(unittest.TestCase):

    def test_generate_avro_dict(self):
        actual = generate_avro_dict(
            avro_prefix="gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro",
            listing="generate_hail_gvs_import_test_files/avro_path_listing.txt"
        )
        args = []
        for k, v in actual.items():
            s = f'{k}: {json.dumps(v, indent=4)}'
            args.append(s)
        a = ',\n'.join(args)
        print(a)
        self.assertEqual(actual, None)

    def test_generate_avro_args(self):
        avro_dict = generate_avro_dict(
            avro_prefix="gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro",
            listing="generate_hail_gvs_import_test_files/avro_path_listing.txt"
        )
        x = generate_avro_args(avro_dict)
        print(x)


