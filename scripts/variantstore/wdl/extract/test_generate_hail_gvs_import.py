import unittest

from hail_gvs_import import generate_avro_args


class TestGenerateAvroDict(unittest.TestCase):

    def test_generate_avro_args(self):
        expected = """
refs=[
    [
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/refs/ref_ranges_001/ref_ranges_001_000000000000.avro",
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/refs/ref_ranges_001/ref_ranges_001_000000000001.avro",
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/refs/ref_ranges_001/ref_ranges_001_000000000002.avro",
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/refs/ref_ranges_001/ref_ranges_001_000000000003.avro",
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/refs/ref_ranges_001/ref_ranges_001_000000000004.avro",
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/refs/ref_ranges_001/ref_ranges_001_000000000005.avro",
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/refs/ref_ranges_001/ref_ranges_001_000000000006.avro",
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/refs/ref_ranges_001/ref_ranges_001_000000000007.avro",
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/refs/ref_ranges_001/ref_ranges_001_000000000008.avro",
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/refs/ref_ranges_001/ref_ranges_001_000000000009.avro",
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/refs/ref_ranges_001/ref_ranges_001_000000000010.avro",
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/refs/ref_ranges_001/ref_ranges_001_000000000011.avro",
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/refs/ref_ranges_001/ref_ranges_001_000000000012.avro",
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/refs/ref_ranges_001/ref_ranges_001_000000000013.avro",
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/refs/ref_ranges_001/ref_ranges_001_000000000014.avro",
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/refs/ref_ranges_001/ref_ranges_001_000000000015.avro",
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/refs/ref_ranges_001/ref_ranges_001_000000000016.avro",
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/refs/ref_ranges_001/ref_ranges_001_000000000017.avro",
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/refs/ref_ranges_001/ref_ranges_001_000000000018.avro",
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/refs/ref_ranges_001/ref_ranges_001_000000000019.avro",
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/refs/ref_ranges_001/ref_ranges_001_000000000020.avro",
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/refs/ref_ranges_001/ref_ranges_001_000000000021.avro",
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/refs/ref_ranges_001/ref_ranges_001_000000000022.avro",
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/refs/ref_ranges_001/ref_ranges_001_000000000023.avro",
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/refs/ref_ranges_001/ref_ranges_001_000000000024.avro",
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/refs/ref_ranges_001/ref_ranges_001_000000000025.avro",
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/refs/ref_ranges_001/ref_ranges_001_000000000026.avro",
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/refs/ref_ranges_001/ref_ranges_001_000000000027.avro",
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/refs/ref_ranges_001/ref_ranges_001_000000000028.avro"
    ]
],
sample_mapping_data=[
    "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/sample_mapping_data/sample_mapping_data_000000000000.avro"
],
site_filtering_data=[
    "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/site_filtering_data/site_filtering_data_000000000000.avro"
],
vets=[
    [
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/vets/vet_001/vet_001_000000000000.avro",
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/vets/vet_001/vet_001_000000000001.avro",
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/vets/vet_001/vet_001_000000000002.avro",
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/vets/vet_001/vet_001_000000000003.avro",
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/vets/vet_001/vet_001_000000000004.avro",
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/vets/vet_001/vet_001_000000000005.avro",
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/vets/vet_001/vet_001_000000000006.avro",
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/vets/vet_001/vet_001_000000000007.avro",
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/vets/vet_001/vet_001_000000000008.avro",
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/vets/vet_001/vet_001_000000000009.avro",
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/vets/vet_001/vet_001_000000000010.avro",
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/vets/vet_001/vet_001_000000000011.avro",
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/vets/vet_001/vet_001_000000000012.avro",
        "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/vets/vet_001/vet_001_000000000013.avro"
    ]
],
vqsr_filtering_data=[
    "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/vqsr_filtering_data/vqsr_filtering_data_000000000000.avro",
    "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/vqsr_filtering_data/vqsr_filtering_data_000000000001.avro",
    "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/vqsr_filtering_data/vqsr_filtering_data_000000000002.avro"
],
vqsr_tranche_data=[
    "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro/vqsr_tranche_data/vqsr_tranche_data_000000000000.avro"
]
"""
        with open('generate_hail_gvs_import_test_files/avro_path_listing.txt', 'r') as listing:
            avro_prefix = "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/submissions/748efa9e-f176-47ac-8a93-eff632fc6b8a/GvsExtractAvroFilesForHail/c9fbac0f-9c49-46f5-92a4-fc468e8fa94f/call-OutputPath/avro"
            avro_dict = generate_avro_args(
                avro_prefix=avro_prefix,
                gcs_listing=listing
            )
            actual = generate_avro_args(avro_dict)
            self.assertEqual(actual.strip(), expected.strip())
