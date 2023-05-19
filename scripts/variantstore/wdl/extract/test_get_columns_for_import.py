import unittest
import json


from get_columns_for_import import get_column_values


class TestBulkIngestGenomes(unittest.TestCase):

    def test_get_column_values(self):
        numSamples = 5
        with open('bulk_ingest_test_files/columns_for_import.json') as columnSamples:
            columnSamplesExpected = json.load(columnSamples)
            expected = {
                'vcf_column' : 'hg38_reblocked_v2_vcf',
                'vcf_column_index': 'hg38_reblocked_v2_vcf_index'
            }

            actual = get_column_values(columnSamplesExpected, numSamples)
            self.assertEqual(actual, expected)

    def test_get_column_quickstart_values(self):
        numSamples = 10
        with open('bulk_ingest_test_files/quickstart_columns_for_import.json') as quickstartColumnSamples:
            columnSamplesExpected = json.load(quickstartColumnSamples)
            expected = {
                'vcf_column' : 'hg38_reblocked_v2_vcf',
                'vcf_column_index': 'hg38_reblocked_v2_vcf_index'
            }

            actual = get_column_values(columnSamplesExpected, numSamples)
            self.assertEqual(actual, expected)

    def test_get_column_aou_values(self):
        numSamples = 50
        ## note that external_sample_names is research_id
        with open('bulk_ingest_test_files/aou_columns_for_import.json') as aouColumnSamples:
            columnSamplesExpected = json.load(aouColumnSamples)

            expected = {
                'vcf_column' : 'reblocked_gvcf',
                'vcf_column_index': 'reblocked_gvcf_index'
            }

            actual = get_column_values(columnSamplesExpected, numSamples)
            self.assertEqual(actual, expected)


    def test_get_column_shriners_values(self):
        numSamples = 20
        with open('bulk_ingest_test_files/shriners_columns_for_import.json') as shrinersColumnSamples:
            columnSamplesExpected = json.load(shrinersColumnSamples)
            expected = {
                'vcf_column' : 'gvcf',
                'vcf_column_index': 'gvcf_index'
            }

            actual = get_column_values(columnSamplesExpected, numSamples)
            self.assertEqual(actual, expected)
