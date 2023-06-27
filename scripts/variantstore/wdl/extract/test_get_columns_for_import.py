import unittest
import json


from get_columns_for_import import get_column_values


class TestBulkIngestGenomes(unittest.TestCase):

    def test_get_column_values(self):
        numSamples = 5
        with open('bulk_ingest_test_files/columns_for_import.json') as columnSamples:
            columnSamplesExpected = json.load(columnSamples)
            expected = ('hg38_reblocked_v2_vcf', 'hg38_reblocked_v2_vcf_index')
            actual = get_column_values(columnSamplesExpected, numSamples, None, None)
            self.assertEqual(actual, expected)

    def test_get_column_user_defined_values(self):
        numSamples = 5
        with open('bulk_ingest_test_files/columns_for_import.json') as columnSamples:
            columnSamplesExpected = json.load(columnSamples)
            expected = ("control_vcf", "control_vcf_index")
            actual = get_column_values(columnSamplesExpected, numSamples, "control_vcf", "control_vcf_index")
            self.assertEqual(actual, expected)

    def test_get_column_bad_user_defined_values(self):
        numSamples = 5
        with open('bulk_ingest_test_files/columns_for_import.json') as columnSamples:
            columnSamplesExpected = json.load(columnSamples)
            with self.assertRaises(ValueError):
                get_column_values(columnSamplesExpected, numSamples, "im_the_problem", "its_me")

    def test_get_column_missing_user_defined_vcf(self):
        numSamples = 5
        with open('bulk_ingest_test_files/columns_for_import.json') as columnSamples:
            columnSamplesExpected = json.load(columnSamples)
            with self.assertRaises(ValueError):
                get_column_values(columnSamplesExpected, numSamples, None, "control_vcf_index")

    def test_get_column_unspecified_user_defined_index(self):
        numSamples = 5
        with open('bulk_ingest_test_files/columns_for_import.json') as columnSamples:
            columnSamplesExpected = json.load(columnSamples)
            expected = ("control_vcf", "control_vcf_index")
            actual = get_column_values(columnSamplesExpected, numSamples, "control_vcf", None)
            self.assertEqual(actual, expected)

    def test_get_column_missing_user_defined_index(self):
        numSamples = 5
        with open('bulk_ingest_test_files/columns_for_import.json') as columnSamples:
            columnSamplesExpected = json.load(columnSamples)
            with self.assertRaises(ValueError):
                get_column_values(columnSamplesExpected, numSamples, "random_vcf", None)


    def test_get_column_quickstart_values(self):
        numSamples = 10
        with open('bulk_ingest_test_files/quickstart_columns_for_import.json') as quickstartColumnSamples:
            columnSamplesExpected = json.load(quickstartColumnSamples)
            expected = ('hg38_reblocked_v2_vcf', 'hg38_reblocked_v2_vcf_index')
            actual = get_column_values(columnSamplesExpected, numSamples, None, None)
            self.assertEqual(actual, expected)

    def test_get_column_aou_values(self):
        numSamples = 50
        ## note that external_sample_names is research_id
        with open('bulk_ingest_test_files/aou_columns_for_import.json') as aouColumnSamples:
            columnSamplesExpected = json.load(aouColumnSamples)
            expected = ('reblocked_gvcf', 'reblocked_gvcf_index')
            actual = get_column_values(columnSamplesExpected, numSamples, None, None)
            self.assertEqual(actual, expected)


    def test_get_column_shriners_values(self):
        numSamples = 20
        with open('bulk_ingest_test_files/shriners_columns_for_import.json') as shrinersColumnSamples:
            columnSamplesExpected = json.load(shrinersColumnSamples)
            expected = ('gvcf','gvcf_index')
            actual = get_column_values(columnSamplesExpected, numSamples, None, None)
            self.assertEqual(actual, expected)


