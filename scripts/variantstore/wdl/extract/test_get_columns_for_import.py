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

## columns like the following are choking the python. Not sure how they would come out of the url tho--maybe the issue is that we never expect an object?
#     {"biobank_id":"A669839009","crai_path":"gs://prod-genomics-data-broad/wgs_sample_raw_data/crams_crais/BI_A669839009_22045006567_0389692910_1.cram.crai","cram_md5_path":"gs://prod-genomics-data-broad/wgs_sample_raw_data/crams_crais/BI_A669839009_22045006567_0389692910_1.cram.md5sum","cram_path":"gs://prod-genomics-data-broad/wgs_sample_raw_data/crams_crais/BI_A669839009_22045006567_0389692910_1.cram","drc_contamination":"1.24479E-5","drc_fp_concordance":"PASS","drc_mean_coverage":"39.631482","drc_sex_concordance":"PASS","gvcf_md5_path":"gs://prod-genomics-data-broad/wgs_sample_raw_data/ss_vcf_research/BI_A669839009_22045006567_0389692910_1.hard-filtered.gvcf.gz.md5sum","gvcf_path":"gs://prod-genomics-data-broad/wgs_sample_raw_data/ss_vcf_research/BI_A669839009_22045006567_0389692910_1.hard-filtered.gvcf.gz","pass_to_research_pipeline":"True","qc_status":"PASS","research_id":"1000093","sample_id":"22045006567","sex_at_birth":"F","site_id":"bi","vcf_hf_index_path":"gs://prod-genomics-data-broad/wgs_sample_raw_data/ss_vcf_clinical/BI_A669839009_22045006567_0389692910_1.hard-filtered.vcf.gz.tbi","vcf_hf_md5_path":"gs://prod-genomics-data-broad/wgs_sample_raw_data/ss_vcf_clinical/BI_A669839009_22045006567_0389692910_1.hard-filtered.vcf.gz.md5sum","vcf_hf_path":"gs://prod-genomics-data-broad/wgs_sample_raw_data/ss_vcf_clinical/BI_A669839009_22045006567_0389692910_1.hard-filtered.vcf.gz","vcf_raw_index_path":"","vcf_raw_md5_path":"","vcf_raw_path":""},
#     {"biobank_id":"A669839009","biobankidsampleid":"A66983900922045006567","crai":"D","cram":"D","cram_basename":"BI_A669839009_22045006567_0389692910_1.cram","cram_md5":"D","cram_md5_hash":"2c8efc8fcbb728892f72483afe143a7a","gvcf":"D","gvcf_basename":"BI_A669839009_22045006567_0389692910_1.hard-filtered.gvcf.gz","gvcf_md5":"D","gvcf_md5_hash":"162b4f7ef0b4768e59259ae0ce471762","sample_id":"22045006567","sex_at_birth":"F","site_id":"bi","vcf_hf":"D","vcf_hf_basename":"BI_A669839009_22045006567_0389692910_1.hard-filtered.vcf.gz","vcf_hf_index":"D","vcf_hf_md5":"D","vcf_hf_md5_hash":"c56557f012f08a62206a30cef2afa05c","vcf_raw":"ND","vcf_raw_basename":"","vcf_raw_index":"ND","vcf_raw_md5":"ND","vcf_raw_md5_hash":""},
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
