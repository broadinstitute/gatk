import argparse

from google.cloud import bigquery
from google.cloud.bigquery.job import QueryJobConfig
import utils

# add labels for DSP Cloud Cost Control Labeling and Reporting
query_labels_map = {'service': 'gvs', 'team': 'variants', 'managedby': 'gvs_process_sample_vcf_headers'}
default_config = QueryJobConfig(labels=query_labels_map, priority="INTERACTIVE", use_query_cache=True)
client = ''

def process_sample_vcf_headers(project_id, dataset_name):
    populate_tables_from_scratch(project_id, dataset_name)
    clean_up_scratch_table(project_id, dataset_name)

def populate_tables_from_scratch(project_id, dataset_name):
    global client
    client = bigquery.Client(project=project_id,
                             default_query_job_config=default_config)

    sql = f"INSERT INTO {project_id}.{dataset_name}.vcf_header_lines (vcf_header_lines_hash, vcf_header_lines, is_expected_unique) SELECT vcf_header_lines_hash, vcf_header_lines, is_expected_unique FROM {project_id}.{dataset_name}.vcf_header_lines_scratch WHERE vcf_header_lines IS NOT NULL"
    utils.execute_with_retry(client, "vcf_header_lines", sql)

    sql = f"INSERT INTO {project_id}.{dataset_name}.sample_vcf_header (sample_id, vcf_header_lines_hash) SELECT sample_id, vcf_header_lines_hash FROM {project_id}.{dataset_name}.vcf_header_lines_scratch"
    utils.execute_with_retry(client, "sample_vcf_header", sql)


def clean_up_scratch_table(project_id, dataset_name):
    global client
    sql = f"DELETE FROM {project_id}.{dataset_name}.vcf_header_lines_scratch WHERE vcf_header_lines_hash IS NOT NULL"
    utils.execute_with_retry(client, "vcf_header_lines_scratch", sql)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False, description='')
    parser.add_argument('--project_id', type=str, help='Google project for the GVS dataset', required=True)
    parser.add_argument('--dataset_name',type=str, help='BigQuery dataset name', required=True)

    args = parser.parse_args()

    process_sample_vcf_headers(args.project_id, args.dataset_name)
