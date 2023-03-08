import argparse

from google.cloud import bigquery
from google.cloud.bigquery.job import QueryJobConfig
import utils


def run_avro_query(call_set_identifier, dataset_name, table_name, project_id, sql):
    # add labels for DSP Cloud Cost Control Labeling and Reporting
    query_labels_map = {'service': 'gvs', 'team': 'variants', 'managedby': 'gvs_extract_avro_files_for_hail'}

    default_config = QueryJobConfig(labels=query_labels_map, priority="INTERACTIVE", use_query_cache=True)
    client = bigquery.Client(project=project_id,
                             default_query_job_config=default_config)
    query_return = utils.execute_with_retry(client, table_name, sql)
    utils.write_job_stats([{'job': query_return['job'], 'label': query_return['label']}], client,
                          f"{project_id}.{dataset_name}", call_set_identifier, 'GvsExtractAvroFilesForHail',
                          'GenerateAvroFiles', table_name, True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False, description='Extract avro files from GVS BigQuery dataset')
    parser.add_argument('--call_set_identifier', type=str,
                        help='callset identifier used to track costs in cost_observability table', default='false')
    parser.add_argument('--dataset_name',type=str, help='BigQuery dataset name', required=True)
    parser.add_argument('--table_name',type=str, help='BigQuery table name', required=True)
    parser.add_argument('--project_id', type=str, help='Google project for the GVS dataset', required=True)
    parser.add_argument('--sql', type=str, help='SQL to run to extract Avro data', required=True)

    args = parser.parse_args()

    run_avro_query(args.call_set_identifier,
                        args.dataset_name,
                        args.table_name,
                        args.project_id,
                        args.sql)
