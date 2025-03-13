# -*- coding: utf-8 -*-
import math
import argparse

from google.cloud import bigquery
from google.cloud.bigquery.job import QueryJobConfig

import run_avro_query
import utils


def get_number_of_partitions(dataset_name, project_id):
    query_labels_map = {
        "id": "construct_sample_info_avro_queries",
        "gvs_tool_name": "gvs_extract_avro_files_for_hail",
        "service": "gvs",
        "team": "variants",
        "managedby": "gvs_extract_avro_files_for_hail"
    }

    # Default QueryJobConfig will be merged into job configs passed in
    # but if a specific default config is being updated (eg labels), new config must be added
    # to the client._default_query_job_config that already exists
    default_config = QueryJobConfig(labels=query_labels_map, priority="INTERACTIVE", use_query_cache=True)
    client = bigquery.Client(project=project_id,
                             default_query_job_config=default_config)

    sql = f"SELECT max(sample_id) / 4000 as max_table_num FROM `{project_id}.{dataset_name}.sample_info`"

    query_result = utils.execute_with_retry(client, f"get max partitioned tabled num", sql)
    max_table_num = [row.get('max_table_num') for row in query_result.get('results')]
    return math.ceil(max_table_num[0])


def construct_sample_info_avro_queries(call_set_identifier, dataset_name, project_id, avro_prefix, new_sample_cutoff=None):
    num_of_tables = get_number_of_partitions(dataset_name, project_id)

    # For Foxtrot, we need to start at the new_sample_cutoff and export everything remaining in that group and then export the sample_info table in 4000 sample chunks to match with the other tables
    ## we want to know what table number to start with if we have a new sample cutoff
    ## and we will want groupings of 4000 samples once we have made it from the new sample cutoff to the next group start
    ## so we will need one avro file of everything from the int(new_sample_cutoff) to the end of the next table, and then the rest of the tables
    if new_sample_cutoff:
        ## set the loop down below to begin at the group after the cutoff
        start_table = ceil(int(new_sample_cutoff) // 4000)
        ## run a query to get the samples from the cutoff to the end of the table
        first_file = start_table - 1
        file_name = f"*.{first_file:03}.avro"
        sql = f"""
            EXPORT DATA OPTIONS(
                uri='{avro_prefix}/sample_mapping/{file_name}', format='AVRO', compression='SNAPPY', overwrite=true) AS
            SELECT sample_id, sample_name, '40',
            'gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list' AS intervals_file
            FROM `{project_id}.{dataset_name}.sample_info`
            WHERE sample_id > {new_sample_cutoff} AND sample_id <= {start_table * 4000}
            AND is_control = false
            ORDER BY sample_id"""

    else:
        start_table = 1


    for i in range(start_table, num_of_tables + 1):
        file_name = f"*.{i:03}.avro"
        id_where_clause = f"sample_id >= {((i -1) * 4000) + 1} AND sample_id <= {i * 4000}"
        sql = f"""
            EXPORT DATA OPTIONS(
                uri='{avro_prefix}/sample_mapping/{file_name}', format='AVRO', compression='SNAPPY', overwrite=true) AS
            SELECT sample_id, sample_name, '40',
            'gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list' AS intervals_file
            FROM `{project_id}.{dataset_name}.sample_info`
            WHERE {id_where_clause} AND is_control = false
            ORDER BY sample_id"""

        print(f"{sql}\n")
        run_avro_query.run_avro_query(call_set_identifier, dataset_name, "sample_info", project_id, sql)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False, description='Construct and run queries to generate sample avro files from GVS BigQuery dataset')
    parser.add_argument('--call_set_identifier', type=str,
                        help='callset identifier used to track costs in cost_observability table', required=True)
    parser.add_argument('--dataset_name',type=str, help='BigQuery dataset name', required=True)
    parser.add_argument('--project_id', type=str, help='Google project for the GVS dataset', required=True)
    parser.add_argument('--avro_prefix', type=str, help='prefix for the Avro file path', required=True)
    parser.add_argument('--new_sample_cutoff', type=str, help='cutoff for sample ids')

    args = parser.parse_args()

    if (args.new_sample_cutoff):
            construct_sample_info_avro_queries(args.call_set_identifier,
                           args.dataset_name,
                           args.project_id,
                           args.avro_prefix,
                           args.new_sample_cutoff)
    else:
        construct_sample_info_avro_queries(args.call_set_identifier,
                   args.dataset_name,
                   args.project_id,
                   args.avro_prefix)
