# -*- coding: utf-8 -*-
import math
import argparse

from google.cloud import bigquery
from google.cloud.bigquery.job import QueryJobConfig

import run_avro_query
import utils

## from create_ranges_cohort_extract_data_table import load_sample_names ## TODO can I just import this?!?!

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

## load_sample_names(sample_names_to_extract, fq_temp_table_dataset) TODO can I just import this?!?!

def load_sample_names(sample_names_to_extract, fq_temp_table_dataset):
    schema = [bigquery.SchemaField("sample_name", "STRING", mode="REQUIRED")]
    fq_sample_table = f"{fq_temp_table_dataset}.{EXTRACT_SAMPLE_TABLE}"

    job_labels = client._default_query_job_config.labels
    job_labels["gvs_query_name"] = "load-sample-names"

    job_config = bigquery.LoadJobConfig(source_format=bigquery.SourceFormat.CSV, skip_leading_rows=0, schema=schema,
                                        labels=job_labels)

    with open(sample_names_to_extract, "rb") as source_file:
        job = client.load_table_from_file(source_file, fq_sample_table, job_config=job_config)

    job.result()  # Waits for the job to complete.

    # setting the TTL needs to be done as a second API call
    table = bigquery.Table(fq_sample_table, schema=schema)
    expiration = datetime.datetime.utcnow() + datetime.timedelta(hours=TEMP_TABLE_TTL_HOURS)
    table.expires = expiration
    client.update_table(table, ["expires"])

    return fq_sample_table


def construct_sample_info_avro_queries(call_set_identifier, dataset_name, project_id, avro_prefix):
    num_of_tables = get_number_of_partitions(dataset_name, project_id)

   # if we have a file of sample names, load it into a temporary table
   if sample_names_to_extract:
       ## TODO note that I haven't defined sample_names_to_extract OR fq_temp_table_dataset
       fq_sample_name_table = load_sample_names(sample_names_to_extract, fq_temp_table_dataset)
   else:
       fq_sample_name_table = f"{project_id}.{dataset_name}.sample_info"

    for i in range(1, num_of_tables + 1):
        file_name = f"*.{i:03}.avro"
        id_where_clause = f"sample_id >= {((i -1) * 4000) + 1} AND sample_id <= {i * 4000}"

        sql = f"""
            EXPORT DATA OPTIONS(
                uri='{avro_prefix}/sample_mapping/{file_name}', format='AVRO', compression='SNAPPY', overwrite=true) AS
            SELECT sample_id, sample_name, '40',
            'gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list' AS intervals_file
            FROM `{fq_sample_name_table}`
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

    args = parser.parse_args()

    construct_sample_info_avro_queries(args.call_set_identifier,
                   args.dataset_name,
                   args.project_id,
                   args.avro_prefix)
