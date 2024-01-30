# -*- coding: utf-8 -*-
import argparse
import ijson
import re

from google.cloud import bigquery
from google.cloud.bigquery.job import QueryJobConfig

import utils



def write_output_files(fq_vat_table, empty_columns, pass_file_output, results_file_output):
    with open(pass_file_output, 'w') as pass_file,  open(results_file_output, 'w') as results_file :
        if empty_columns != []:
            pass_file.write('false')
            results_file.write(f"The following columns in {fq_vat_table} are empty: ({', '.join(empty_columns)})")
        else:
            pass_file.write('true')
            results_file.write(f"No nullable {fq_vat_table} columns are completely empty")
    pass_file.close()
    results_file.close()


def get_fields_to_check(schema_file):
    fields_to_check = []
    with open(schema_file, 'r') as input_file:
        items = ijson.items(input_file, '', use_float=True)
        schema = items.__next__();
        for field in schema:
            if field.get('mode') == "Nullable":
                fields_to_check.append(field.get('name'))
    return fields_to_check


def check_for_null_columns(fq_vat_table, query_project, schema_file_input, pass_file_output, results_file_output):
    query_labels_map = {
        "id": "check_for_null_columns",
        "gvs_tool_name": "gvs_validate_vat"
    }

    # add labels for DSP Cloud Cost Control Labeling and Reporting
    query_labels_map.update({'service': 'gvs', 'team': 'variants', 'managedby': 'prepare_ranges_callset'})

    # Default QueryJobConfig will be merged into job configs passed in
    # but if a specific default config is being updated (eg labels), new config must be added
    # to the client._default_query_job_config that already exists
    default_config = QueryJobConfig(labels=query_labels_map, priority="INTERACTIVE", use_query_cache=True)
    client = bigquery.Client(project=query_project,
                             default_query_job_config=default_config)

    fields_to_check = get_fields_to_check(schema_file_input)

    select_statements = [f"(SELECT '{field}' AS column, (SELECT count({field}) FROM `{fq_vat_table}` WHERE {field} is NOT NULL) AS number)"
                         for field in fields_to_check]
    all_fields = " UNION ALL ".join(select_statements)
    sql = f"SELECT column FROM ({all_fields}) AS all_fields where number = 0"

    query_return = utils.execute_with_retry(client, f"get null fields from VAT", sql)

    empty_columns = [row.get('column') for row in query_return.get('results')]
    write_output_files(fq_vat_table, empty_columns, pass_file_output, results_file_output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False,
                                     description='Check variant annotations table for empty columns')
    parser.add_argument('--fq_vat_table', type=str, required=True,
                        help='fully-qualified BigQuery VAT table')
    parser.add_argument('--query_project', type=str, required=True,
                        help='Google project to run the BigQuery queries with')
    parser.add_argument('--schema_file_input', type=str, required=True,
                        help='latest VAT schema JSON')
    parser.add_argument('--pass_file_output', type=str, required=True,
                        help='location to write file with pass/fail into')
    parser.add_argument('--results_file_output', type=str, required=True,
                        help='location to write query results into')

    args = parser.parse_args()

    check_for_null_columns(args.fq_vat_table,
                           args.query_project,
                           args.schema_file_input,
                           args.pass_file_output,
                           args.results_file_output)
