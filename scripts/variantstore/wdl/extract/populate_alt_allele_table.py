import os
import argparse

from google.cloud import bigquery
from google.cloud.bigquery.job import QueryJobConfig
from google.oauth2 import service_account
from pathlib import Path

import utils

client = None


def populate_alt_allele_table(query_project, vet_table_name, fq_dataset, sa_key_path):
    global client
    default_config = QueryJobConfig(priority="INTERACTIVE", use_query_cache=True)

    if sa_key_path:
        credentials = service_account.Credentials.from_service_account_file(
            sa_key_path, scopes=["https://www.googleapis.com/auth/cloud-platform"],
        )

        client = bigquery.Client(credentials=credentials,
                                 project=query_project,
                                 default_query_job_config=default_config)
    else:
        client = bigquery.Client(project=query_project,
                                 default_query_job_config=default_config)

    os.chdir(os.path.dirname(__file__))
    alt_allele_temp_function = Path('alt_allele_temp_function.sql').read_text()
    alt_allele_positions = Path('alt_allele_positions.sql').read_text()
    fq_vet_table = f"{fq_dataset}.{vet_table_name}"
    query_with = f"""INSERT INTO `{fq_dataset}.alt_allele`
                WITH 
                  position1 as (select * from `{fq_vet_table}` WHERE call_GT IN ('0/1', '1/0', '1/1', '0|1', '1|0', '1|1', '0/2', '0|2','2/0', '2|0')),
                  position2 as (select * from `{fq_vet_table}` WHERE call_GT IN ('1/2', '1|2', '2/1', '2|1'))"""

    sql = alt_allele_temp_function + query_with + alt_allele_positions
    result = utils.execute_with_retry(client, f"into alt allele from {vet_table_name}", sql)
    return result

if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False, description='Populate an alt_allele table for the BigQuery Variant Store')

    parser.add_argument('--query_project',type=str, help='Google project where query should be executed', required=True)
    parser.add_argument('--vet_table_name',type=str, help='vet table name to ingest', required=True)
    parser.add_argument('--fq_dataset',type=str, help='project and dataset for data', required=True)
    parser.add_argument('--sa_key_path',type=str, help='Path to json key file for SA', required=False)


    # Execute the parse_args() method
    args = parser.parse_args()

    populate_alt_allele_table(args.query_project, args.vet_table_name, args.fq_dataset, args.sa_key_path)
