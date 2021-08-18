from google.cloud import bigquery
from google.cloud.bigquery.job import QueryJobConfig
from google.oauth2 import service_account
from pathlib import Path

import argparse

client = None

def execute_with_retry(label, sql):
    retry_delay = [30, 60, 90] # 3 retries with incremental backoff
    start = time.time()
    while len(retry_delay) > 0:
        try:
            query_label = label.replace(" ","-").strip().lower()

            existing_labels = client._default_query_job_config.labels
            job_labels = existing_labels
            job_labels["gvs_query_name"] = query_label
            job_config = bigquery.QueryJobConfig(labels=job_labels)
            query = client.query(sql, job_config=job_config)

            print(f"STARTING - {label}")
            JOB_IDS.add((label, query.job_id))
            results = query.result()
            print(f"COMPLETED ({time.time() - start} s, {3-len(retry_delay)} retries) - {label}")
            return results
        except Exception as err:
            # if there are no retries left... raise
            if (len(retry_delay) == 0):
                raise err
            else:
                t = retry_delay.pop(0)
                print(f"Error {err} running query {label}, sleeping for {t}")
                time.sleep(t)

def make_alt_allele_table(query_project, vet_table_name, fq_dataset, sa_key_path):
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

    alt_allele_temp_function = Path('alt_allele_temp_function.sql').read_text()
    alt_allele_positions = Path('alt_allele_positions.sql').read_text()
    first = True if vet_table_name == "vet_001" else False

    query_beginning = f"CREATE OR REPLACE TABLE {fq_dataset}.alt_allele PARTITION BY \
                    RANGE_BUCKET(location, GENERATE_ARRAY(0, 25000000000000, 1000000000000)) \
                    CLUSTER BY location, sample_id AS \n"
    if not first:
        query_beginning = f"INSERT INTO {fq_dataset}.alt_allele \m"
    fq_vet_table = f"~{fq_dataset}.{vet_table_name}"
    query_with = f"""WITH 
                  position1 as (select * from {fq_vet_table} WHERE call_GT IN ('0/1', '1/0', '1/1', '0|1', '1|0', '1|1', '0/2', '0|2','2/0', '2|0')), 
                  position2 as (select * from {fq_vet_table} WHERE call_GT IN ('1/2', '1|2', '2/1', '2|1'))"""

    sql = alt_allele_temp_function + query_beginning + query_with + alt_allele_positions
    result = execute_with_retry(f"into alt allele from {vet_table_name}", sql)
    return result

if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False, description='Populate an alt_allele table for the BigQuery Variant Store')

    parser.add_argument('--query_project',type=str, help='Google project where query should be executed', required=True)
    parser.add_argument('--vet_table_name',type=str, help='vet table name to ingest', required=True)
    parser.add_argument('--fq_dataset',type=str, help='project and dataset for data', required=True)
    parser.add_argument('--sa_key_path',type=str, help='Path to json key file for SA', required=False)


    # Execute the parse_args() method
    args = parser.parse_args()

    make_alt_allele_table(args.query_project, args.vet_table_name, args.fq_dataset, args.sa_key_path)
