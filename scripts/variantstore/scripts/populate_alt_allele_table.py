import os
import argparse
import sys

from google.cloud import bigquery
from google.cloud.bigquery.job import QueryJobConfig
from pathlib import Path

import utils

client = None


def populate_alt_allele_table(call_set_identifier, query_project, vet_table_name, fq_dataset, max_sample_id):
    global client
    # add labels for DSP Cloud Cost Control Labeling and Reporting to default_config
    default_config = QueryJobConfig(priority="INTERACTIVE", use_query_cache=True, use_legacy_sql=False,
                                    labels={'service': 'gvs', 'team': 'variants', 'managedby': 'create_alt_allele'})

    client = bigquery.Client(project=query_project,
                             default_query_job_config=default_config)

    os.chdir(os.path.dirname(__file__))
    alt_allele_temp_function = Path('alt_allele_temp_function.sql').read_text()
    alt_allele_positions = Path('alt_allele_positions.sql').read_text()
    fq_vet_table = f"{fq_dataset}.{vet_table_name}"
    fq_sample_info = f"{fq_dataset}.sample_info"
    fq_alt_allele = f"{fq_dataset}.alt_allele"

    query_with = f"""

    INSERT INTO `{fq_dataset}.alt_allele`
    WITH
        position1 AS (
            SELECT * FROM `{fq_vet_table}` WHERE
                call_GT IN ('0/1', '1/0', '1/1', '0|1', '1|0', '1|1', '0/2', '0|2','2/0', '2|0', '1', '2') AND
                sample_id > {max_sample_id} AND
                sample_id IN (SELECT sample_id from {fq_sample_info} WHERE withdrawn is NULL)
        ),
        position2 AS (
            SELECT * FROM `{fq_vet_table}` WHERE
                call_GT IN ('1/2', '1|2', '2/1', '2|1') AND
                sample_id > {max_sample_id} AND
                sample_id IN (SELECT sample_id from {fq_sample_info} WHERE withdrawn is NULL)
        )
    """

    sql = alt_allele_temp_function + query_with + alt_allele_positions
    query_return = utils.execute_with_retry(client, f"into alt allele from {vet_table_name}", sql)
    utils.write_job_stats([{'job': query_return['job'], 'label': query_return['label']}], client, fq_dataset,
                          call_set_identifier, 'CreateAltAlleles', 'PopulateAltAlleleTable', vet_table_name)

    # Verify that the number of distinct samples is correct
    verify_distinct_samples(fq_alt_allele, fq_sample_info)

    return query_return['results']


def verify_distinct_samples(fq_alt_allele_table, fq_sample_info_table):
    """
    Verify that the number of distinct samples in the alt_allele table
    matches the number of samples in the sample_info table for whom withdrawn is NULL.

    Args:
        fq_alt_allele_table: Fully qualified alt_allele table name
        fq_sample_info_table: Fully qualified sample_info table name
    """
    # Query to count distinct samples in alt_allele
    alt_allele_query = f"""
    SELECT COUNT(DISTINCT sample_id) AS alt_allele_sample_count
    FROM `{fq_alt_allele_table}`
    """

    # Query to count non-withdrawn samples in sample_info
    sample_info_query = f"""
    SELECT COUNT(sample_id) AS sample_info_count
    FROM `{fq_sample_info_table}`
    WHERE withdrawn IS NULL
    """

    # Execute queries
    alt_allele_count_job = client.query(alt_allele_query)
    sample_info_count_job = client.query(sample_info_query)

    # Get results
    alt_allele_count = list(alt_allele_count_job)[0].alt_allele_sample_count
    sample_info_count = list(sample_info_count_job)[0].sample_info_count

    # Compare counts and raise error if they don't match
    if alt_allele_count != sample_info_count:
        print(f"ERROR: Sample count mismatch! Alt allele table contains {alt_allele_count} distinct samples, "
              f"but sample_info table contains {sample_info_count} non-withdrawn samples.")
        sys.exit(1)
    else:
        print(f"Verification passed: Both tables contain {sample_info_count} samples.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False,
                                     description='Populate an alt_allele table for the BigQuery Variant Store')

    parser.add_argument('--call_set_identifier', type=str,
                        help='callset identifier used to track costs in cost_observability table', default='false')
    parser.add_argument('--query_project', type=str, help='Google project where query should be executed',
                        required=True)
    parser.add_argument('--vet_table_name', type=str, help='vet table name to ingest', required=True)
    parser.add_argument('--fq_dataset', type=str, help='project and dataset for data', required=True)
    parser.add_argument('--max_sample_id',type=str, help='Maximum value of sample_id already loaded', required=True)

    args = parser.parse_args()

    populate_alt_allele_table(args.call_set_identifier,
                              args.query_project,
                              args.vet_table_name,
                              args.fq_dataset,
                              args.max_sample_id)
