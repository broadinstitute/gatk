import os
import argparse

from google.cloud import bigquery
from google.cloud.bigquery.job import QueryJobConfig
from pathlib import Path
from string import Template

import utils

client = None


def alt_allele_positions_sql(template_text, include_vrs_allele_id):
    """Fill the ${vrs_allele_id_altN} placeholders in alt_allele_positions.sql.

    alt_allele has a vrs_allele_id column only in a dataset whose vet tables have vrs_allele_ids, and the INSERT is
    positional, so the column is selected only then. It is the vet.vrs_allele_ids entry for the row's ALT, at the same
    offset the row reads from SPLIT(alt, ','). NULLIF turns the empty-string placeholder that ingest writes for a `*`
    ALT, or for an ALT on a contig with no refget accession, into NULL.
    """
    def column(offset):
        if not include_vrs_allele_id:
            return ""
        return f", NULLIF(vrs_allele_ids[SAFE_OFFSET({offset})], '') as vrs_allele_id"

    return Template(template_text).substitute(vrs_allele_id_alt1=column(0), vrs_allele_id_alt2=column(1))


def populate_alt_allele_table(call_set_identifier, query_project, vet_table_name, fq_dataset, max_sample_id,
                              include_vrs_allele_id):
    global client
    # add labels for DSP Cloud Cost Control Labeling and Reporting to default_config
    default_config = QueryJobConfig(priority="INTERACTIVE", use_query_cache=True, use_legacy_sql=False,
                                    labels={'service': 'gvs', 'team': 'variants', 'managedby': 'create_alt_allele'})

    client = bigquery.Client(project=query_project,
                             default_query_job_config=default_config)

    os.chdir(os.path.dirname(__file__))
    alt_allele_temp_function = Path('alt_allele_temp_function.sql').read_text()
    alt_allele_positions = alt_allele_positions_sql(Path('alt_allele_positions.sql').read_text(), include_vrs_allele_id)
    fq_vet_table = f"{fq_dataset}.{vet_table_name}"
    fq_sample_info = f"{fq_dataset}.sample_info"

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
    return query_return['results']


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
    parser.add_argument('--include_vrs_allele_id', action='store_true',
                        help='Populate alt_allele.vrs_allele_id from vet.vrs_allele_ids. Set only for a dataset whose vet '
                             'tables have vrs_allele_ids, in which case alt_allele has the column.')

    args = parser.parse_args()

    populate_alt_allele_table(args.call_set_identifier,
                              args.query_project,
                              args.vet_table_name,
                              args.fq_dataset,
                              args.max_sample_id,
                              args.include_vrs_allele_id)
