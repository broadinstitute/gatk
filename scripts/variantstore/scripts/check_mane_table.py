# -*- coding: utf-8 -*-
# Validates the loaded MANE (MANE.GRCh38.summary) BigQuery table before it is used to annotate the VAT.
#
# The MANE -> VAT join keys on the transcript (Ensembl_nuc, version-insensitive) and filters on
# MANE_status, so this checks a row-count floor (which catches an empty or partially loaded table --
# it would otherwise pass "table exists" and silently annotate nothing) plus the two columns the join
# depends on: Ensembl_nuc must be a populated Ensembl transcript id, and MANE_status must be one of the
# two expected values.
import argparse
import sys

from google.cloud import bigquery
from google.cloud.bigquery.job import QueryJobConfig

import utils


def check_mane_table(fq_mane_table, query_project, min_rows):
    query_labels_map = {
        "id": "check_mane_table",
        "gvs_tool_name": "gvs_validate_mane"
    }
    # add labels for DSP Cloud Cost Control Labeling and Reporting
    query_labels_map.update({'service': 'gvs', 'team': 'variants', 'managedby': 'create_vat'})

    default_config = QueryJobConfig(labels=query_labels_map, priority="INTERACTIVE", use_query_cache=True, use_legacy_sql=False)
    client = bigquery.Client(project=query_project,
                             default_query_job_config=default_config)

    sql = f"""
        SELECT
            COUNT(*) AS total_rows,
            COUNTIF(Ensembl_nuc IS NULL OR NOT STARTS_WITH(Ensembl_nuc, 'ENST')) AS bad_transcripts,
            COUNTIF(MANE_status NOT IN ('MANE Select', 'MANE Plus Clinical')) AS bad_statuses
        FROM `{fq_mane_table}`
    """

    query_return = utils.execute_with_retry(client, "validate mane table", sql)
    row = next(iter(query_return.get('results')))
    total_rows = row.get('total_rows')
    bad_transcripts = row.get('bad_transcripts')
    bad_statuses = row.get('bad_statuses')

    print(f"MANE table {fq_mane_table}: total_rows={total_rows}, bad_transcripts={bad_transcripts}, "
          f"bad_statuses={bad_statuses}")

    errors = []
    if total_rows < min_rows:
        errors.append(f"{total_rows} rows is below the expected floor of {min_rows} "
                      f"(table may be empty or only partially loaded)")
    if bad_transcripts > 0:
        errors.append(f"{bad_transcripts} rows have a NULL or non-Ensembl (not ENST...) Ensembl_nuc "
                      f"(the transcript join key)")
    if bad_statuses > 0:
        errors.append(f"{bad_statuses} rows have an unexpected MANE_status "
                      f"(expected 'MANE Select' or 'MANE Plus Clinical')")

    if errors:
        print("MANE table validation FAILED:", file=sys.stderr)
        for e in errors:
            print(f"  - {e}", file=sys.stderr)
        sys.exit(1)

    print(f"MANE table validation passed: {total_rows} rows, all Ensembl_nuc are ENST, all MANE_status valid.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False,
                                     description='Validate the loaded MANE (MANE.GRCh38.summary) BigQuery table')
    parser.add_argument('--fq_mane_table', type=str, required=True,
                        help='fully-qualified BigQuery MANE table (project.dataset.table)')
    parser.add_argument('--query_project', type=str, required=True,
                        help='Google project to run the BigQuery queries with')
    parser.add_argument('--min_rows', type=int, required=False, default=15000,
                        help='minimum expected row count; MANE 1.4 is ~19.4k rows')

    args = parser.parse_args()

    check_mane_table(args.fq_mane_table,
                     args.query_project,
                     args.min_rows)
