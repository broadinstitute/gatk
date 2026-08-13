# -*- coding: utf-8 -*-
# Validates the loaded Entrez (gene2ensembl) BigQuery table before it is used to annotate the VAT.
#
# The Entrez -> VAT join keys on the gene (Ensembl_gene_identifier -> GeneID), so this checks the
# two columns that join depends on, plus a row-count floor that catches an empty or partially loaded
# table (which would otherwise pass "table exists" + "0 bad rows" and silently annotate nothing).
#
# Note: we intentionally do NOT check Ensembl_rna_identifier for NULLs. gene2ensembl legitimately
# leaves the transcript columns blank ("-", loaded as NULL) for a large fraction of rows, and the
# gene-level join does not use them.
import argparse
import sys

from google.cloud import bigquery
from google.cloud.bigquery.job import QueryJobConfig

import utils


def check_entrez_table(fq_entrez_table, query_project, min_rows):
    query_labels_map = {
        "id": "check_entrez_table",
        "gvs_tool_name": "gvs_validate_entrez"
    }

    # add labels for DSP Cloud Cost Control Labeling and Reporting
    query_labels_map.update({'service': 'gvs', 'team': 'variants', 'managedby': 'create_vat'})

    default_config = QueryJobConfig(labels=query_labels_map, priority="INTERACTIVE", use_query_cache=True, use_legacy_sql=False)
    client = bigquery.Client(project=query_project,
                             default_query_job_config=default_config)

    sql = f"""
        SELECT
            COUNT(*) AS total_rows,
            COUNTIF(SAFE_CAST(GeneID AS INT64) IS NULL) AS bad_gene_ids,
            COUNTIF(Ensembl_gene_identifier IS NULL OR Ensembl_gene_identifier = '-') AS bad_gene_identifiers
        FROM `{fq_entrez_table}`
    """

    query_return = utils.execute_with_retry(client, "validate entrez table", sql)
    row = next(iter(query_return.get('results')))
    total_rows = row.get('total_rows')
    bad_gene_ids = row.get('bad_gene_ids')
    bad_gene_identifiers = row.get('bad_gene_identifiers')

    print(f"Entrez table {fq_entrez_table}: total_rows={total_rows}, bad_gene_ids={bad_gene_ids}, "
          f"bad_gene_identifiers={bad_gene_identifiers}")

    errors = []
    if total_rows < min_rows:
        errors.append(f"{total_rows} rows is below the expected floor of {min_rows} "
                      f"(table may be empty or only partially loaded)")
    if bad_gene_ids > 0:
        errors.append(f"{bad_gene_ids} rows have a NULL or non-numeric GeneID")
    if bad_gene_identifiers > 0:
        errors.append(f"{bad_gene_identifiers} rows have a NULL or placeholder ('-') Ensembl_gene_identifier "
                      f"(the gene-level join key)")

    if errors:
        print("Entrez table validation FAILED:", file=sys.stderr)
        for e in errors:
            print(f"  - {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Entrez table validation passed: all GeneIDs numeric, all Ensembl_gene_identifiers populated, "
          f"{total_rows} rows.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False,
                                     description='Validate the loaded Entrez (gene2ensembl) BigQuery table')
    parser.add_argument('--fq_entrez_table', type=str, required=True,
                        help='fully-qualified BigQuery Entrez table (project.dataset.table)')
    parser.add_argument('--query_project', type=str, required=True,
                        help='Google project to run the BigQuery queries with')
    parser.add_argument('--min_rows', type=int, required=False, default=50000,
                        help='minimum expected row count; gene2ensembl_human is ~88k rows')

    args = parser.parse_args()

    check_entrez_table(args.fq_entrez_table,
                       args.query_project,
                       args.min_rows)
