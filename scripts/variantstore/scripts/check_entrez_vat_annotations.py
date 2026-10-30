# -*- coding: utf-8 -*-
# Post-join validation of the VAT's entrez_gene_id annotation. Runs after the VAT is built.
#
# Two checks:
#   1. Coverage  - of rows that carry a gene_id, at least `min_coverage` fraction receive a non-empty
#                  entrez_gene_id array. entrez_gene_id is a repeated field, so an unmatched row is an
#                  empty array [] (ARRAY_LENGTH = 0), not NULL. This guards against a regression back to
#                  the sparse (~26%) transcript-keyed join, or to 0% if the join silently matches nothing.
#   2. Consistency - entrez_gene_id is a gene-level attribute, so every row sharing a gene_id must carry
#                  the same set of GeneIDs. Fails if any gene_id maps to more than one distinct array.
import argparse

from google.cloud import bigquery
from google.cloud.bigquery.job import QueryJobConfig

import utils


def write_output_files(passed, message, pass_file_output, results_file_output):
    with open(pass_file_output, 'w') as pass_file, open(results_file_output, 'w') as results_file:
        pass_file.write('true' if passed else 'false')
        results_file.write(message)


def check_entrez_vat_annotations(fq_vat_table, query_project, min_coverage, pass_file_output, results_file_output):
    query_labels_map = {
        "id": "check_entrez_vat_annotations",
        "gvs_tool_name": "gvs_validate_vat"
    }

    # add labels for DSP Cloud Cost Control Labeling and Reporting
    query_labels_map.update({'service': 'gvs', 'team': 'variants', 'managedby': 'create_vat'})

    default_config = QueryJobConfig(labels=query_labels_map, priority="INTERACTIVE", use_query_cache=True, use_legacy_sql=False)
    client = bigquery.Client(project=query_project,
                             default_query_job_config=default_config)

    # Coverage: fraction of gene-bearing rows that received at least one Entrez GeneID.
    coverage_sql = f"""
        SELECT
            COUNTIF(gene_id IS NOT NULL) AS gene_rows,
            COUNTIF(gene_id IS NOT NULL AND ARRAY_LENGTH(entrez_gene_id) > 0) AS annotated_rows
        FROM `{fq_vat_table}`
    """
    coverage_row = next(iter(utils.execute_with_retry(client, "entrez coverage", coverage_sql).get('results')))
    gene_rows = coverage_row.get('gene_rows')
    annotated_rows = coverage_row.get('annotated_rows')
    coverage = (annotated_rows / gene_rows) if gene_rows else 0.0

    # Consistency: each gene_id must map to a single distinct entrez_gene_id array. The array is
    # collapsed to an order-independent string (sorted, with empty arrays becoming '') so that
    # arrays with the same members in a different order are treated as equal, while an empty array
    # is still distinguished from a populated one.
    consistency_sql = f"""
        WITH per_gene AS (
            SELECT
                gene_id,
                COUNT(DISTINCT IFNULL(
                    (SELECT STRING_AGG(CAST(e AS STRING), ',' ORDER BY e) FROM UNNEST(entrez_gene_id) AS e),
                    '')) AS distinct_arrays
            FROM `{fq_vat_table}`
            WHERE gene_id IS NOT NULL
            GROUP BY gene_id
        )
        SELECT COUNTIF(distinct_arrays > 1) AS inconsistent_genes, COUNT(*) AS total_genes FROM per_gene
    """
    consistency_row = next(iter(utils.execute_with_retry(client, "entrez consistency", consistency_sql).get('results')))
    inconsistent_genes = consistency_row.get('inconsistent_genes')
    total_genes = consistency_row.get('total_genes')

    coverage_ok = coverage >= min_coverage
    consistency_ok = inconsistent_genes == 0
    passed = coverage_ok and consistency_ok

    coverage_msg = (f"entrez_gene_id coverage {coverage:.4f} ({annotated_rows}/{gene_rows} gene-bearing rows) "
                    f"{'>=' if coverage_ok else '<'} threshold {min_coverage}")
    consistency_msg = (f"{inconsistent_genes} of {total_genes} gene_ids map to more than one distinct "
                       f"entrez_gene_id array")
    message = f"{fq_vat_table}: {coverage_msg}; {consistency_msg}."

    write_output_files(passed, message, pass_file_output, results_file_output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False,
                                     description='Check VAT entrez_gene_id coverage and gene-consistency')
    parser.add_argument('--fq_vat_table', type=str, required=True,
                        help='fully-qualified BigQuery VAT table')
    parser.add_argument('--query_project', type=str, required=True,
                        help='Google project to run the BigQuery queries with')
    parser.add_argument('--min_coverage', type=float, required=False, default=0.85,
                        help='minimum fraction of gene-bearing rows that must carry an Entrez GeneID')
    parser.add_argument('--pass_file_output', type=str, required=True,
                        help='location to write file with pass/fail into')
    parser.add_argument('--results_file_output', type=str, required=True,
                        help='location to write query results into')

    args = parser.parse_args()

    check_entrez_vat_annotations(args.fq_vat_table,
                                 args.query_project,
                                 args.min_coverage,
                                 args.pass_file_output,
                                 args.results_file_output)
