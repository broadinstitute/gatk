# -*- coding: utf-8 -*-
# Validates that the MANE join populated the VAT with the coverage it should (VS-1970).
#
# MANE is sparse (~19k transcripts, one MANE Select per gene), so a plain "fraction of rows populated"
# check is meaningless -- most rows legitimately have no MANE flag. Instead this compares, per
# MANE_status, the transcripts that SHOULD be flagged (MANE transcripts whose version-stripped base is
# present in the VAT) against those that ARE flagged (a VAT row with that base has the mane_*_name
# column set). A correct version-insensitive join flags ~100% of the present transcripts; a regression
# to exact-version matching (the VS-1970 bug) flags only ~6% and fails this check. The check re-derives
# the expected set with its own version-stripping, so a matching bug on the WDL side also shows up as a
# divergence. It cannot catch an error in the shared premise that version-insensitive matching is
# correct -- but that premise is well established (base ENST accessions are stable).
import argparse
import sys

from google.cloud import bigquery
from google.cloud.bigquery.job import QueryJobConfig

import utils


def write_output_files(passed, message, pass_file_output, results_file_output):
    with open(pass_file_output, 'w') as pass_file, open(results_file_output, 'w') as results_file:
        pass_file.write('true' if passed else 'false')
        results_file.write(message)


def check_mane_coverage(fq_vat_table, fq_mane_table, query_project, min_coverage,
                        pass_file_output, results_file_output):
    query_labels_map = {
        "id": "check_mane_coverage",
        "gvs_tool_name": "gvs_validate_mane"
    }
    # add labels for DSP Cloud Cost Control Labeling and Reporting
    query_labels_map.update({'service': 'gvs', 'team': 'variants', 'managedby': 'create_vat'})

    default_config = QueryJobConfig(labels=query_labels_map, priority="INTERACTIVE",
                                    use_query_cache=True, use_legacy_sql=False)
    client = bigquery.Client(project=query_project, default_query_job_config=default_config)

    # Per MANE_status: how many MANE transcripts are present in the callset (base match), and how many
    # of those are actually flagged. The `vat` CTE collapses to one row per base transcript, recording
    # whether any VAT row with that base carries each mane_*_name.
    sql = f"""
        WITH vat AS (
            SELECT
                SPLIT(transcript, '.')[OFFSET(0)] AS base,
                LOGICAL_OR(mane_select_name IS NOT NULL)        AS has_select,
                LOGICAL_OR(mane_plus_clinical_name IS NOT NULL) AS has_plus_clinical
            FROM `{fq_vat_table}`
            WHERE transcript IS NOT NULL
            GROUP BY base
        ),
        mane AS (
            SELECT DISTINCT MANE_status, SPLIT(Ensembl_nuc, '.')[OFFSET(0)] AS base
            FROM `{fq_mane_table}`
            WHERE Ensembl_nuc IS NOT NULL AND Ensembl_nuc != '-'
              AND MANE_status IN ('MANE Select', 'MANE Plus Clinical')
        )
        SELECT
            m.MANE_status AS mane_status,
            COUNTIF(v.base IS NOT NULL) AS present_in_callset,
            COUNTIF(v.base IS NOT NULL AND (
                (m.MANE_status = 'MANE Select'        AND v.has_select) OR
                (m.MANE_status = 'MANE Plus Clinical' AND v.has_plus_clinical)
            )) AS flagged
        FROM mane m
        LEFT JOIN vat v ON v.base = m.base
        GROUP BY m.MANE_status
        ORDER BY m.MANE_status
    """

    rows = list(utils.execute_with_retry(client, "validate mane coverage", sql).get('results'))

    lines = []
    failures = []
    for row in rows:
        status = row.get('mane_status')
        present = row.get('present_in_callset')
        flagged = row.get('flagged')
        if present == 0:
            lines.append(f"{status}: 0 present in callset (skipped)")
            continue
        coverage = flagged / present
        lines.append(f"{status}: {flagged}/{present} flagged ({coverage:.1%})")
        if coverage < min_coverage:
            failures.append(f"{status} coverage {coverage:.1%} is below the {min_coverage:.0%} floor "
                            f"({flagged} of {present} present transcripts flagged); the MANE join may "
                            f"have regressed to exact-version matching")

    if not rows:
        failures.append("no MANE rows found; the mane_annotations table may be missing or empty")

    passed = not failures
    prefix = "MANE coverage passed" if passed else "MANE coverage FAILED"
    message = f"{prefix}. " + "; ".join(lines)
    if failures:
        message += " || " + " || ".join(failures)

    print(message)
    write_output_files(passed, message, pass_file_output, results_file_output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False,
                                     description='Validate MANE annotation coverage in the VAT')
    parser.add_argument('--fq_vat_table', type=str, required=True,
                        help='fully-qualified VAT table (project.dataset.table)')
    parser.add_argument('--fq_mane_table', type=str, required=True,
                        help='fully-qualified MANE annotations table (project.dataset.mane_annotations)')
    parser.add_argument('--query_project', type=str, required=True,
                        help='Google project to run the BigQuery queries with')
    parser.add_argument('--min_coverage', type=float, required=False, default=0.95,
                        help='minimum fraction of present MANE transcripts that must be flagged')
    parser.add_argument('--pass_file_output', type=str, required=True,
                        help='path to write the pass/fail boolean')
    parser.add_argument('--results_file_output', type=str, required=True,
                        help='path to write the human-readable result')
    args = parser.parse_args()

    check_mane_coverage(args.fq_vat_table,
                        args.fq_mane_table,
                        args.query_project,
                        args.min_coverage,
                        args.pass_file_output,
                        args.results_file_output)
