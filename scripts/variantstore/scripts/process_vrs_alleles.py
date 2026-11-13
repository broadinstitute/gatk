import argparse

# add labels for DSP Cloud Cost Control Labeling and Reporting
query_labels_map = {'service': 'gvs', 'team': 'variants', 'managedby': 'gvs_process_vrs_alleles'}

SCRATCH_TABLE = "vrs_allele_scratch"
TARGET_TABLE = "vrs_allele"

# Every vrs_allele column. The scratch table has these plus sample_id.
VRS_ALLELE_COLUMNS = [
    "vrs_allele_id",
    "vrs_location_id",
    "refget_accession",
    "ref_genome_coord_start",
    "ref_genome_coord_end",
    "state_type",
    "state_sequence",
    "state_length",
    "state_repeat_subunit_length",
]

# What ingest writes to vet.vrs_allele_ids for a spanning-deletion `*` ALT (VetCreator.NO_VRS_ALLELE_ID).
NO_VRS_ALLELE_ID = ""

# How many offending ids to show in a failure message.
EXAMPLE_LIMIT = 10


# --- SQL builders (pure functions, no BigQuery client needed -- unit-testable) ---------------------

def _sample_id_list(sample_ids):
    """Render sample ids as a SQL literal list. Literals, not a subquery, so BigQuery prunes the
    sample_id-partitioned vet tables to these samples instead of scanning every partition."""
    if not sample_ids:
        raise ValueError("sample_ids must not be empty")
    return ", ".join(str(int(sid)) for sid in sorted(sample_ids))


def scratch_sample_ids_sql(project_id, dataset_name):
    """The samples whose VRS alleles are waiting in the scratch table: this run's, scoping every
    later statement so a concurrent ingest's rows are left alone."""
    return f"""
        SELECT DISTINCT sample_id
        FROM `{project_id}.{dataset_name}.{SCRATCH_TABLE}`
    """


def conflicting_ids_sql(project_id, dataset_name, sample_ids):
    """vrs_allele_ids whose rows disagree on any other column, within the scratch rows or against
    the rows already in vrs_allele.

    The id is a content hash of the location and state, so a conflict means two different alleles
    hashed to one id -- a bug in the ID computation, not in the data. Merging would keep an arbitrary
    one of them, so this must find nothing before anything is merged.
    """
    columns = ", ".join(VRS_ALLELE_COLUMNS)
    return f"""
        WITH candidate AS (
            SELECT {columns}
            FROM `{project_id}.{dataset_name}.{SCRATCH_TABLE}`
            WHERE sample_id IN ({_sample_id_list(sample_ids)})
            UNION ALL
            SELECT {columns}
            FROM `{project_id}.{dataset_name}.{TARGET_TABLE}`
        )
        SELECT vrs_allele_id
        FROM candidate
        GROUP BY vrs_allele_id
        HAVING COUNT(DISTINCT TO_JSON_STRING(STRUCT({columns}))) > 1
        ORDER BY vrs_allele_id
        LIMIT {EXAMPLE_LIMIT}
    """


def vrs_allele_insert_sql(project_id, dataset_name, sample_ids):
    """Anti-join INSERT that adds this run's distinct alleles to vrs_allele.

    Deduplicating and idempotent: the same allele arrives once per sample that carries it, and
    re-running (task retry, workflow resume, or re-ingest) inserts nothing already present. All rows
    for one id are identical once conflicting_ids_sql has found nothing, so ANY_VALUE is safe.
    """
    columns = ", ".join(VRS_ALLELE_COLUMNS)
    any_values = ",\n            ".join(f"ANY_VALUE(s.{c})" for c in VRS_ALLELE_COLUMNS[1:])
    return f"""
        INSERT INTO `{project_id}.{dataset_name}.{TARGET_TABLE}` ({columns})
        SELECT
            s.vrs_allele_id,
            {any_values}
        FROM `{project_id}.{dataset_name}.{SCRATCH_TABLE}` s
        LEFT JOIN `{project_id}.{dataset_name}.{TARGET_TABLE}` t USING (vrs_allele_id)
        WHERE s.sample_id IN ({_sample_id_list(sample_ids)}) AND t.vrs_allele_id IS NULL
        GROUP BY s.vrs_allele_id
    """


def missing_vet_ids_sql(project_id, dataset_name, sample_ids):
    """vrs_allele_ids that this run's vet rows reference but vrs_allele does not hold.

    This is the check that ties the two halves together: vet.vrs_allele_ids is the per-sample
    association, vrs_allele the lookup, and a vet id with no lookup row is a dangling reference.
    The table wildcard is restricted to the numbered vet tables. A spanning-deletion `*` ALT has the
    NO_VRS_ALLELE_ID placeholder instead of an id, so it is not a dangling reference.
    """
    return f"""
        SELECT v.sample_id, vet_id AS vrs_allele_id
        FROM `{project_id}.{dataset_name}.vet_*` v
        CROSS JOIN UNNEST(v.vrs_allele_ids) AS vet_id
        LEFT JOIN `{project_id}.{dataset_name}.{TARGET_TABLE}` a ON a.vrs_allele_id = vet_id
        WHERE REGEXP_CONTAINS(_TABLE_SUFFIX, r'^[0-9]+$')
          AND v.sample_id IN ({_sample_id_list(sample_ids)})
          AND vet_id != '{NO_VRS_ALLELE_ID}'
          AND a.vrs_allele_id IS NULL
        ORDER BY v.sample_id, vrs_allele_id
        LIMIT {EXAMPLE_LIMIT}
    """


def duplicate_ids_sql(project_id, dataset_name):
    """vrs_allele_ids appearing more than once in vrs_allele. The table is meant to be keyed on the
    id; concurrent merges into one dataset are the realistic way to break that."""
    return f"""
        SELECT vrs_allele_id
        FROM `{project_id}.{dataset_name}.{TARGET_TABLE}`
        GROUP BY vrs_allele_id
        HAVING COUNT(*) > 1
        ORDER BY vrs_allele_id
        LIMIT {EXAMPLE_LIMIT}
    """


def clean_up_scratch_sql(project_id, dataset_name, sample_ids):
    """DELETE this run's rows from the scratch table (naturally idempotent -- re-delete is a no-op)."""
    return (f"DELETE FROM `{project_id}.{dataset_name}.{SCRATCH_TABLE}` "
            f"WHERE sample_id IN ({_sample_id_list(sample_ids)})")


# --- Execution (imports the BigQuery client lazily so the SQL builders stay dependency-free) --------

def _make_client(project_id):
    from google.cloud import bigquery
    from google.cloud.bigquery.job import QueryJobConfig
    default_config = QueryJobConfig(labels=query_labels_map, priority="INTERACTIVE",
                                    use_query_cache=False, use_legacy_sql=False)
    return bigquery.Client(project=project_id, default_query_job_config=default_config)


def _execute(client, label, sql):
    import utils
    return utils.execute_with_retry(client, label, sql)


def _rows(client, label, sql):
    return list(_execute(client, label, sql)['results'])


def process_vrs_alleles(project_id, dataset_name, client=None):
    """Merge this run's VRS alleles from the scratch table into vrs_allele, check the result, and
    empty the scratch rows.

    Each check raises before the scratch rows are deleted, so a failure leaves them in place for
    inspection, and a re-run after a fix picks up where this one stopped: the merge is idempotent.
    """
    if client is None:
        client = _make_client(project_id)

    sample_ids = [row.sample_id for row in
                  _rows(client, "vrs_allele_scratch samples", scratch_sample_ids_sql(project_id, dataset_name))]
    if not sample_ids:
        print(f"No rows in {SCRATCH_TABLE}; nothing to merge.")
        return
    print(f"Merging VRS alleles for {len(sample_ids)} sample(s) from {SCRATCH_TABLE} into {TARGET_TABLE}.")

    conflicting = _rows(client, "vrs_allele conflicting ids",
                        conflicting_ids_sql(project_id, dataset_name, sample_ids))
    if conflicting:
        raise RuntimeError(
            f"{TARGET_TABLE} merge aborted: these vrs_allele_ids have rows that disagree on location or "
            f"state, which means two different alleles got the same ID (showing up to {EXAMPLE_LIMIT}): "
            f"{[row.vrs_allele_id for row in conflicting]}. Nothing was merged; {SCRATCH_TABLE} is unchanged.")

    _execute(client, TARGET_TABLE, vrs_allele_insert_sql(project_id, dataset_name, sample_ids))

    missing = _rows(client, "vet vrs_allele_ids missing from vrs_allele",
                    missing_vet_ids_sql(project_id, dataset_name, sample_ids))
    if missing:
        raise RuntimeError(
            f"vet rows reference vrs_allele_ids that are not in {TARGET_TABLE} after the merge (showing up to "
            f"{EXAMPLE_LIMIT} (sample_id, vrs_allele_id) pairs): "
            f"{[(row.sample_id, row.vrs_allele_id) for row in missing]}. {SCRATCH_TABLE} is unchanged.")

    duplicates = _rows(client, "vrs_allele duplicate ids", duplicate_ids_sql(project_id, dataset_name))
    if duplicates:
        raise RuntimeError(
            f"{TARGET_TABLE} holds more than one row for these vrs_allele_ids (showing up to {EXAMPLE_LIMIT}): "
            f"{[row.vrs_allele_id for row in duplicates]}. Was another ingest merging into this dataset at "
            f"the same time? {SCRATCH_TABLE} is unchanged.")

    _execute(client, SCRATCH_TABLE, clean_up_scratch_sql(project_id, dataset_name, sample_ids))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False, description='')
    parser.add_argument('--project_id', type=str, help='Google project for the GVS dataset', required=True)
    parser.add_argument('--dataset_name', type=str, help='BigQuery dataset name', required=True)

    args = parser.parse_args()

    process_vrs_alleles(args.project_id, args.dataset_name)
