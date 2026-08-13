import argparse

# add labels for DSP Cloud Cost Control Labeling and Reporting
query_labels_map = {'service': 'gvs', 'team': 'variants', 'managedby': 'gvs_process_sample_vcf_headers'}


# --- SQL builders (pure functions, no BigQuery client needed -- unit-testable) ---------------------

def vcf_header_lines_insert_sql(project_id, dataset_name):
    """Anti-join INSERT that populates vcf_header_lines from the scratch table.

    Deduplicating and idempotent: the Parquet ingest path writes the full header text for every chunk
    of every sample (no write-time dedup), and re-running (task retry, workflow resume, or re-ingest)
    inserts nothing already present -- the LEFT JOIN ... IS NULL skips hashes already in the target.
    See the VS-1968 design doc. This also fixes a latent double-insert bug on the legacy BQ path.

    The source subquery collapses to one row per hash; all texts for a given hash are identical by
    construction (the hash is an MD5 of the text), so ANY_VALUE is safe.
    """
    return f"""
        INSERT INTO `{project_id}.{dataset_name}.vcf_header_lines` (vcf_header_lines_hash, vcf_header_lines, is_expected_unique)
        SELECT
            s.vcf_header_lines_hash,
            ANY_VALUE(s.vcf_header_lines),
            ANY_VALUE(s.is_expected_unique)
        FROM `{project_id}.{dataset_name}.vcf_header_lines_scratch` s
        LEFT JOIN `{project_id}.{dataset_name}.vcf_header_lines` t USING (vcf_header_lines_hash)
        WHERE s.vcf_header_lines IS NOT NULL AND t.vcf_header_lines_hash IS NULL
        GROUP BY s.vcf_header_lines_hash
    """


def sample_vcf_header_insert_sql(project_id, dataset_name):
    """Anti-join INSERT for the sample->hash association map.

    Idempotent: only inserts (sample_id, hash) pairs not already present.
    """
    return f"""
        INSERT INTO `{project_id}.{dataset_name}.sample_vcf_header` (sample_id, vcf_header_lines_hash)
        SELECT DISTINCT s.sample_id, s.vcf_header_lines_hash
        FROM `{project_id}.{dataset_name}.vcf_header_lines_scratch` s
        LEFT JOIN `{project_id}.{dataset_name}.sample_vcf_header` t
            USING (sample_id, vcf_header_lines_hash)
        WHERE t.sample_id IS NULL
    """


def clean_up_scratch_sql(project_id, dataset_name):
    """DELETE the migrated rows from the scratch table (naturally idempotent -- re-delete is a no-op)."""
    return f"DELETE FROM `{project_id}.{dataset_name}.vcf_header_lines_scratch` WHERE vcf_header_lines_hash IS NOT NULL"


# --- Execution (imports the BigQuery client lazily so the SQL builders stay dependency-free) --------

def _make_client(project_id):
    from google.cloud import bigquery
    from google.cloud.bigquery.job import QueryJobConfig
    default_config = QueryJobConfig(labels=query_labels_map, priority="INTERACTIVE",
                                    use_query_cache=True, use_legacy_sql=False)
    return bigquery.Client(project=project_id, default_query_job_config=default_config)


def process_sample_vcf_headers(project_id, dataset_name, client=None):
    import utils
    if client is None:
        client = _make_client(project_id)
    populate_tables_from_scratch(project_id, dataset_name, client)
    clean_up_scratch_table(project_id, dataset_name, client)


def populate_tables_from_scratch(project_id, dataset_name, client):
    import utils
    utils.execute_with_retry(client, "vcf_header_lines", vcf_header_lines_insert_sql(project_id, dataset_name))
    utils.execute_with_retry(client, "sample_vcf_header", sample_vcf_header_insert_sql(project_id, dataset_name))


def clean_up_scratch_table(project_id, dataset_name, client):
    import utils
    utils.execute_with_retry(client, "vcf_header_lines_scratch", clean_up_scratch_sql(project_id, dataset_name))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False, description='')
    parser.add_argument('--project_id', type=str, help='Google project for the GVS dataset', required=True)
    parser.add_argument('--dataset_name',type=str, help='BigQuery dataset name', required=True)

    args = parser.parse_args()

    process_sample_vcf_headers(args.project_id, args.dataset_name)
