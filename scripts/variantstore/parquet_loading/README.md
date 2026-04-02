# Parquet to BigQuery Loading System

Robust, fault-tolerant system for loading hundreds of thousands of Parquet files from Google Cloud Storage into BigQuery tables.

## Overview

This implementation loads Parquet files into BigQuery with:
- **Idempotency**: Safe to run multiple times — already-loaded data is detected by inspecting BigQuery directly
- **Fault Tolerance**: Survives VM preemptions and failures
- **Deterministic Job IDs**: Can resume interrupted BigQuery load jobs
- **Verification**: Confirms all files were successfully loaded by querying BigQuery partition data

## Architecture

### Components

1. **`parse_and_group_files.py`**: Discovers Parquet files in GCS, groups by target table, filters already-loaded files by querying `INFORMATION_SCHEMA.PARTITIONS` and the `sample_chromosome_ploidy` table directly
2. **`load_parquet_to_bq.py`**: Loads Parquet files with deterministic job IDs and preemption recovery
3. **`verify_all_loaded.py`**: Verifies all files were loaded successfully by comparing GCS-derived `(table_name, sample_id)` pairs against BigQuery partition data

**Note**: All Python scripts are located in `scripts/variantstore/scripts/` and are packaged into the variants Docker image at `/app/`.

### Idempotency via BigQuery Partition Inspection

Idempotency is determined by querying BigQuery directly rather than maintaining a separate tracking table:

- **Superpartitioned tables** (`vet_%`, `ref_ranges_%`): `INFORMATION_SCHEMA.PARTITIONS` is queried for non-empty partitions whose partition ID parses as a valid `sample_id`.
- **Regular tables** (`sample_chromosome_ploidy`): The table itself is queried for distinct `sample_id` values.

If a `(table_name, sample_id)` pair is already represented in BigQuery, all Parquet files for that pair are skipped without attempting a reload. This approach is more reliable than a tracking table because it reflects the actual committed state of BigQuery.

### GCS Directory Structure

Individual Parquet files are named like

```
<data type>_<superpartition>_<sample id>_<VCF filename>.parquet
```

where `<data type>` is `vet` or `ref_ranges`.

Sample GCS directory structure:

```
gs://bucket/output_dir/
├── vet/
│   ├── 001/
│   │   ├── vet_001_1_sample1.parquet
│   │   ├── vet_001_2_sample1.parquet
...
│   └── 002/
│       ├── vet_002_4001_sample1.parquet
...
└── ref_ranges/
│   ├── 001/
│   │   ├── ref_ranges_001_2_sample1.parquet
...
```

**Mapping**: `vet/001/` → BigQuery table `vet_001`, `ref_ranges/042/` → BigQuery table `ref_ranges_042`

## Prerequisites

- Python 3.8+
- Google Cloud SDK
- BigQuery API enabled
- Service account with BigQuery permissions:
  - `bigquery.tables.create`
  - `bigquery.tables.get`
  - `bigquery.tables.getData`
  - `bigquery.tables.updateData`
  - `bigquery.jobs.create`
- GCS read permissions for source buckets

## Docker Image

The scripts are packaged in the **variants Docker image** built from `scripts/variantstore/scripts/Dockerfile`.

### Building the Image

```bash
cd scripts/variantstore/scripts/
./build_docker.sh
```

**Build Process:**
1. Creates Docker buildx builder for linux/amd64
2. Builds image with all Python scripts copied to `/app/`
3. Runs unit tests from `test/test_*.py`
4. Tags image: `YYYY-MM-DD-alpine-<image-id>`
5. Pushes to Google Artifact Registry: `us-central1-docker.pkg.dev/broad-dsde-methods/gvs/variants:<tag>`

**Prerequisites for Building:**
- Docker configured for GAR: `gcloud auth configure-docker us-central1-docker.pkg.dev`
- Access to `broad-dsde-methods` GCP project

### Using Pre-built Images

The WDL workflow uses the variants Docker image. Update the `python_bq_docker` input to use your built image:

```json
{
  "LoadParquetsToBigQuery.python_bq_docker": "us-central1-docker.pkg.dev/broad-dsde-methods/gvs/variants:2025-10-28-alpine-abc123"
}
```

## Installation

```bash
pip install google-cloud-bigquery google-cloud-storage
```

## Usage

### 1. Run Standalone Scripts (Local Testing)

**Note**: For production use, the scripts should be run from the Docker container. These examples show local execution for testing.

#### Discover and Group Files

```bash
# List all Parquet files
gcloud storage ls --recursive gs://my-bucket/output_dir/ > all_objects.txt
grep '\.parquet$' all_objects.txt > all_files.txt

# Group by table (already-loaded (table_name, sample_id) pairs are skipped automatically)
python3 scripts/variantstore/scripts/parse_and_group_files.py \
  --input all_files.txt \
  --output-dir grouped_files \
  --project-id my-project \
  --dataset my_dataset \
  --superpartitioned-table-prefixes vet ref_ranges \
  --regular-table-prefixes sample_chromosome_ploidy
```

#### Load Files for a Table

```bash
python3 scripts/variantstore/scripts/load_parquet_to_bq.py \
  --project-id my-project \
  --dataset-name my_dataset \
  --table-name vet_001 \
  --files-fofn grouped_files/vet_001.fofn \
  --schema-path schemas/vet_001_schema.json \
  --pending-jobs-path pending_jobs_vet_001.json \
  --batch-size 1000 \
  --output-stats stats_vet_001.json
```

#### Verify All Files Loaded

```bash
python3 scripts/variantstore/scripts/verify_all_loaded.py \
  --project-id my-project \
  --dataset-name my_dataset \
  --gcs-files-list all_files.txt \
  --output-dir verification_output
```

### 3. Run via Cromwell Workflow (Recommended)

Create an inputs JSON file:

```json
{
  "LoadParquetsToBigQuery.output_gcs_dir": "gs://my-bucket/output_dir",
  "LoadParquetsToBigQuery.project_id": "my-project",
  "LoadParquetsToBigQuery.dataset_name": "my_dataset",
  "LoadParquetsToBigQuery.table_prefixes": ["vet", "ref_ranges"],
  "LoadParquetsToBigQuery.batch_size": 10000,
  "LoadParquetsToBigQuery.python_bq_docker": "us-central1-docker.pkg.dev/broad-dsde-methods/gvs/variants:2025-10-28-alpine-abc123"
}
```

**Note**: Replace the `python_bq_docker` value with your actual built image tag from the build process.

Run the workflow:

```bash
cromwell run scripts/variantstore/wdl/LoadParquetsToBigQuery.wdl \
  -i inputs.json
```

## Configuration

### Batch Size

- **Default**: 10,000 files per BigQuery load job
- **Recommend**: 500-1,000 for faster retries on large tables
- **Max**: 10,000 (BigQuery limit)

### Preemptible VMs

All tasks use preemptible VMs (3-5 attempts) to save ~70% on compute costs. The deterministic job ID system ensures safe resumption after preemptions.

### Parallelization

The workflow loads different tables in parallel. Control concurrency with `max_parallel_loads` (default: 50) to stay within BigQuery's ~100 concurrent load job limit per project.

## Fault Tolerance

### Preemption Recovery

When a VM is preempted during loading:

1. **Job IDs are deterministic**: Hash of table name + batch contents
2. **State persisted**: `pending_jobs.json` tracks submitted jobs
3. **On restart**: Call `client.get_job(job_id)` to check status
4. **Resume or skip**: If job completed, record success; if failed/missing, resubmit

### Idempotency on Workflow Restart

On workflow restart, `parse_and_group_files.py` queries BigQuery (`INFORMATION_SCHEMA.PARTITIONS` for `vet_%`/`ref_ranges_%` tables and `sample_chromosome_ploidy` directly) to determine which `(table_name, sample_id)` pairs are already loaded. Files belonging to already-loaded pairs are skipped. No separate tracking table is consulted or required.

## Monitoring

### Check Load Status

```sql
-- Overall progress per table (via partition metadata)
SELECT
  table_name,
  COUNT(*) as partitions_loaded,
  SUM(total_logical_bytes) as total_bytes
FROM `my-project.my_dataset.INFORMATION_SCHEMA.PARTITIONS`
WHERE REGEXP_CONTAINS(table_name, "^vet_[0-9]+$|^ref_ranges_[0-9]+$")
  AND total_logical_bytes > 0
  AND NOT STARTS_WITH(partition_id, '__')
GROUP BY table_name
ORDER BY table_name;

-- Ploidy samples loaded
SELECT COUNT(DISTINCT sample_id) as samples_with_ploidy
FROM `my-project.my_dataset.sample_chromosome_ploidy`;
```

## Troubleshooting

### No files discovered

- Check GCS path and permissions
- Verify files match pattern `*.parquet`
- Confirm table prefixes are correct

### Load job failures

- Check BigQuery quotas (concurrent jobs per project)
- Verify schema compatibility
- Review stderr logs for detailed errors

### Verification fails

- Re-run the workflow (idempotent)
- Check `missing_files.txt` for paths
- Investigate load failures by querying `INFORMATION_SCHEMA.PARTITIONS`

### Quota exceeded errors

- Reduce `max_parallel_loads`
- Reduce `batch_size` to spread load over time
- Request quota increase from Google Cloud

## Performance

### Typical Scale

- **100,000 files across 25 tables**:
  - Discovery: ~1 minute
  - Loading: ~10 minutes (parallel)
  - Verification: ~1 minute
  - Cost: ~$0.50 (with preemptibles) + BigQuery storage

### Optimization

- **Parallel loading**: Different tables load concurrently
- **Batching**: 500-10,000 files per BigQuery job
- **Preemptible VMs**: 70% cost savings
- **Free BigQuery loads**: Loads from GCS have no API charges

## Design Document

See `docs/parquet_to_bigquery_loading_design.md` for detailed architecture, edge cases, and implementation notes.

## Known Limitations

- BigQuery concurrent load jobs: ~100 per project/region (soft limit)
- Maximum 10,000 URIs per load job
- Tables and schemas must be pre-created
