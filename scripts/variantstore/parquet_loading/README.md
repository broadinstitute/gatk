# Parquet to BigQuery Loading System

Robust, fault-tolerant system for loading hundreds of thousands of Parquet files from Google Cloud Storage into BigQuery tables.

## Overview

This implementation loads Parquet files into BigQuery with:
- **Idempotency**: Safe to run multiple times
- **Fault Tolerance**: Survives VM preemptions and failures
- **Deterministic Job IDs**: Can resume interrupted BigQuery load jobs
- **Tracking**: Records which files have been loaded in a BigQuery table
- **Verification**: Confirms all files were successfully loaded

## Architecture

### Components

1. **`create_tracking_table.py`**: Creates the BigQuery tracking table
2. **`parse_and_group_files.py`**: Discovers Parquet files in GCS, groups by target table, filters already-loaded files
3. **`load_parquet_to_bq.py`**: Loads Parquet files with deterministic job IDs and preemption recovery
4. **`verify_all_loaded.py`**: Verifies all files were loaded successfully
5. **`LoadParquetsToBigQuery.wdl`**: Cromwell workflow orchestrating all tasks

**Note**: All Python scripts are located in `scripts/variantstore/scripts/` and are packaged into the variants Docker image at `/app/`.

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
│   │   ├── vet_001_2_sample1.parquet
...
```

**Mapping**: `vet/001/` → BigQuery table `vet_001`, `ref_ranges/042/` → table `ref_042`

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

### 1. Create Tracking Table

```bash
python3 scripts/variantstore/parquet_loading/create_tracking_table.py \
  --project-id my-project \
  --dataset-name my_dataset
```

### 2. Run Standalone Scripts (Local Testing)

**Note**: For production use, the scripts should be run from the Docker container. These examples show local execution for testing.

#### Discover and Group Files

```bash
# List all Parquet files
gcloud storage ls --recursive gs://my-bucket/output_dir/ > all_objects.txt
grep '\.parquet$' all_objects.txt > all_files.txt

# Group by table
python3 scripts/variantstore/scripts/parse_and_group_files.py \
  --input all_files.txt \
  --output-dir grouped_files \
  --project-id my-project \
  --dataset my_dataset \
  --table-prefixes vet ref_ranges
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

### Tracking Table

The `parquet_load_status` table records:
- `file_path` (PRIMARY KEY)
- `table_name`
- `load_timestamp`
- `load_job_id`
- `file_size_bytes`
- `rows_loaded`

Files are only inserted after successful load, preventing duplicate loads.

## Monitoring

### Check Load Status

```sql
-- Overall progress
SELECT 
  table_name,
  COUNT(*) as files_loaded,
  SUM(rows_loaded) as total_rows,
  MIN(load_timestamp) as first_load,
  MAX(load_timestamp) as last_load
FROM `my-project.my_dataset.parquet_load_status`
GROUP BY table_name
ORDER BY table_name;

-- Recent loads
SELECT *
FROM `my-project.my_dataset.parquet_load_status`
WHERE load_timestamp > TIMESTAMP_SUB(CURRENT_TIMESTAMP(), INTERVAL 1 HOUR)
ORDER BY load_timestamp DESC;

-- Check for duplicates (should be 0)
SELECT file_path, COUNT(*) as count
FROM `my-project.my_dataset.parquet_load_status`
GROUP BY file_path
HAVING COUNT(*) > 1;
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
- Investigate load failures in tracking table

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
- Tracking table must exist before loading

## Future Enhancements

- Pub/Sub triggers for incremental loading
- Schema validation before loading
- Automatic retry with exponential backoff
- BigQuery table partitioning support
- Cost tracking and reporting
