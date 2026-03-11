# BigQuery Quota Handling Implementation

## Overview

This document describes the quota error handling and retry logic implemented in the Parquet-to-BigQuery loading system.

## Changes Made

### 1. Import Retry Logic (✅ Implemented)

Added exponential backoff retry logic directly into `load_parquet_to_bq.py` for handling BigQuery quota and rate limit errors.

### 2. Catch Specific Exceptions (✅ Implemented)

The code now specifically catches and handles these retryable BigQuery errors:
- `google.api_core.exceptions.TooManyRequests` - Quota exceeded errors
- `google.api_core.exceptions.ServiceUnavailable` - Temporary service unavailability
- `google.api_core.exceptions.InternalServerError` - Transient internal errors

### 3. Add Exponential Backoff (✅ Implemented)

Implemented three-tier exponential backoff strategy:
- **Attempt 1**: Wait 30 seconds before retry
- **Attempt 2**: Wait 60 seconds before retry
- **Attempt 3**: Wait 120 seconds before retry
- **After 3 retries**: Fail the batch and continue with next batch

## Implementation Details

### New Helper Functions

#### `_submit_load_job_with_retry()`
Handles quota errors when submitting BigQuery load jobs:
```python
def _submit_load_job_with_retry(client, batch, table_id, job_config, job_id, location, max_retries=3)
```

**Features:**
- Retries on `TooManyRequests`, `ServiceUnavailable`, `InternalServerError`
- Handles `Conflict` exception (job already exists) by fetching existing job
- Fails immediately on non-retryable errors
- Logs retry attempts with delay information

#### `_wait_for_job_with_retry()`
Handles quota errors while waiting for job completion:
```python
def _wait_for_job_with_retry(load_job, max_retries=3)
```

**Features:**
- Retries transient errors during job polling
- Reloads job state after each retry
- Distinguishes between transient and actual job failures

#### `_execute_query_with_retry()`
Handles quota errors for tracking table queries:
```python
def _execute_query_with_retry(client, query, job_config=None, max_retries=3)
```

**Features:**
- Used for inserting tracking records
- Same retry strategy as load jobs
- Prevents tracking table updates from failing due to quota

### Retry Strategy

**Exponential Backoff Schedule:**
| Attempt | Delay | Cumulative Wait |
|---------|-------|-----------------|
| 1       | 0s    | 0s              |
| 2       | 30s   | 30s             |
| 3       | 60s   | 90s             |
| 4       | 120s  | 210s            |

**Total maximum wait per batch**: ~210 seconds (3.5 minutes)

### Error Handling Flow

```
Load Batch
    ↓
Try Submit Job
    ↓
Quota Error? ─[Yes]→ Wait & Retry (up to 3 times)
    ↓ [No]
Job Submitted
    ↓
Wait for Completion
    ↓
Transient Error? ─[Yes]→ Wait & Retry (up to 3 times)
    ↓ [No]
Job Complete
    ↓
Insert Tracking Records
    ↓
Quota Error? ─[Yes]→ Wait & Retry (up to 3 times)
    ↓ [No]
Success!
```

### Behavioral Changes

**Before:**
- Any exception during load → batch fails immediately
- No distinction between quota errors and other errors
- Generic error message
- Batch marked as failed, continue to next

**After:**
- Quota/rate limit errors → automatic retry with backoff
- Up to 3 retries per operation (submit, wait, tracking insert)
- Detailed logging of retry attempts and delays
- Only fails batch after exhausting all retries
- Non-retryable errors still fail immediately

## Benefits

1. **Improved Resilience**: Automatic recovery from temporary quota issues
2. **Better Visibility**: Clear logging shows when retries are happening
3. **Cost Efficiency**: Reduces need for manual re-runs
4. **Fault Isolation**: Quota errors don't necessarily fail entire task
5. **Smart Retry**: Only retries appropriate error types

## Configuration

The retry behavior can be adjusted by modifying the `max_retries` parameter in the helper functions:

```python
# Default is 3 retries
load_job = _submit_load_job_with_retry(..., max_retries=3)

# Can be increased for more aggressive retry:
load_job = _submit_load_job_with_retry(..., max_retries=5)
```

The retry delays are configured in the `retry_delays` list:

```python
retry_delays = [30, 60, 120]  # Current: 30s, 60s, 120s
```

## Monitoring

When quota errors occur, you'll see messages like:

```
  Quota/rate limit error: 429 POST https://bigquery.googleapis.com/...: Quota exceeded
  Retrying in 30 seconds (attempt 1/3)...
```

This allows you to:
- Track frequency of quota issues
- Adjust batch_size or WDL parallelism if needed
- Request quota increases from Google if persistent

## Future Enhancements (Not Yet Implemented)

### 4. WDL-Level Throttling
Add `maxConcurrency` to the scatter in `GvsImportGenomes.wdl`:

```wdl
scatter (fofn in DiscoverParquetFiles.file_fofns) {
  call LoadParquetFilesToBQ {
    input:
      project_id = project_id,
      dataset_name = dataset_name,
      fofn_file = fofn,
      batch_size = 10000,
      variants_docker = effective_variants_docker,
  }
}

# Would add:
# maxConcurrency: 50  # Limit to 50 parallel tables at a time
```

### 5. Adaptive Batch Sizing
Dynamically adjust `batch_size` based on quota errors encountered.

## Testing Recommendations

1. **Unit tests**: Mock BigQuery client to simulate quota errors
2. **Integration tests**: Use small dataset with intentionally low quotas
3. **Load tests**: Monitor retry frequency in production-scale runs
4. **Monitoring**: Track "batches_failed" in output statistics

## Related Files

- `scripts/variantstore/scripts/load_parquet_to_bq.py` - Main implementation
- `scripts/variantstore/wdl/GvsImportGenomes.wdl` - WDL workflow
- `docs/parquet_to_bigquery_loading_design.md` - Original design document
