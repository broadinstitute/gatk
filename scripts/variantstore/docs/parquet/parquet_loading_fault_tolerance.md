# Fault Tolerance Mechanisms in Parquet Loading Pipeline

```
┌─────────────────────────────────────────────────────────────────────────────────────┐
│                           FAULT TOLERANCE MECHANISMS                                │
└─────────────────────────────────────────────────────────────────────────────────────┘

┌──────────────────┐    ┌──────────────────┐    ┌──────────────────┐
│  VM PREEMPTION   │    │  QUOTA ERRORS    │    │  PARTIAL FAILS   │
│     RECOVERY     │    │     HANDLING     │    │    RECOVERY      │
└──────────────────┘    └──────────────────┘    └──────────────────┘
         │                        │                        │
         ▼                        ▼                        ▼

┌─────────────────────────────────────────────────────────────────────────────────────┐
│                              LAYER 1: WORKFLOW LEVEL                               │
├─────────────────────────────────────────────────────────────────────────────────────┤
│                                                                                     │
│  WDL Runtime Configuration:                                                         │
│  ┌─────────────────────────────────────────────────────────────────────────────┐   │
│  │  runtime {                                                                  │   │
│  │    preemptible: 5        ◄── Up to 5 preemption attempts                  │   │
│  │    maxRetries: 3         ◄── 3 additional retries for other failures      │   │
│  │    disks: "local-disk"   ◄── Fresh disk on each retry                     │   │
│  │  }                                                                          │   │
│  └─────────────────────────────────────────────────────────────────────────────┘   │
│                                                                                     │
│  Idempotent Task Design:                                                            │
│  • DiscoverParquetFiles: Queries BigQuery directly for already-loaded pairs ✓      │
│  • LoadParquetFilesToBQ: Deterministic job IDs for resume ✓                       │
│                                                                                     │
└─────────────────────────────────────────────────────────────────────────────────────┘
                                        │
                                        ▼
┌─────────────────────────────────────────────────────────────────────────────────────┐
│                            LAYER 2: BIGQUERY JOB LEVEL                             │
├─────────────────────────────────────────────────────────────────────────────────────┤
│                                                                                     │
│  Deterministic Job IDs:                                                             │
│  ┌─────────────────────────────────────────────────────────────────────────────┐   │
│  │  job_id = f"load_{table_name}_{sha1_hash_of_batch}"                        │   │
│  │                                                                             │   │
│  │  VM Dies → New VM → Same job_id → BigQuery says:                          │   │
│  │  "Conflict: Job already exists" → Fetch existing job ✓                    │   │
│  └─────────────────────────────────────────────────────────────────────────────┘   │
│                                                                                     │
│  Load Job Retry with Exponential Backoff:                                          │
│  ┌─────────────────────────────────────────────────────────────────────────────┐   │
│  │  try:                                                                       │   │
│  │    submit_job()                                                             │   │
│  │  except TooManyRequests:    ◄── Quota exceeded                            │   │
│  │    sleep(30s) → retry                                                       │   │
│  │  except ServiceUnavailable: ◄── Temporary unavailable                     │   │
│  │    sleep(60s) → retry                                                       │   │
│  │  except InternalServerError: ◄── Transient error                          │   │
│  │    sleep(120s) → retry                                                      │   │
│  └─────────────────────────────────────────────────────────────────────────────┘   │
│                                                                                     │
└─────────────────────────────────────────────────────────────────────────────────────┘
                                        │
                                        ▼
┌─────────────────────────────────────────────────────────────────────────────────────┐
│                            LAYER 3: LOCAL STATE TRACKING                           │
├─────────────────────────────────────────────────────────────────────────────────────┤
│                                                                                     │
│  Pending Jobs File (pending_jobs.json):                                            │
│  ┌─────────────────────────────────────────────────────────────────────────────┐   │
│  │  {                                                                          │   │
│  │    "load_vet_001_abc123": {                                                 │   │
│  │      "files": ["gs://bucket/file1.parquet", ...],                          │   │
│  │      "location": "us-central1",                                             │   │
│  │      "table_name": "vet_001"                                                │   │
│  │    }                                                                        │   │
│  │  }                                                                          │   │
│  └─────────────────────────────────────────────────────────────────────────────┘   │
│                                                                                     │
│  State Lifecycle:                                                                   │
│  Submit Job → Write to File → Job Completes → Remove from File                     │
│       │              │              │              │                              │
│       │              ▼              │              ▼                              │
│  On VM Restart:  File Exists?   Job Still     File Cleanup                        │
│       │              │         Running?           │                              │
│       ▼              ▼              ▼              ▼                              │
│  Fresh Disk     [File Lost]    Check BigQuery   Continue                          │
│  No File    →   [Use Job IDs]  →   Status    →   Processing                       │
│                                                                                     │
└─────────────────────────────────────────────────────────────────────────────────────┘
                                        │
                                        ▼
┌─────────────────────────────────────────────────────────────────────────────────────┐
│                        LAYER 4: PERSISTENT STATE (BIGQUERY DATA)                   │
├─────────────────────────────────────────────────────────────────────────────────────┤
│                                                                                     │
│  BigQuery as Source of Truth:                                                      │
│  ┌─────────────────────────────────────────────────────────────────────────────┐   │
│  │  SELECT table_name, SAFE_CAST(partition_id AS INT64) AS sample_id          │   │
│  │  FROM INFORMATION_SCHEMA.PARTITIONS                                        │   │
│  │  WHERE REGEXP_CONTAINS(table_name, "^vet_[0-9]+$|^ref_ranges_[0-9]+$")    │   │
│  │    AND total_logical_bytes > 0                                              │   │
│  │                                                                             │   │
│  │  UNION ALL                                                                  │   │
│  │                                                                             │   │
│  │  SELECT DISTINCT "sample_chromosome_ploidy", sample_id                     │   │
│  │  FROM sample_chromosome_ploidy                                             │   │
│  └─────────────────────────────────────────────────────────────────────────────┘   │
│                                                                                     │
│  Already-loaded (table_name, sample_id) pairs are skipped on workflow restart.     │
│  No separate tracking table is written to or consulted.                            │
│                                                                                     │
└─────────────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────────────────┐
│                               FAILURE SCENARIOS                                    │
├─────────────────────────────────────────────────────────────────────────────────────┤
│                                                                                     │
│  ┌─────────────────┐    ┌─────────────────┐    ┌─────────────────┐                │
│  │  VM PREEMPTED   │    │  QUOTA EXCEEDED │    │  BATCH FAILURE  │                │
│  └─────────────────┘    └─────────────────┘    └─────────────────┘                │
│           │                        │                        │                     │
│           ▼                        ▼                        ▼                     │
│  • New VM starts            • Wait 30s-120s         • Mark batch failed          │
│  • Fresh disk               • Retry same job        • Continue next batch        │
│  • No pending_jobs.json     • Exponential backoff   • Don't fail entire task    │
│  • Same job IDs generated   • Up to 3 attempts      • Partial success possible   │
│  • BigQuery: "Job exists"   • Then fail batch       • Report in final stats     │
│  • Fetch existing job       • Continue pipeline     • Retry via workflow rerun   │
│  • Resume from there ✓      • Other batches OK ✓    • Already loaded = skipped ✓ │
│                                                                                     │
└─────────────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────────────────┐
│                                  RECOVERY FLOW                                     │
├─────────────────────────────────────────────────────────────────────────────────────┤
│                                                                                     │
│  Workflow Restart → Query BigQuery Partitions → Skip Already Loaded Pairs          │
│         │                    │                         │                          │
│         ▼                    ▼                         ▼                          │
│  Generate FOFNs → Filter Out Loaded → Only Load Missing → Verify via BQ           │
│                                                                                     │
│  Result: 100% Idempotent - Safe to retry any number of times ✓                    │
│                                                                                     │
└─────────────────────────────────────────────────────────────────────────────────────┘
```

## Summary

This diagram illustrates the comprehensive fault tolerance mechanisms implemented in the parquet loading pipeline. The system provides resilience at four distinct layers:

### Layer 1: Workflow Level
- **WDL Configuration**: Up to 5 preemptible VM attempts plus 3 additional retries
- **Idempotent Design**: Tasks can be safely rerun without side effects

### Layer 2: BigQuery Job Level
- **Deterministic Job IDs**: SHA1-based IDs allow job recovery across VM restarts
- **Exponential Backoff**: Automatic retry for quota and transient errors (30s → 60s → 120s)

### Layer 3: Local State Tracking
- **Pending Jobs File**: Tracks in-flight jobs on local disk for within-VM recovery
- **Graceful Degradation**: When file is lost (VM preemption), falls back to job ID mechanism

### Layer 4: Persistent State (BigQuery Data)
- **BigQuery as Source of Truth**: `INFORMATION_SCHEMA.PARTITIONS` and the `sample_chromosome_ploidy` table itself provide an authoritative view of what has been loaded
- **Cross-Run Recovery**: Already-loaded `(table_name, sample_id)` pairs are filtered out on workflow restart — no separate tracking table required

## Key Benefits

1. **Complete Idempotency**: Safe to restart workflow any number of times
2. **Granular Recovery**: Can recover from failures at file, batch, or workflow level
3. **Quota Resilience**: Automatic handling of BigQuery rate limits and quotas; elimination of the tracking table removes a significant source of MERGE-based quota pressure
4. **Cost Efficiency**: Uses preemptible VMs while maintaining reliability
5. **Progress Preservation**: Never loses work due to transient failures
6. **Authoritative State**: Idempotency is based on what BigQuery actually contains, not what a secondary table records
