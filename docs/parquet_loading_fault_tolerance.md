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
│  • GetUningestedSampleIds: Auto-cleanup of temp tables ✓                          │
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
│                          LAYER 4: PERSISTENT STATE TRACKING                        │
├─────────────────────────────────────────────────────────────────────────────────────┤
│                                                                                     │
│  BigQuery Tracking Table (parquet_load_status):                                    │
│  ┌─────────────────────────────────────────────────────────────────────────────┐   │
│  │  CREATE TABLE parquet_load_status (                                        │   │
│  │    file_path STRING PRIMARY KEY,    ◄── Prevents duplicates               │   │
│  │    table_name STRING,                                                       │   │
│  │    load_timestamp TIMESTAMP,                                                │   │
│  │    load_job_id STRING,                                                      │   │
│  │    rows_loaded INT64                                                        │   │
│  │  )                                                                          │   │
│  └─────────────────────────────────────────────────────────────────────────────┘   │
│                                                                                     │
│  MERGE Logic for Idempotency:                                                      │
│  ┌─────────────────────────────────────────────────────────────────────────────┐   │
│  │  MERGE parquet_load_status T                                               │   │
│  │  USING new_records S                                                        │   │
│  │  ON T.file_path = S.file_path                                              │   │
│  │  WHEN NOT MATCHED THEN INSERT (...)  ◄── Only insert new files           │   │
│  └─────────────────────────────────────────────────────────────────────────────┘   │
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
│  Workflow Restart → Check Tracking Table → Skip Already Loaded Files               │
│         │                    │                         │                          │
│         ▼                    ▼                         ▼                          │
│  Generate FOFNs → Filter Out Loaded → Only Load Missing → Update Tracking         │
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

### Layer 4: Persistent State Tracking
- **BigQuery Tracking Table**: Permanent record of loaded files with MERGE-based idempotency
- **Cross-Run Recovery**: Filters out already-loaded files on workflow restart

## Key Benefits

1. **Complete Idempotency**: Safe to restart workflow any number of times
2. **Granular Recovery**: Can recover from failures at file, batch, or workflow level  
3. **Quota Resilience**: Automatic handling of BigQuery rate limits and quotas
4. **Cost Efficiency**: Uses preemptible VMs while maintaining reliability
5. **Progress Preservation**: Never loses work due to transient failures

## Failure Handling

The system gracefully handles three primary failure modes:
- **VM Preemption**: Deterministic job IDs enable seamless recovery on new VMs
- **Quota Exceeded**: Exponential backoff prevents cascade failures while respecting limits
- **Partial Batch Failure**: Individual batch failures don't stop overall progress

This multi-layered approach ensures the pipeline can handle the scale and complexity of loading hundreds of thousands of parquet files reliably in a cloud environment.