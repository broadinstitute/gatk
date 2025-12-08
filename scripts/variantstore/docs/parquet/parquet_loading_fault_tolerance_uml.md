# Fault Tolerance Diagrams for Parquet Loading

## Mermaid Diagrams (GitHub/VS Code Compatible)


### 1. Architecture Diagram - Complete Workflow

```mermaid
flowchart TB
    subgraph Cromwell["Cromwell Workflow Engine"]
        WDL[GvsImportGenomes.wdl]
    end
    
    subgraph LoadDataTask["LoadData Tasks (Parallel)"]
        LD1[Create Parquet Files<br/>Sample 1-N]
        LD2[Create Parquet Files<br/>Sample N+1-2N]
        LD3[...]
    end
    
    subgraph GCS["GCS Storage"]
        VET["vet/000/, vet/001/, ..."]
        REF["ref_ranges/000/, ref_ranges/001/, ..."]
    end
    
    subgraph FaultTolerance["Fault Tolerance Layer"]
        Lifecycle[Configure Lifecycle<br/>14 day cleanup]
        Track[Create Tracking Table<br/>parquet_load_status]
        Discover[Discover Parquet Files<br/>Filter already loaded]
    end
    
    subgraph LoadTasks["LoadParquetFilesToBQ Tasks (Parallel)"]
        LT1[Load vet_000]
        LT2[Load vet_001]
        LT3[Load ref_ranges_000]
        LT4[...]
    end
    
    subgraph RetryLogic["Retry Mechanisms"]
        Det[Deterministic Job IDs<br/>SHA1 hash of batch]
        Quota[Quota Retry<br/>Exponential backoff]
        Local[Local State<br/>pending_jobs.json]
        Persist[Persistent State<br/>BigQuery tracking table]
    end
    
    subgraph BigQuery["BigQuery"]
        Tables[vet_XXX, ref_ranges_XXX<br/>tables]
        TrackTable[parquet_load_status<br/>table]
    end
    
    Verify[Verify All Loaded]
    
    WDL --> LoadDataTask
    LoadDataTask --> GCS
    GCS --> FaultTolerance
    FaultTolerance --> LoadTasks
    LoadTasks --> RetryLogic
    RetryLogic --> BigQuery
    BigQuery --> Verify
    
    style FaultTolerance fill:#9cf
    style RetryLogic fill:#fcf
    style BigQuery fill:#ffc
```


### 2. Complete Workflow - Activity Diagram

```mermaid
flowchart TD
    Start([GvsImportGenomes<br/>Workflow Starts]) --> CreateParquet[LoadData Task:<br/>Create Parquet Files]
    
    CreateParquet --> SetLoaded[SetIsLoadedColumn:<br/>Update sample_info]
    
    SetLoaded --> ConfigLife[ConfigureParquetLifecycle:<br/>Set 14-day deletion rule]
    
    ConfigLife --> CreateTrack[CreateParquetTrackingTable:<br/>Create parquet_load_status]
    
    CreateTrack --> Discover[DiscoverParquetFiles:<br/>List & group files]
    
    Discover --> CheckFiles{Files found?}
    
    CheckFiles -->|Yes| ParallelLoad[LoadParquetFilesToBQ Tasks<br/>Parallel scatter]
    CheckFiles -->|No| VerifyDone
    
    ParallelLoad --> LoadTable1[Load vet_000]
    ParallelLoad --> LoadTable2[Load vet_001]
    ParallelLoad --> LoadTable3[Load ref_ranges_000]
    ParallelLoad --> LoadTable4[...]
    
    LoadTable1 --> Batch1{Process batches<br/>in table}
    LoadTable2 --> Batch2{Process batches<br/>in table}
    LoadTable3 --> Batch3{Process batches<br/>in table}
    LoadTable4 --> Batch4{Process batches<br/>in table}
    
    Batch1 -->|Each batch| Submit1[Submit with<br/>deterministic job_id]
    Batch2 -->|Each batch| Submit2[Submit with<br/>deterministic job_id]
    Batch3 -->|Each batch| Submit3[Submit with<br/>deterministic job_id]
    Batch4 -->|Each batch| Submit4[Submit with<br/>deterministic job_id]
    
    Submit1 --> Track1[Update tracking table]
    Submit2 --> Track2[Update tracking table]
    Submit3 --> Track3[Update tracking table]
    Submit4 --> Track4[Update tracking table]
    
    Track1 --> Sync1{All batches<br/>complete?}
    Track2 --> Sync2{All batches<br/>complete?}
    Track3 --> Sync3{All batches<br/>complete?}
    Track4 --> Sync4{All batches<br/>complete?}
    
    Sync1 -->|Yes| Done1[Table complete]
    Sync2 -->|Yes| Done2[Table complete]
    Sync3 -->|Yes| Done3[Table complete]
    Sync4 -->|Yes| Done4[Table complete]
    
    Sync1 -->|No| Batch1
    Sync2 -->|No| Batch2
    Sync3 -->|No| Batch3
    Sync4 -->|No| Batch4
    
    Done1 --> VerifyDone
    Done2 --> VerifyDone
    Done3 --> VerifyDone
    Done4 --> VerifyDone
    
    VerifyDone[VerifyParquetLoading:<br/>Check all files loaded]
    
    VerifyDone --> CheckVerify{All files<br/>loaded?}
    
    CheckVerify -->|Yes| Success([SUCCESS])
    CheckVerify -->|No| Partial([PARTIAL:<br/>Report missing files])
    
    style ConfigLife fill:#9cf
    style CreateTrack fill:#9cf
    style Submit1 fill:#fcf
    style Submit2 fill:#fcf
    style Submit3 fill:#fcf
    style Submit4 fill:#fcf
```

---


### 3. Sequence Diagram - VM Preemption Recovery

```mermaid
sequenceDiagram
    participant C as Cromwell
    participant V1 as VM-1
    participant L1 as LoadParquet Task
    participant BQ as BigQuery
    participant V2 as VM-2
    participant L2 as LoadParquet Task (Retry)

    Note over C,BQ: Initial Run
    C->>V1: Start LoadParquetFilesToBQ
    activate V1
    V1->>L1: Execute Python Script
    activate L1
    
    L1->>L1: Generate deterministic job_id<br/>load_vet_001_abc123def
    L1->>L1: Write pending_jobs.json
    L1->>BQ: submit_load_job(job_id="load_vet_001_abc123def")
    activate BQ
    Note right of BQ: Job submitted<br/>and starts processing
    
    Note over V1,L1: VM Preemption!
    V1--xV1: VM destroyed<br/>Local disk lost<br/>pending_jobs.json gone
    deactivate L1
    deactivate V1
    
    Note over C,V2: Cromwell Retry on New VM
    C->>V2: Retry LoadParquetFilesToBQ
    activate V2
    V2->>L2: Execute Python Script (fresh start)
    activate L2
    
    L2->>L2: Generate same job_id<br/>load_vet_001_abc123def<br/>(deterministic)
    L2->>L2: No pending_jobs.json<br/>(fresh disk)
    
    L2->>BQ: submit_load_job(job_id="load_vet_001_abc123def")
    BQ-->>L2: Conflict Exception<br/>"Job already exists"
    Note right of BQ: BigQuery prevents<br/>duplicate submission
    
    L2->>BQ: get_job(job_id="load_vet_001_abc123def")
    BQ-->>L2: Job Status: RUNNING/DONE
    L2->>L2: Resume from existing job
    L2->>BQ: wait_for_completion()
    BQ-->>L2: Job completed successfully
    deactivate BQ
    
    L2->>BQ: INSERT INTO parquet_load_status
    L2->>C: Task completed successfully
    deactivate L2
    deactivate V2
```


### 4. Flowchart - Quota Error Handling with Exponential Backoff

```mermaid
flowchart TD
    Start([Start: Load Batch]) --> GenID[Generate deterministic job_id]
    GenID --> Submit[Submit BigQuery load job]
    
    Submit --> CheckSubmit{Job submission<br/>successful?}
    
    CheckSubmit -->|Yes| Submitted[Job submitted]
    CheckSubmit -->|Exception| CheckError{Exception type?}
    
    CheckError -->|TooManyRequests<br/>ServiceUnavailable<br/>InternalServerError| CheckRetry{Retry attempts<br/>< 3?}
    CheckError -->|Conflict:<br/>Job exists| FetchJob[Fetch existing job<br/>from BigQuery]
    CheckError -->|Other error| MarkFailed[Mark batch as failed]
    
    CheckRetry -->|Yes| Wait[Wait exponential backoff<br/>30s → 60s → 120s]
    CheckRetry -->|No| LogQuota[Log quota error]
    LogQuota --> MarkFailed
    
    Wait --> IncRetry[Increment retry count]
    IncRetry --> Submit
    
    FetchJob --> Submitted
    Submitted --> WaitComplete[Wait for job completion]
    
    WaitComplete --> CheckWait{Job wait<br/>successful?}
    
    CheckWait -->|Yes| JobDone[Job completed]
    CheckWait -->|Transient error| CheckWaitRetry{Wait retry<br/>attempts < 3?}
    CheckWait -->|Job failed| MarkFailed
    
    CheckWaitRetry -->|Yes| WaitRetry[Wait exponential backoff]
    CheckWaitRetry -->|No| MarkFailed
    
    WaitRetry --> Reload[Reload job state]
    Reload --> WaitComplete
    
    JobDone --> CheckErrors{Job has<br/>errors?}
    CheckErrors -->|Yes| MarkFailed
    CheckErrors -->|No| RecordSuccess[Record success in<br/>tracking table]
    
    RecordSuccess --> UpdatePending[Update pending_jobs.json]
    UpdatePending --> NextBatch[Continue to next batch]
    MarkFailed --> NextBatch
    
    NextBatch --> End([End])
    
    style CheckRetry fill:#ff9
    style CheckWaitRetry fill:#ff9
    style Wait fill:#9f9
    style WaitRetry fill:#9f9
    style RecordSuccess fill:#9f9
```


### 5. State Diagram - BigQuery Load Job Lifecycle

```mermaid
stateDiagram-v2
    [*] --> PENDING: create_job_id()
    
    PENDING --> SUBMITTING: submit_job()
    PENDING --> RESUMING: job_exists_in_bigquery()
    
    SUBMITTING --> SUBMITTED: success
    SUBMITTING --> QUOTA_RETRY: quota_error
    SUBMITTING --> FAILED: non_retryable_error
    
    QUOTA_RETRY --> WAITING: sleep(backoff_delay)
    WAITING --> SUBMITTING: retry_attempt < max_retries
    WAITING --> FAILED: retry_attempt >= max_retries
    
    RESUMING --> SUBMITTED: fetch_existing_job()
    
    SUBMITTED --> RUNNING: job.state = RUNNING
    SUBMITTED --> DONE: job.state = DONE
    
    RUNNING --> WAITING_COMPLETION: wait_for_result()
    WAITING_COMPLETION --> DONE: job_completed
    WAITING_COMPLETION --> TRANSIENT_RETRY: transient_error
    WAITING_COMPLETION --> FAILED: job_failed
    
    TRANSIENT_RETRY --> WAITING_RETRY: sleep(backoff_delay)
    WAITING_RETRY --> WAITING_COMPLETION: retry_attempt < max_retries
    WAITING_RETRY --> FAILED: retry_attempt >= max_retries
    
    DONE --> SUCCESS: no_job_errors
    DONE --> FAILED: job_has_errors
    
    SUCCESS --> TRACKED: insert_tracking_record()
    TRACKED --> [*]
    FAILED --> [*]
    
    note right of QUOTA_RETRY
        Exponential backoff:
        30s → 60s → 120s
    end note
    
    note right of TRACKED
        File recorded in
        parquet_load_status table
    end note
```


## Summary

These Mermaid diagrams provide comprehensive visualization of the fault tolerance mechanisms in the parquet loading system:

1. **Architecture Diagram** - Shows the four-layer fault tolerance architecture
2. **Activity Diagram** - Depicts the complete workflow from VCF processing through verification
3. **Sequence Diagram** - Shows VM preemption scenario with deterministic job ID recovery
4. **Flowchart** - Illustrates the complete retry and recovery decision flow
5. **State Diagram** - Tracks the lifecycle of a load job through all possible states

They render natively in GitHub and VS Code with the Markdown Preview Mermaid extension, making them ideal for documentation purposes.
