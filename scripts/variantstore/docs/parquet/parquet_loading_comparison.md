# LoadData Task Comparison: Old vs New Approach

This document compares the original BigQuery Write API approach with the new Parquet-based loading approach for the GVS ingest process.

## Old LoadData Task (BigQuery Write API)

```mermaid
flowchart TD
    Start([Start Original LoadData Task]) --> Localize[Localize gvcf file]
    Localize --> CheckStatus[Check and update status tables]
    CheckStatus --> ProcessHeader[Process header info]
    ProcessHeader --> StreamRows[Stream each row of the GVCF<br/>separating into two streams<br/>computing approximate allele-specific QC values]
    StreamRows --> WriteAPI[Write data to BigQuery<br/>using Write API<br/>as it's being processed]
    WriteAPI --> Finalize[Finish, finalizing the streams<br/>so the full write is committed]
    Finalize --> FinalizeHeader[Finalize header and ploidy data]
    FinalizeHeader --> UpdateStatus[Update status tables]
    UpdateStatus --> End([End])
    
    style WriteAPI fill:#ffcccc
    style Finalize fill:#ffcccc
    style Start fill:#e1f5e1
    style End fill:#e1f5e1
```

**Key Characteristics:**
- **Direct streaming** to BigQuery using the Write API
- **Synchronous processing** - data written as rows are processed
- **Stream finalization** required to commit the full write
- **Resource intensive** - Write API incurs costs and quota usage
- **Single-phase operation** - ingest and load happen together

---

## New LoadData Task (Parquet-based)

```mermaid
flowchart TD
    Start([Start New LoadData Task]) --> Localize[Localize files]
    Localize --> CheckStatus[Check and update status tables]
    CheckStatus --> ProcessHeader[Process header info]
    ProcessHeader --> StreamRows[Stream each row of the GVCF<br/>separating into two streams<br/>computing approximate allele-specific QC values]
    StreamRows --> WriteParquet[Write data to local parquet files on VM]
    WriteParquet --> FinalizeParquet[Finish, finalizing the parquet files]
    FinalizeParquet --> FinalizeHeader[Finalize header and ploidy data<br/>*not implemented yet, trivial to do]
    FinalizeHeader --> UpdateStatus[Update status tables]
    UpdateStatus --> CopyGCS[Copy parquet files to temporary GCS bucket<br/>passed in by user]
    CopyGCS --> End([End])
    
    style WriteParquet fill:#ccffcc
    style FinalizeParquet fill:#ccffcc
    style CopyGCS fill:#ccffcc
    style FinalizeHeader fill:#ffffcc
    style Start fill:#e1f5e1
    style End fill:#e1f5e1
```

**Key Characteristics:**
- **Local file writing** - data written to parquet files on VM disk
- **Asynchronous processing** - ingest separated from BigQuery loading
- **File-based output** - parquet format enables flexible loading strategies
- **Cost efficient** - avoids expensive Write API, uses batch loading instead
- **Two-phase operation** - ingest produces files, separate task loads to BigQuery
- **Temporary storage** - files staged in GCS bucket for loading

**Note:** Header and ploidy data finalization is not yet implemented in the new flow, but is a trivial addition.

