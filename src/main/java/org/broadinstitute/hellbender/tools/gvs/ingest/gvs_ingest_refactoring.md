# GVS Ingest Package Refactoring

## Overview

The `org.broadinstitute.hellbender.tools.gvs.ingest` package has been refactored to support multiple output formats through a consistent interface-based architecture. This refactoring separates concerns by format type and introduces polymorphism for easier addition of new output formats.

## Package Structure

```
src/main/java/org/broadinstitute/hellbender/tools/gvs/ingest/
├── VetWriter.java                    (interface)
├── RefRangesWriter.java              (interface)
├── SamplePloidyWriter.java           (interface)
├── bq/
│   ├── AbstractBQWriter.java         (abstract base)
│   ├── VetBQWriter.java
│   ├── RefRangesBQWriter.java
│   └── SamplePloidyBQWriter.java
└── parquet/
    ├── AbstractParquetFileWriter.java (abstract base)
    ├── VetParquetFileWriter.java
    ├── RefRangesParquetFileWriter.java
    ├── SamplePloidyParquetFileWriter.java
    └── HeaderParquetFileWriter.java
```

## Key Design Elements

### Writer Interfaces

Three core interfaces define the write operations for GVS data types:

- **`VetWriter`** - Writes variant information (location, sample ID, variant context)
- **`RefRangesWriter`** - Writes reference ranges (location, length, state, compressed data)
- **`SamplePloidyWriter`** - Writes chromosome ploidy information

All interfaces extend `Closeable` and include a default `commitData()` no-op method that implementations can override.

### Abstract Base Classes

Two abstract base classes provide shared infrastructure for their respective output formats:

- **`AbstractBQWriter`** - Handles BigQuery Write API operations via `PendingBQWriter`
  - Provides `write(JSONObject)` for adding rows
  - Implements `commitData()` to flush and commit write streams

- **`AbstractParquetFileWriter`** - Handles Parquet file operations
  - Uses composition pattern with `ParquetWriteSupport`
  - Provides `write(JSONObject)` for writing records
  - Manages compression and schema configuration

### Format-Specific Implementations

#### BigQuery (`bq/`)
- `VetBQWriter`, `RefRangesBQWriter`, `SamplePloidyBQWriter`
- Extend `AbstractBQWriter` and implement respective interfaces
- Write data directly to BigQuery via the Write API

#### Parquet (`parquet/`)
- `VetParquetFileWriter`, `RefRangesParquetFileWriter`, `SamplePloidyParquetFileWriter`
- Extend `AbstractParquetFileWriter` and implement respective interfaces
- Write data to local Parquet files for later bulk loading


## Benefits

1. **Separation of Concerns**: Output format logic is isolated in dedicated packages
2. **Extensibility**: New formats can be added by implementing the writer interfaces
3. **Code Reuse**: Abstract base classes eliminate duplication across similar implementations
4. **Type Safety**: Interface contracts ensure consistent APIs across formats
5. **Testability**: Each implementation can be tested independently

## Note

VCF header writing does not yet have a dedicated interface; that should be added when VS-1803 (support writing VCF header data with Parquet) is implemented.
