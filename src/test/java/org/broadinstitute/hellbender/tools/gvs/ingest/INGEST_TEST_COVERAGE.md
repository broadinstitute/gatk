# GVS Ingest Test Coverage

This document maps the existing tests for `CreateVariantIngestFiles` and its
subordinate `*Creator` / `*Writer` classes, and identifies known gaps.

---

## Existing tests

### `CreateVariantIngestFilesIntegrationTest`

The only integration-level test for the ingest pipeline. Three test methods run
the full `CreateVariantIngestFiles` tool end-to-end against a single reblocked
NA12878 GVCF on chr20, exercising the **PARQUET** output path:

| Test method | What differs |
|---|---|
| `testCreateDefaultVariantIngestFiles` | Default GQ filtering (ignore GQ60+) |
| `testIgnoreGQ0CreateVariantIngestFiles` | `--ref-block-gq-to-ignore ZERO` |
| `testIgnoreGQ40AndAboveCreateVariantIngestFiles` | `--ref-block-gq-to-ignore FORTY --ignore-above-gq-threshold true` |

All three methods validate:
- Exact output filenames for VET, ref_ranges, and ploidy Parquet files
- That the output directory contains **exactly** those three files (plus Hadoop
  `.crc` sidecars which are excluded from the check)
- Parquet content converted to TSV and compared against golden TSV files

`VetCreator`, `RefRangesCreator`, and `SamplePloidyCreator` are exercised
indirectly through this test, but only along the PARQUET happy path for a
diploid autosomal chr20 sample.

---

### `VetCreatorUnitTest`

Two unit tests targeting `VetCreator` directly:

| Test method | What it covers |
|---|---|
| `testParquetOutputFileNaming` | `VetCreator.getOutputFileName()` produces the expected `vet_<table>_<id>_<vcf>.parquet` string |
| `testErrorFile` | Constructing a `VetCreator` when the output file already exists throws a wrapped `FileAlreadyExistsException` |

---

### `VetFieldEnumUnitTest`

Five unit tests for `VetFieldEnum.getColumnValue()` — the core annotation
extraction logic used when writing VET data. This is the most thorough
unit-level coverage in the ingest package:

| Test method | What it covers |
|---|---|
| `testAlleleSpecificOneAlt` | AS annotations with a single alt allele (`0/1`) |
| `testAlleleSpecificTwoAlts` | AS annotation trimming for het-var (`1/2`) genotypes |
| `testNonAlleleSpecificSingleAlt` | Non-AS fallback with a single alt allele |
| `testNonAlleleSpecificTwoAlts` | Non-AS fallback with two alt alleles |
| `testForceNonASLoading` | `forceLoadingFromNonAlleleSpecific=true` overrides AS annotations |

---

### `PendingBQWriterTest`

Tests the underlying BigQuery Write API infrastructure used by all BQ writers.
Tagged `groups = {"cloud"}` — requires live GCP credentials and is **not** run
in normal CI. There are no offline unit tests for any of the BQ writer classes
(`VetBQWriter`, `RefRangesBQWriter`, `SamplePloidyBQWriter`).

---

## Coverage gaps

| Class | Unit tests | Integration coverage | Notes |
|---|---|---|---|
| `RefRangesCreator` | **None** | ✓ (happy path only) | See detail below |
| `SamplePloidyCreator` | **None** | ✓ (happy path only) | Mixed-ploidy logic untested |
| `VcfHeaderLineScratchCreator` | **None** | **None** | `enableVCFHeaders` is `false` in all integration tests |
| `VetParquetFileWriter` | **None** | ✓ via integration test | |
| `RefRangesParquetFileWriter` | **None** | ✓ via integration test | |
| `SamplePloidyParquetFileWriter` | **None** | ✓ via integration test | |
| `AbstractParquetFileWriter` | **None** | ✓ via integration test | |
| `VetBQWriter` | **None** | Cloud-only, not in CI | |
| `RefRangesBQWriter` | **None** | Cloud-only, not in CI | |
| `SamplePloidyBQWriter` | **None** | Cloud-only, not in CI | |
| Output filename logic (`RefRangesCreator`) | **None** | ✓ via integration test | Contrast: `VetCreator.getOutputFileName` has a dedicated unit test |
| Output filename logic (`SamplePloidyCreator`) | **None** | ✓ via integration test | |

### `RefRangesCreator` — notable untested branches

`RefRangesCreator` has the most complex logic in the ingest package and the
weakest direct test coverage. The following code paths are exercised only
through the three integration test cases (single diploid autosomal chr20 sample)
or not at all:

- **Compressed reference storage** (`storeCompressedReferences=true`) — the
  `writeCompressed` path in both `apply()` and `writeMissingPositions()` is
  never invoked by any test.
- **PAR region handling** — `PloidyUtils.doesVariantOverlapPAR()` returns
  `false` for all chr20 positions; PAR-overlap suppression of ploidy recording
  is untested.
- **`dropAboveGqThreshold` / `getGQStateEnumGreaterThan`** — the
  `testIgnoreGQ40AndAboveCreateVariantIngestFiles` integration test exercises
  this path, but only at the tool level with a fixed input file.
- **Adjacent-interval merging in `setCoveredInterval`** — the `overlapping`
  branch is exercised by the integration test input but has no targeted unit
  test.
- **`writeMissingPositions` chunking** — blocks longer than
  `MAX_REFERENCE_BLOCK_BASES` (1000 bp) are split into multiple writes; no
  test constructs a gap large enough to trigger more than one iteration.
- **`writeMissingIntervals` with a fully uncovered chromosome** — not exercised
  by any test.

### `SamplePloidyCreator` — notable untested branches

- **Mixed-ploidy normalization** — the block starting at `if
  (ploidiesWithCounts.size() != 1)` (choosing a dominant ploidy among competing
  ploidies and logging a warning) is never triggered; the integration test
  input is a clean diploid sample throughout.
- **Mixed-ploidy error threshold** — the `UserException` thrown when the
  second-best ploidy exceeds `MIXED_PLOIDY_ERROR_CUTOFF` (5%) is untested.

