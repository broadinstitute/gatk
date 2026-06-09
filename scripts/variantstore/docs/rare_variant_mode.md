# Rare Variant Mode

## Overview

Rare variant mode is an optional operational mode for GVS designed to support
research on rare variants. Standard GVS extracts do not emit read depth (DP)
for reference genotypes, and GQ0 reference blocks are emitted as `./. ` no-calls.
Both of these defaults can make it difficult to distinguish true low-confidence
reference calls from missing data in downstream rare variant analyses.

Rare variant mode addresses this by storing DP in the `ref_ranges` tables during
ingest and propagating it through to the extracted VCF, and by emitting GQ0
reference blocks as genotype calls rather than no-calls — using DP presence as
the mechanism to distinguish real GQ0 calls from true missing data.

## Enabling Rare Variant Mode

Set `rare_variant_mode = true` in the top-level workflow input. This single flag
controls all three pipeline stages described below. It is supported in:

- `GvsJointVariantCalling.wdl` (Beta callsets)
- `GvsBulkIngestGenomes.wdl` / `GvsImportGenomes.wdl` (ingest stage)
- `GvsPrepareRangesCallset.wdl` (prepare stage)
- `GvsExtractCallset.wdl` (extract stage)

The flag can also be applied to individual stages independently using the
lower-level flags described below, though this requires understanding the
dependency described in the GQ0 disambiguation section.

## Pipeline Behavior

### Stage 1: Ingest

**Flag:** `--include-ref-ranges-dp` on `CreateVariantIngestFiles`
(set automatically when `rare_variant_mode = true` in `GvsImportGenomes.wdl`)

During ingest, the DP value from each reference block genotype in the source
GVCF is read and stored in the `ref_ranges_%` BigQuery tables as an optional
`dp` column. The `dp` column is present in the `ref_ranges` schema for all GVS
datasets (the column is always defined as `NULLABLE`), but it is only populated
when this flag is set.

Notes:
- If the source GVCF does not include DP in a reference block genotype, the `dp`
  field for that row is stored as `null` even in rare variant mode.
- Synthetic reference rows added to fill coverage gaps (positions not covered by
  any reference block in the GVCF) always have `null` dp.
- True no-calls (positions with missing genotype data) are stored as GQ0 reference
  blocks with `null` dp. This is an existing GVS convention and is not changed by
  rare variant mode. See GQ0 Disambiguation below.

### Stage 2: Prepare

**Flag:** `--include_ref_ranges_dp` on `create_ranges_cohort_extract_data_table.py`
(set automatically when `rare_variant_mode = true` in `GvsPrepareRangesCallset.wdl`)

During the prepare stage, the `dp` column is included in the cohort extract
reference table that is materialized for the extract stage. Without this flag,
`dp` is not selected from the `ref_ranges` source tables and is not available
to the extract tool.

### Stage 3: Extract

**Flags:** `--emit-ref-ranges-dp` and `--emit-gq0-ref-blocks` on `ExtractCohort`
(both set automatically when `rare_variant_mode = true` in `GvsExtractCallset.wdl`)

During extract, two behaviors are enabled:

1. **DP emission**: When `--emit-ref-ranges-dp` is set, the `DP` FORMAT field is
   included in the VCF header and populated for reference genotypes where a
   non-null DP value is available from the prepare table.

2. **GQ0 block emission**: When `--emit-gq0-ref-blocks` is set, GQ0 reference
   entries are emitted as `0/0` genotype calls with GQ=0, rather than as `./. `
   no-calls. This only applies to entries where DP is non-null. See GQ0
   Disambiguation below.

## GQ0 Disambiguation

This is the most subtle aspect of rare variant mode and is critical to understand.

GVS represents true no-calls (positions with missing genotype data in the source
GVCF) as GQ0 reference blocks during ingest. This is an existing convention.
The consequence is that at extract time, GQ0 entries in the `ref_ranges` tables
can represent one of two things:

1. A **true GQ0 reference block**: a real genotype call with low confidence (GQ
   between 0 and 10), which came with a valid DP value from the source GVCF.
2. A **no-call transcribed as GQ0**: a position where the source GVCF had missing
   genotype data, stored as GQ0 by GVS convention. These entries do not have a
   valid DP.

In rare variant mode, **DP presence is the mechanism that distinguishes these two
cases**. GQ0 entries with non-null DP are emitted as `0/0` GQ0 genotype calls.
GQ0 entries with null DP remain as `./. ` no-calls. Downstream analyses can
therefore use DP presence to determine whether a GQ0 entry represents a real call.

This design was validated with the rare variant research team, who confirmed that
this distinction is sufficient for their use case.

## Constraints and Caveats

**`--emit-gq0-ref-blocks` requires ingest-time DP storage.** If a callset was
ingested without `--include-ref-ranges-dp`, all DP values in the prepare table
will be null, and `--emit-gq0-ref-blocks` will have no effect — GQ0 entries will
all remain as `./. ` no-calls. The `rare_variant_mode` composite flag ensures
these are always set together; if using the lower-level flags independently,
ensure the callset was ingested with `--include-ref-ranges-dp` before enabling
`--emit-gq0-ref-blocks` at extract time.

**DP may be null even in rare variant mode.** A non-null DP in the extract output
means the source GVCF had a valid DP for that reference block. A null DP means
either (a) the source GVCF did not include DP, (b) the row is a gap-fill synthetic
reference row, or (c) the row represents a no-call.

**There is currently no queryable metadata recording whether a dataset was
ingested in rare variant mode.** The `dp` column is always present in the schema
but is all-null for non-rare-variant datasets. Querying for non-null dp values is
not a reliable indicator of rare variant mode because DP may be null even in
rare-variant-mode datasets (see above). Until explicit metadata tracking is added,
the ingest workflow inputs or job history should be consulted to determine the
ingest-time configuration.

**Pre-existing datasets cannot retroactively gain DP.** Adding DP to a callset
that was already ingested requires re-ingesting the affected samples with
`--include-ref-ranges-dp`. The `dp` column in the BigQuery table will accept the
data without a schema change (the column is already defined as `NULLABLE`).

**`GvsExtractCallsetPgenMerged.wdl` does not currently support rare variant mode.**
The `rare_variant_mode` flag has not been plumbed into the PGEN merged extract
workflow.
