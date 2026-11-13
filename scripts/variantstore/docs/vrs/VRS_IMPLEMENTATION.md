# VRS Allele ID Implementation in GATK Variant Store

## Overview

This document describes the VRS (Variant Representation Specification) allele ID feature integrated into GATK's variant store ingest pipeline. VRS IDs provide standardized, globally-unique identifiers for genomic alleles, enabling interoperability across bioinformatics tools and databases.

**Current status**: Implemented from the Java ingest tool through WDL orchestration to BigQuery loading, behind `enable_vrs_ids` (default `false`). The IDs have not yet been checked against vrs-python beyond a single normalizer case (VS-2059), and `alt_allele` does not carry them yet (VS-2061).

---

## Feature Description

The VRS ID feature computes standardized identifiers (ga4gh:VA.* format) for each variant allele during ingestion. When enabled, the pipeline:

1. **Normalizes** each allele according to VRS 2.0.1 specification (left-alignment, decomposition, state resolution)
2. **Computes** unique identifiers from the normalized representation using SHA-512 digest
3. **Stores** IDs in both the VET (variant extraction table) and a dedicated vrs_allele table
4. **Supports** heterozygous variants (1/2 genotypes) with multiple alleles via array fields

This enables downstream queries to identify and retrieve specific alleles across samples using standardized identifiers.

---

## Key Implementation Classes

### Java Components (`src/main/java/org/broadinstitute/hellbender/tools/gvs/`)

| Class                        | Purpose                                                                                                                                                                     |
|------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **CreateVariantIngestFiles** | Entry point for variant ingest; accepts `--enable-vrs-ids` flag (default: false); instantiates VrsAlleleCreator when enabled                                                |
| **VetCreator**               | Extracts variants and computes VRS IDs via `VrsNormalizer` and `VrsIdComputer` during `apply()`; stores IDs in VET table and delegates allele records to `VrsAlleleCreator` |
| **VrsAlleleCreator**         | Orchestrates writing of VrsAlleleRecord objects, tagged with the sample id, to vrs_allele_scratch output (Parquet or BigQuery); one instance per sample                     |
| **VrsNormalizer**            | Implements VRS 2.0.1 normalization steps (trimming, rolling, expansion, state resolution) to produce canonical allele representation                                        |
| **VrsIdComputer**            | Generates ga4gh:VA.* and ga4gh:SL.* identifiers from normalized alleles using RFC 8785 canonical JSON + SHA-512 digest                                                      |
| **VrsAlleleRecord**          | One `vrs_allele` row: the allele and location IDs, refget accession, normalized start and end, and the state                                                                |

### Test Classes

- `VrsAlleleCreatorUnitTest`: file naming, Parquet output, no-op scenarios
- `VrsNormalizerTest`: trimming, rolling, expansion and state resolution for SNPs, indels and substitutions
- `VrsIdComputerTest`: deterministic ID generation and format; it does not compare IDs with known values (VS-2059)
- `VrsAlleleRecordTest`, `RefGetAccessionProviderTest`
- `VetCreatorUnitTest`: a spanning-deletion `*` ALT gets the placeholder and keeps `vrs_allele_ids` aligned with `alt`
- `CreateVariantIngestFilesIntegrationTest`: end-to-end VRS ID computation from a gVCF. It needs a local reference at a hardcoded path and skips without it.

---

## WDL Integration

VRS ID flow is wired through three WDL workflows with consistent opt-in semantics:

| Workflow                     | Role                                                                                                                                                                                                                                  |
|------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **GvsCreateTables.wdl**      | Creates the vet_* tables, whose schema now always includes `vrs_allele_ids` (empty with VRS off), and ref_ranges; the VRS tables are created in GvsAssignIds                                                                          |
| **GvsAssignIds.wdl**         | Orchestrates sample ID assignment; conditionally creates the vrs_allele_scratch and vrs_allele tables when `enable_vrs_ids=true`                                                                                                      |
| **GvsImportGenomes.wdl**     | Main ingestion engine; passes `enable_vrs_ids` to ingest tool, conditionally uploads vrs_allele_scratch parquet files to GCS, includes vrs_allele_scratch in parquet discovery, loading and verification, then runs ProcessVrsAlleles |
| **GvsBulkIngestGenomes.wdl** | High-level bulk orchestrator; accepts and forwards `enable_vrs_ids` input (default: false) to both downstream workflows                                                                                                               |

When `enable_vrs_ids=true`, the parquet discovery logic includes the `vrs_allele_scratch` prefix, ensuring those files are discovered, loaded into BigQuery and verified.

### Loading and verification

`vrs_allele` is a lookup table with one row per distinct allele. The per-sample association lives in `vet_*.vrs_allele_ids`. Ingest does not write `vrs_allele` directly. It follows the same pattern as `vcf_header_lines_scratch`:

1. Each sample's allele records go to `vrs_allele_scratch`, with its `sample_id`. The same allele appears once per sample that carries it.
2. `VerifyParquetLoading` treats `vrs_allele_scratch` like any other per-sample table. The footer row count comparison checks each sample's rows against its Parquet, and the cross-family check requires a VRS file for every sample with vet data.
3. `ProcessVrsAlleles` (`process_vrs_alleles.py`) then works on the samples in the scratch table:
   - It fails if any `vrs_allele_id` has rows that disagree on location or state, before merging anything.
   - It inserts the distinct alleles not already in `vrs_allele`. This is idempotent, so a re-run inserts nothing new.
   - It fails if any `vrs_allele_ids` entry in those samples' vet rows has no `vrs_allele` row, or if `vrs_allele` holds an id more than once. The empty-string placeholder for a `*` ALT is exempt.
   - Only then does it delete those samples' scratch rows, so a failed check leaves them for inspection.
4. `SetIsLoadedColumn` waits for `ProcessVrsAlleles`.

Like `ProcessVCFHeaders`, this assumes one ingest at a time per dataset. The steps are scoped to the samples in the scratch table, and the duplicate check catches two concurrent merges inserting the same allele.

---

## How to Use

### Enable VRS ID Computation

When invoking `CreateVariantIngestFiles` or the bulk ingest workflow:

```bash
# Java tool
java ... CreateVariantIngestFiles --enable-vrs-ids true ...

# WDL (Terra)
# Pass enable_vrs_ids=true to GvsBulkIngestGenomes workflow input
```

Default is `false` for backward compatibility; existing ingests unaffected.

### Output Files

When enabled, each sample produces:

1. **vet_*.parquet**: the standard vet file with one new column, `vrs_allele_ids`: REPEATED STRING, one ID per entry in `alt` (a comma-separated string). A spanning-deletion `*` ALT, which VRS cannot represent, gets an empty string (`VetCreator.NO_VRS_ALLELE_ID`), so the two stay aligned.

2. **vrs_allele_scratch_*.parquet**: the sample's allele records containing:
   - `sample_id`
   - `vrs_allele_id`: standardized identifier (ga4gh:VA.*)
   - `vrs_location_id`: location identifier (ga4gh:SL.*)
   - Normalized coordinates and state in VRS 0-based interbase format

These load into BigQuery as:
   - `vet_*`, with the `vrs_allele_ids` column
   - `vrs_allele_scratch`, emptied after each ingest (see Loading and verification)
   - `vrs_allele` (shared across all samples), one row per distinct allele

### Handling Heterozygous Variants

For het-vars with multiple alternate alleles (e.g., GT=1/2 with alts ["T", "ATG"]):
- Both alleles receive independent, distinct VRS IDs
- Array order matches alt field order; preserves genotype semantics
- Single row in vet_* table with arrays holding one ID per allele
- Location information stored in canonical vrs_allele table (no duplication)

---

## Architecture Notes

**Opt-In Semantics**: Default `false` preserves existing behavior; enabling VRS IDs incurs:
- CPU cost for normalization and SHA-512 digest computation
- Storage overhead: one additional column in the vet tables, one additional parquet file per sample
- BigQuery load time: one additional table to load (vrs_allele_scratch), plus the ProcessVrsAlleles merge

**Deduplication**: Allele records are duplicated across samples in `vrs_allele_scratch`. `ProcessVrsAlleles` deduplicates them on `vrs_allele_id` when merging into `vrs_allele`.

**Coordinate Systems**: Input is 1-based inclusive (VCF/gVCF standard); VRS normalization and IDs use 0-based interbase coordinates per VRS 2.0.1 specification.

---

## External References

- [VRS 2.0.1 Specification](https://vrs.ga4gh.org/)
- [RFC 8785 (JCS - JSON Canonicalization Scheme)](https://tools.ietf.org/html/rfc8785)
- [GA4GH Identifiers](https://github.com/ga4gh/ga4gh-identifiers)

