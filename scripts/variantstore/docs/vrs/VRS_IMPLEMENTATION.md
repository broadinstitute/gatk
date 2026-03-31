# VRS Allele ID Implementation in GATK Variant Store

## Overview

This document describes the VRS (Variant Representation Specification) allele ID feature integrated into GATK's variant store ingest pipeline. VRS IDs provide standardized, globally-unique identifiers for genomic alleles, enabling interoperability across bioinformatics tools and databases.

**Current Status**: ✅ Fully implemented and integrated end-to-end, from Java ingest tool through WDL orchestration to BigQuery loading.

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

| Class | Purpose |
|-------|---------|
| **CreateVariantIngestFiles** | Entry point for variant ingest; accepts `--enable-vrs-ids` flag (default: false); instantiates VrsAlleleCreator when enabled |
| **VetCreator** | Extracts variants and computes VRS IDs via `VrsNormalizer` and `VrsIdComputer` during `apply()`; stores IDs in VET table and delegates allele records to `VrsAlleleCreator` |
| **VrsAlleleCreator** | Orchestrates writing of VrsAlleleRecord objects to vrs_allele output (TSV, Parquet, or BigQuery); one instance per sample |
| **VrsNormalizer** | Implements VRS 2.0.1 normalization steps (trimming, rolling, expansion, state resolution) to produce canonical allele representation |
| **VrsIdComputer** | Generates ga4gh:VA.* and ga4gh:SL.* identifiers from normalized alleles using RFC 8785 canonical JSON + SHA-512 digest |
| **VrsAlleleRecord** | POJO transport object holding `vrs_allele_id` and `vrs_location_id` for each alternate allele; supports het-vars via array storage |

### Test Classes

- `VrsAlleleCreatorUnitTest` (12 tests): File naming, TSV/Parquet I/O, no-op scenarios
- `VrsNormalizerTest` (7 tests): Trimming, rolling, expansion, state resolution for SNPs, indels, substitutions
- `VrsIdComputerTest` (17 tests): Deterministic ID generation, format validation, coordinate/sequence variation
- **Integration test** (in CreateVariantIngestFilesIntegrationTest): End-to-end VRS ID computation from GVCFs

**Current test status**: 133 unit tests passing, integration test included

---

## WDL Integration

VRS ID flow is wired through three WDL workflows with consistent opt-in semantics:

| Workflow | Role |
|----------|------|
| **GvsCreateTables.wdl** | Creates mandatory BigQuery tables (vet_*, ref_ranges); does NOT include vrs_allele (optional) |
| **GvsAssignIds.wdl** | Orchestrates sample ID assignment; conditionally creates optional vrs_allele table when `enable_vrs_ids=true` |
| **GvsImportGenomes.wdl** | Main ingestion engine; passes `enable_vrs_ids` to ingest tool, conditionally uploads vrs_allele parquet files to GCS, includes vrs_allele in parquet discovery for BigQuery loading |
| **GvsBulkIngestGenomes.wdl** | High-level bulk orchestrator; accepts and forwards `enable_vrs_ids` input (default: false) to both downstream workflows |

When `enable_vrs_ids=true`, the parquet discovery logic includes the `vrs_allele` prefix, ensuring those files are discovered and loaded into BigQuery.

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

When enabled, three output formats are produced (per sample):

1. **vet_*.tsv** (or parquet): Standard variant extraction with two new columns:
   - `vrs_allele_ids`: REPEATED STRING, one ID per alternate allele
   - `vrs_location_ids`: REPEATED STRING, one location ID per alternate allele

2. **vrs_allele_*.parquet**: Canonical allele records containing:
   - `vrs_allele_id`: Standardized identifier (ga4gh:VA.*)
   - `vrs_location_id`: Location identifier (ga4gh:SL.*)
   - Normalized coordinates and sequence in VRS 0-based interbase format

3. **BigQuery tables**:
   - `vet_*` with vrs_allele_ids, vrs_location_ids columns
   - `vrs_allele` (shared across all samples) for allele deduplication and canonical storage

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
- Storage overhead: two additional columns in VET table, one additional parquet file per sample
- BigQuery load time: one additional table to load (vrs_allele)

**Deferred Deduplication**: Allele records may be duplicated across samples (same allele in different samples). Deduplication is handled via BigQuery view or post-load MERGE on vrs_allele_id, enabling efficient discovery without per-sample uniqueness constraints.

**Coordinate Systems**: Input is 1-based inclusive (VCF/gVCF standard); VRS normalization and IDs use 0-based interbase coordinates per VRS 2.0.1 specification.

---

## External References

- [VRS 2.0.1 Specification](https://vrs.ga4gh.org/)
- [RFC 8785 (JCS - JSON Canonicalization Scheme)](https://tools.ietf.org/html/rfc8785)
- [GA4GH Identifiers](https://github.com/ga4gh/ga4gh-identifiers)

