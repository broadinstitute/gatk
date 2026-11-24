# VRS Allele ID Implementation in GATK Variant Store

## Overview

This document describes the VRS (Variant Representation Specification) allele ID feature integrated into GATK's variant store ingest pipeline. VRS IDs provide standardized, globally-unique identifiers for genomic alleles, enabling interoperability across bioinformatics tools and databases.

**Current status**: Implemented from the Java ingest tool through WDL orchestration to BigQuery loading, behind `enable_vrs_ids`. A dataset either carries VRS IDs or does not, and an ingest follows the dataset (see Enable VRS ID Computation). The IDs have not yet been checked against vrs-python beyond a single normalizer case (VS-2059). `GvsPopulateAltAllele` copies them into `alt_allele` in a dataset that carries them (VS-2061).

---

## Feature Description

The VRS ID feature computes standardized identifiers (ga4gh:VA.* format) for each variant allele during ingestion. When enabled, the pipeline:

1. **Normalizes** each allele according to the VRS 2.0.1 specification: fully-justified normalization (trim, roll, expand across the whole ambiguous region) and state resolution. This is not left-alignment: an indel in a repeat gets an interval spanning the whole repeat, in 0-based interbase coordinates with no VCF anchor base, so it does not line up with the single position of a left-aligned vid. Multi-allelic sites are not decomposed: each ALT is normalized on its own. Only a pure insertion or deletion (one side empty after trimming) is rolled. A substitution, and a delins where both trimmed sides are non-empty but differ in length (`GACA`→`GT`, which trims to `ACA`→`T`), stays at its trimmed interval with a literal state, as in vrs-python. For an insertion, GVS decides whether the expanded allele is reference-derived, and its `repeatSubunitLength`, with its own search (`findReferenceDerivedRepeatSubunitLength`) rather than vrs-python's fixed-cycle-start check. We kept it because no case has been found where the two disagree: it matched vrs-python 2.3.3 on 20,000 random insertions into repetitive sequence, and on every repeat case in the VS-2059 concordance test. Both give `repeatSubunitLength` as the inserted length (4 for two inserted CA units), not the repeat's period. A case where they differ belongs in that test's fixture.
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
- `VrsNormalizerTest`: trimming, rolling, expansion and state resolution for SNPs, indels, substitutions and delins
- `VrsIdComputerTest`: deterministic ID generation and format; it does not compare IDs with known values (VS-2059)
- `VrsAlleleRecordTest`, `RefGetAccessionProviderTest`
- `VetCreatorUnitTest`: a spanning-deletion `*` ALT, and every ALT on a contig with no refget accession, gets the placeholder and keeps `vrs_allele_ids` aligned with `alt`
- `CreateVariantIngestFilesIntegrationTest`: end-to-end VRS ID computation from a gVCF, against the checked-in chr20/21 subset of hg38.

---

## WDL Integration

VRS ID flow is wired through five WDL workflows, with the decision made once per ingest by `ResolveEnableVrsIds` (`GvsUtils.wdl`):

| Workflow                     | Role                                                                                                                                                                                                                                                                                              |
|------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **GvsCreateTables.wdl**      | Creates the vet_* tables, with `vrs_allele_ids` only when VRS IDs are on, and ref_ranges; the VRS tables are created in GvsAssignIds                                                                                                                                                              |
| **GvsAssignIds.wdl**         | Orchestrates sample ID assignment; resolves `enable_vrs_ids` against the dataset before assigning any sample ID, failing on a contradiction, then creates the vet tables to match and, with VRS on, the vrs_allele_scratch and vrs_allele tables. Outputs the resolved value                      |
| **GvsImportGenomes.wdl**     | Main ingestion engine; resolves `enable_vrs_ids` again before loading, so a direct run is checked too, and passes the result to the ingest tool. With VRS on, uploads vrs_allele_scratch parquet files, includes them in parquet discovery, loading and verification, then runs ProcessVrsAlleles |
| **GvsBulkIngestGenomes.wdl** | High-level bulk orchestrator; passes an optional `enable_vrs_ids` to GvsAssignIds, and the value GvsAssignIds resolved to the data pass of GvsImportGenomes                                                                                                                                       |
| **GvsPopulateAltAllele.wdl** | Reads the dataset's VRS state with `ResolveEnableVrsIds` (flag unset), creates `alt_allele` with `vrs_allele_id` only when the vet tables have `vrs_allele_ids`, fails if an existing `alt_allele` disagrees, and fills the column (see alt_allele below)                                         |

With VRS IDs on, the parquet discovery logic includes the `vrs_allele_scratch` prefix, ensuring those files are discovered, loaded into BigQuery and verified.

### Loading and verification

`vrs_allele` is a lookup table with one row per distinct allele. The per-sample association lives in `vet_*.vrs_allele_ids`. Ingest does not write `vrs_allele` directly. It follows the same pattern as `vcf_header_lines_scratch`:

1. Each sample's allele records go to `vrs_allele_scratch`, with its `sample_id`. The same allele appears once per sample that carries it.
2. `VerifyParquetLoading` treats `vrs_allele_scratch` like any other per-sample table. The footer row count comparison checks each sample's rows against its Parquet, and the cross-family check requires a VRS file for every sample with vet data.
3. `ProcessVrsAlleles` (`process_vrs_alleles.py`) then works on the samples in the scratch table:
   - It fails if any `vrs_allele_id` has rows that disagree on location or state, before merging anything.
   - It inserts the distinct alleles not already in `vrs_allele`. This is idempotent, so a re-run inserts nothing new.
   - It fails if any `vrs_allele_ids` entry in those samples' vet rows has no `vrs_allele` row, or if `vrs_allele` holds an id more than once. The empty-string placeholder (see Output Files) is exempt.
   - Only then does it delete those samples' scratch rows, so a failed check leaves them for inspection.
4. `SetIsLoadedColumn` waits for `ProcessVrsAlleles`.

Like `ProcessVCFHeaders`, this assumes one ingest at a time per dataset. The steps are scoped to the samples in the scratch table, and the duplicate check catches two concurrent merges inserting the same allele.

If two merges do race, a rerun keeps failing the duplicate check, because the insert never removes a row. The duplicate rows are identical, since each id is a hash of its row's content and the conflict check runs first, so recover by keeping one of each, then rerun `ProcessVrsAlleles` to finish the merge. Rewrite the table in place, as below, rather than with `CREATE OR REPLACE TABLE ... AS SELECT`, which turns the `REQUIRED` columns `NULLABLE`:

```sql
CREATE TEMP TABLE deduped AS SELECT DISTINCT * FROM `<project>.<dataset>.vrs_allele`;
TRUNCATE TABLE `<project>.<dataset>.vrs_allele`;
INSERT INTO `<project>.<dataset>.vrs_allele` SELECT * FROM deduped;
```

---

## How to Use

### Enable VRS ID Computation

When invoking `CreateVariantIngestFiles` or the bulk ingest workflow:

```bash
# Java tool
java ... CreateVariantIngestFiles --enable-vrs-ids true -R Homo_sapiens_assembly38.fasta --ref-version 38 ...

# WDL (Terra)
# Set enable_vrs_ids=true on GvsBulkIngestGenomes to start a new dataset with VRS IDs
```

The Java tool fails at startup if `--enable-vrs-ids` is set without a reference, or without `--ref-version 38`: the refget accessions are GRCh38 only, as is GVS. The WDL always passes both when the flag is on.

Whether an ingest computes VRS IDs follows the dataset. Its vet tables record it: an ingest with VRS IDs on creates them with the `vrs_allele_ids` column, and one with them off creates them without it. So no dataset has a column of empty arrays, and every sample in a dataset has IDs or none does. `ResolveEnableVrsIds` reads the vet tables and decides:

| Dataset                            | `enable_vrs_ids` unset | `true`          | `false` |
|------------------------------------|------------------------|-----------------|---------|
| New: no vet tables yet             | off                    | on              | off     |
| Vet tables with `vrs_allele_ids`   | on                     | on              | fails   |
| Vet tables without it              | off                    | fails (VS-2064) | off     |
| Some vet tables with, some without | fails                  | fails           | fails   |

A contradiction fails rather than being overridden, so an operator never gets the opposite of what they asked for without noticing. The error names the dataset and the tables. Off into a VRS dataset would otherwise load silently, with empty `vrs_allele_ids`. On into a dataset without the column would fail at the Parquet load, which cannot add it. Adding the column and computing IDs for an existing dataset is a separate workflow, not yet written (VS-2064), so to get VRS IDs, ingest into a new dataset.

`GvsAssignIds` resolves before assigning any sample ID, and outputs the result, which `GvsBulkIngestGenomes` passes on rather than its own input. `GvsImportGenomes` resolves again before loading, so a direct run of it is held to the same rule.

### Output Files

When enabled, each sample produces:

1. **vet_*.parquet**: the standard vet file with one new column, `vrs_allele_ids`: REPEATED STRING, one ID per entry in `alt` (a comma-separated string). An ALT with no ID gets an empty string (`VetCreator.NO_VRS_ALLELE_ID`), so the two stay aligned. That covers a spanning-deletion `*`, which VRS cannot represent, and every ALT on a contig outside chr1–22, X, Y and M (alt, decoy, HLA), which has no refget accession. The default calling intervals cover only those standard chromosomes.

2. **vrs_allele_scratch_*.parquet**: the sample's allele records containing:
   - `sample_id`
   - `vrs_allele_id`: standardized identifier (ga4gh:VA.*)
   - `vrs_location_id`: location identifier (ga4gh:SL.*)
   - Normalized coordinates and state in VRS 0-based interbase format

These load into BigQuery as:
   - `vet_*`, with the `vrs_allele_ids` column
   - `vrs_allele_scratch`, emptied after each ingest (see Loading and verification)
   - `vrs_allele` (shared across all samples), one row per distinct allele

### alt_allele

`alt_allele` follows the vet tables, as an ingest does. In a dataset whose vet tables have `vrs_allele_ids`, `GvsPopulateAltAllele` creates `alt_allele` with a nullable STRING `vrs_allele_id` as its last column and fills it. In a dataset without them, `alt_allele` has no such column, and the query is the same as before VS-2061. `GvsPopulateAltAllele` does not take `enable_vrs_ids`: it calls `ResolveEnableVrsIds` with the flag unset, which reports the dataset's state and fails if the vet tables disagree.

`alt_allele` has one row per allele, so each row takes one entry of `vet.vrs_allele_ids`: the entry at `allele_pos - 1`, the offset `alt_allele_positions.sql` uses for `alt`. The empty-string placeholder becomes NULL, so a `*` row, or a row on a contig with no refget accession, has a NULL ID. Every non-NULL ID has a `vrs_allele` row.

An existing `alt_allele` that disagrees with the vet tables fails the workflow before anything is written: one built before the dataset had VRS IDs lacks the column, and the column would be empty if it were added. `alt_allele` is built entirely from the vet tables, so the recovery the error gives is to drop it and rerun `GvsPopulateAltAllele`. Converting an existing dataset to VRS IDs (VS-2064) has to fill `alt_allele.vrs_allele_id` as well as `vet.vrs_allele_ids`, or the dataset fails this check afterwards.

The ID does not depend on where the caller placed an indel, so it can join `alt_allele` directly to any VRS-keyed source. That could replace the windowed synonym search that maps left-aligned vids to `alt_allele` rows (`GvsMapUnmappedVIDs`), which is not part of VS-2061.

### Handling Heterozygous Variants

For het-vars with multiple alternate alleles (e.g., GT=1/2 with alts ["T", "ATG"]):
- Both alleles receive independent, distinct VRS IDs
- Array order matches alt field order; preserves genotype semantics
- Single row in vet_* table with arrays holding one ID per allele
- Location information stored in canonical vrs_allele table (no duplication)

---

## Architecture Notes

**Opt-In Semantics**: A new dataset is off unless `enable_vrs_ids` is true, and an existing one stays as it is, so existing behavior is preserved. Enabling VRS IDs incurs:
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

