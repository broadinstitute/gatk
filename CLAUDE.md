# Overview of the Genomic Variant Store (GVS)

The Genomic Variant Store (GVS) is a scalable cloud-based platform designed to
efficiently store, manage, and process large-scale genomic variant data. GVS
enables joint calling across hundreds of thousands of samples and provides fast
cohort extraction capabilities for population genomics studies. For a high level
overview see `gvs-product-sheet.pdf`.

GVS was created as a branch of the Genome Analysis Toolkit (GATK), though years
later it is not clear if or when GVS will be merged back into the main GATK
branch. The principal GVS branch in git is called `ah_var_store`. Although GVS
does build upon some of the GATK codebase, the majority of the code in this
repository is not related to GVS.

## Key Features

- **Scalable Storage**: Leverages Google BigQuery for petabyte-scale genomic
  data storage with high availability
- **Joint Calling**: Performs population-scale joint genotyping across large
  cohorts
- **Cohort Extraction**: Rapid extraction of variant subsets for specific sample
  cohorts
- **Cloud-Native**: Designed for Terra platform with WDL-based workflows

## Usage Scenarios

GVS supports two primary usage scenarios as described below:

### GVS Beta Callsets

GVS Beta is intended for end users to be able to run their own joint calling
using the `GvsJointVariantCalling.wdl` as a top-level workflow. See the
documentation in the `beta_docs` directory for more information on GVS Beta
usage.

### AoU Callsets

GVS is used to produce joint callsets for the All of Us Research Program at the
scale of 400,000+ samples. These callsets push the scale limits of GVS and are
not suited for the `GvsJointVariantCalling.wdl` workflow. Instead, more narrowly
scoped GVS workflows are run by the engineers of the Variants team following an
AoU-specific protocol as described in `AOU_DELIVERABLES.md`. Rather than the
collection of VCF files that is produced from GVS Beta, the core deliverable of
an AoU callset is a Hail VDS (Variant Dataset). More AoU-specific documentation
can be found under `scripts/variantstore/docs/aou`.

Conceptually an AoU callset is similar to a GVS Beta callset up to the point of
extraction. While a Beta callset would run `GvsExtractCallset.wdl` to extract a
set of VCF files, an AoU callset will run `GvsExtractAvroFilesForHail.wdl` and
then `GvsCreateVDS.wdl` to generate a Hail VDS. The data within the VDS is used
to generate an ancestry file for creating the VAT (Variant Annotations Table)
and interval lists for the "small callsets" described later in this document.

## Key Concepts

### Hard filtered versus soft filtered variants

During the GVS joint calling process variants are evaluated for quality at both
the site and allele level. By default GVS evaluates variant quality scores using
a machine learning model called VETS (Variant Extract-Train-Score). GVS formerly
used a machine learning model called VQSR (Variant Quality Score Recalibration)
for this purpose. While still supported in GVS, VQSR is rarely used in making
GVS callsets today.

GVS stores all variants and their associated quality scores in its BigQuery
schema, regardless of whether these variants pass or fail quality filters. GVS
artifacts (VCF, PGEN, VDS, etc), are classified as being either "hard filtered"
or "soft filtered" depending on how they treat variants with failing quality
scores. "Hard filtered" artifacts completely exclude variants with failing
quality scores, while "soft filtered" artifacts include all variants along with
their quality scores, leaving it to the user's discretion how to handle these
variants in downstream analysis.

# High-Level Architecture (WDLs, BigQuery, Terra)

## BigQuery locations

GVS leverages BigQuery partitioning and clustering to optimize query performance
and cost. The core tables are either partitioned by sample ID or location,
depending on the access pattern for which they are optimized. GVS uses the term
"location" in a very specific way: the chromosome number is multiplied by
1,000,000,000,000 and then the position on the chromosome is added. For example,
chromosome 1 position 1000 becomes `1000000000000 + 1000 = 1000000001000`. GVS
locations use a chromosome number of 23 for the X chromosome and 24 for the Y
chromosome. This location encoding allows for efficient storage and querying in
BigQuery using a single partitionable value. References in GVS code to
"position" usually refer to position on a given chromosome and do not mean the
same thing as "location".

## Core BigQuery Tables

- **`sample_info`**: Holds sample name, control / non-control status, and a
  withdrawal date if applicable.
- **`sample_load_status`**: This table tracks which data has been loaded for
  each sample. As described in the LoadStatusValue enum in LoadStatus.java, the
  values currently used include REFERENCES_LOADED, VARIANTS_LOADED,
  HEADERS_LOADED. For historical reasons some legacy values are also recognized.
- **`sample_chromosome_ploidy`**: This table contains a per-sample,
  per-chromosome of the ploidy observed in its input GVCF file.
- **`vet_%`**: Primary variant information (sample id, location, alleles,
  annotations). These tables are partitioned by sample id. Due to BigQuery's
  4000 partition limit there can be many `vet_%` tables in a GVS dataset, named
  `vet_001`, `vet_002`, etc. Schema information for vet tables can be found in
  the `CreateBQTables.wdl` workflow in the variable `vet_schema_json`.
- **`ref_ranges_%`**: Reference genome locations, lengths, and qualities.
  Similar to vet, this table is partitioned by sample id and can hold at most
  4000 samples, so there will be multiple `ref_ranges_%` for larger datasets.
  Note that `ref_ranges_%` tables have two alternate schemas, compressed and
  uncompressed, to accommodate different data storage needs. Both schemas are
  described in the `CreateBQTables.wdl` workflow.
- **`vcf_header_lines`**: Contains all unique VCF header seen while ingesting
  all samples.
- **`sample_vcf_header`**: Contains header information on a per-sample basis,
  referencing enries in the `vcf_header_lines` table.
- **`vcf_header_lines_scratch`**: Scratch table for temporary storage of VCF
  header lines during ingestion. This table is used to avoid duplicate entries
  in the `vcf_header_lines` table.
- **`alt_allele`**: This table contains allele information similar to that in
  the `vet_%` tables, but is partitioned by location to support creation of a
  filter model in joint calling. The `alt_allele` table also splits
  multi-allelic (aka "hetvar" or "1/2") genotypes into separate rows for each
  allele.
- **`filter_set_sites`**: Contains filtering information at a site level only,
  no alleles.
- **`filter_set_info`**: Contains filtering information at a position and allele
  level, including the VQSR model and other metadata.
- **`cost_observability`**: Operation cost tracking, written to by the various
  tools that comprise the GVS pipeline.

## Variant Annotation Table (VAT)

### Overview

The Variant Annotation Table (VAT) is an artifact that is currently only
delivered as part of AoU callsets. See the documentation in
`scripts/variantstore/variant-annotations-table/README.md` for more information
on the VAT.

### VIDs

VIDs, or “Variant IDs” are the primary key of the VAT table. Here are three
example VIDs:

1. 1-668638-G-GA
2. 16-88461103-AGGCCCTGGTGAGGACCCAGGGATGGGGCGGGAGACCTGGTGAGGACCGAGG-A
3. Y-56694638-G-A

A VID consists of four components separated by dashes. From left to right these
are chromosome, position, reference allele, and variant allele.

In the first example above, the variant described by the VID is on chromosome 1
(chr1 in hg38 terms), position 668638, with a reference allele of G and a
variant allele of GA. Because the variant allele is longer than the reference
allele, this variant represents an insertion.

In the second example the variant is on chromosome 16 (chr16) at position
88461103 with a reference allele of
AGGCCCTGGTGAGGACCCAGGGATGGGGCGGGAGACCTGGTGAGGACCGAGG and a variant allele of A.
Because the variant allele is shorter than the reference allele, this represents
a deletion.

In the third example, the variant is on chromosome Y (chrY) at position 56694638
with a reference allele of G and a variant allele of A. This VID represents a
SNP (single nucleotide polymorphism).

### VID to Participant ID Mapping Table.

See the documentation in `AOU_DELIVERABLES.md` for information on how to
generate the VID to Participant ID Mapping Table.

Note that for the past Echo callset there may need to be an additional step to
patch this table for a set of VIDs that did not have corresponding Participant
IDs. See the directory `pseudo_vids_only_in_vat` for more information on
unmatched VIDs what were discovered in the VATs of the Delta and Echo callsets.
