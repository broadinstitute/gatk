# Overview of the Genomic Variant Store (GVS)

The Genomic Variant Store (GVS) is a scalable cloud-based platform designed to
efficiently store, manage, and process large-scale genomic variant data. GVS
enables joint calling across hundreds of thousands of samples and provides fast
cohort extraction capabilities for population genomics studies. For a high level
overview see `gvs-product-sheet.pdf`.

GVS was created as a branch of the Genome Analysis Toolkit (GATK). The principal
GVS branch in git is called `ah_var_store`. While this GVS branch has been
updated from GATK `master` periodically, it is not clear if GVS will ever be
merged back into the main GATK branch. While GVS does build upon some of the
GATK codebase, the majority of the code in the GATK repository (even on the
`ah_var_store` branch) is not related to GVS.

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
scale of 400,000+ samples. AoU callsets push the scale limits of GVS and are not
suited for the `GvsJointVariantCalling.wdl` workflow. For AoU, more narrowly
scoped GVS workflows are run by the engineers of the Variants team following an
AoU-specific protocol as described in `AOU_DELIVERABLES.md`. Rather than the
collection of VCF files produced from GVS Beta, the core deliverable of an AoU
callset is a Hail VDS (Variant Dataset). More AoU-specific documentation can be
found under `scripts/variantstore/docs/aou`.

Conceptually an AoU callset is similar to a GVS Beta callset up to the point of
extraction. While a Beta callset would run `GvsPrepareCallset.wdl` and then
`GvsExtractCallset.wdl` to extract a set of VCF files, an AoU callset will run
`GvsExtractAvroFilesForHail.wdl` and then `GvsCreateVDS.wdl` to generate a Hail
VDS. The data within the VDS is used to generate an ancestry file for creating
the VAT (Variant Annotations Table) and interval lists for the "small callsets"
described in `AOU_DELIVERABLES.md`.

## Key Concepts

### Hard filtered versus soft filtered variants

During the GVS joint calling process variants are evaluated for quality at both
the site and allele level. By default GVS evaluates variant quality scores using
a machine learning model called VETS (Variant Extract-Train-Score). GVS formerly
used a machine learning model called VQSR (Variant Quality Score Recalibration)
for this purpose. While still supported in GVS, VQSR is rarely used in making
GVS callsets today. See
`scripts/variantstore/docs/release_notes/VETS_Release.pdf` for more information
on VETS and the transition away from VQSR.

GVS stores all variants and their associated quality scores in its BigQuery
schema, regardless of whether these variants pass or fail quality filters. GVS
artifacts (VCF, PGEN, VDS, etc), are classified as being either "hard filtered"
or "soft filtered" depending on how they treat variants with failing quality
scores. "Hard filtered" artifacts completely exclude variants with failing
quality scores, while "soft filtered" artifacts include all variants along with
their quality scores, leaving it to the user's discretion how to handle these
variants in downstream analysis.

# High-Level Architecture (WDLs, BigQuery, Terra)

## Java Artifacts

GVS is written in a mixture of Java, Python, WDL (Workflow Description
Language), and shell scripts. The Java code can be built from the repository
root using the `./gradlew clean shadowJar` command. This compilation requires a
JDK of at least version 17. For local compilation, if
[SDKMAN!](https://sdkman.io/) is available, one of the Java 17 versions can be
used (most GATK builds seem to use Temurin).

## Docker Versioning

The `GvsUtils.wdl` workflow's `GetToolVersions` task is used to retrieve the
version of various Docker images used by GVS. This task is usually run at the
start of a GVS workflow to ensure that the correct versions of tools are being
used.

## Reference Versioning

The `GvsUtils.wdl` workflow's `GetReference` task is used to retrieve the
canonical paths to references used by GVS. This task is usually run at the start
of a GVS workflow to ensure that the correct references are being used.

## Docker artifacts

GVS uses multiple Docker images to run its workflows. Some of these are specific
to GVS as documented below.

### GATK Docker Image

This builds a GATK Docker image based on the GVS branch of the GATK codebase.
Instructions for building this image can be found in the
`scripts/variantstore/docs/Build Docker from VM/building_gatk_docker_on_a_vm.md`
file. An x86 cloud-based VM is used to build this x86 Docker image, which is
then pushed to the Broad's GAR (Google Artifact Registry) repository.

### Variants Docker Image

This is a Google Cloud SDK Alpine-based Docker image to which many useful tools,
Python packages, and scripts have been added. This can be built on a developer's
local machine using the `scripts/variantstore/scripts/build_docker.sh` script.

### Nirvana Docker Image

This is an Illumina Nirvana-based Docker image build from
`scripts/variantstore/variant-annotations-table/build_docker.sh`. Since Nirvana
is no longer being updated this image hopefully will also not need to be updated
in the future.

## "location" versus "position"

GVS leverages BigQuery partitioning and clustering to optimize query performance
and cost. The core tables are either partitioned by sample ID or location,
depending on the access pattern for which they are optimized. GVS uses the term
"location" in a very specific way: the chromosome number is multiplied by
1,000,000,000,000 and then the position on the chromosome is added. For example,
chromosome 1 position 1000 becomes `1000000000000 + 1000 = 1000000001000`. GVS
locations use a chromosome number of 23 for the X chromosome and 24 for the Y
chromosome. This location encoding allows for efficient storage and querying in
BigQuery using a single partitionable value. GVS references to "position"
usually refer to position on a given chromosome and do not mean the same thing
as "location".

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
  referencing entries in the `vcf_header_lines` table.
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

## Key Workflows

### GvsJointVariantCalling.wdl

This is the top-level workflow for GVS Beta joint calling. It orchestrates the
joint calling process, including sample ingestion, variant calling, and
filtering. The workflow is designed to be run on the Terra platform and can be
customized for specific datasets. This workflow is _not_ used for AoU callsets,
which use a different set of workflows as described in `AOU_DELIVERABLES.md`.

### GvsBulkIngest.wdl

This workflow is used to ingest large sets of samples into GVS and is
extensively documented in `gvs-bulk-ingest-details.md`. This workflow is called
by `GvsJointVariantCalling.wdl` for Beta callsets and is called directly for AoU
callsets. One significant change for the Foxtrot callset is that this workflow
will now be directly supplied with a "FOFN" (File of File Names) TSV rather than
reading from Terra data tables. This change was motivated by hitting the scale
limits of Terra data tables and being unable to read AoU-sized sample sets, but
providing this workflow with a FOFN TSV directly is also much simpler than the
previous TSV to Terra data table to TSV approach.

### GvsPopulateAltAllele.wdl

This workflow populates the `alt_allele` table with allele information from the
`vet_%` tables. This table will enable the creation of a filter model in joint
calling in subsequent steps. The `alt_allele` table is partitioned by location
with split multi-allelics to support efficient querying and filtering.

### GvsCreateFilterSet.wdl

This workflow reads the `alt_allele` tables, scattering widely over locations to
generate a filter model and populate the `filter_set_sites` and
`filter_set_info` tables.

### GvsPrepareCallset.wdl

This workflow prepares a GVS callset for VCF or PGEN extraction. Sample,
reference, and variant data is copied to a set of "prepare tables" that are read
by GVS extract tools.

### GvsExtractCallset.wdl / GvsExtractCallsetPgenMerged.wdl

These workflows extract a GVS callset into VCF or PGEN files. They read the
prepare tables created by `GvsPrepareCallset.wdl` and generate the requested
output files. Output files are partitioned at a minimum by chromosome, but for
"wide" extracts (callsets with many large numbers of samples), there may be
hundreds of output files per chromosome.

## Variant Annotation Table (VAT)

### Overview

The Variant Annotation Table (VAT) is an artifact that is currently only
delivered as part of AoU callsets. See the documentation in
`scripts/variantstore/variant-annotations-table/README.md` for more information
on the VAT. The core workflow for generating the VAT is
`GvsCreateVATFromVDS.wdl`, which is run after the AoU callset has been created.

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
