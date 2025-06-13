# Overview of the Genomic Variant Store (GVS)

The Genomic Variant Store (GVS) is a scalable cloud-based platform designed to
efficiently store, manage, and process large-scale genomic variant data. Built
as an extension of the Genome Analysis Toolkit (GATK), GVS enables joint calling
across thousands of samples and provides fast cohort extraction capabilities for
population genomics studies.

For a very high level overview see `gvs-product-sheet.pdf`.

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

### GVS Beta

GVS Beta is intended for end users to be able to run their own joint calling
using the `GvsJointVariantCalling.wdl` as a top-level orchestrating workflow.
See all the documentation in the `beta_docs` directory for more information on
GVS Beta usage.

### GVS for AoU

The All of Us Research Program implementation represents a production-scale
deployment optimized for the specific requirements of this national precision
medicine initiative.

# High-Level Architecture (WDLs, BigQuery, Terra)

## Architecture Components

### 1. Terra Platform

- **Workflow Orchestration**: All GVS operations run as WDL workflows on Terra
- **Compute Resources**: Utilizes Google Cloud Platform for scalable compute
- **Data Management**: Secure storage and access control for genomic data

### 2. Google BigQuery Backend

- **Primary Storage**: All variant data stored in BigQuery tables
- **Query Performance**: Optimized for fast analytical queries across large
  datasets
- **Schema Design**: Purpose-built tables for variants, samples, and annotations
- **Cost Tracking**: Built-in cost observability for operations

### 3. WDL Workflow Architecture

#### Core Workflows

- **`GvsJointVariantCalling.wdl`**: Master workflow orchestrating the complete
  joint calling pipeline
- **`GvsBulkIngestGenomes.wdl`**: Bulk ingestion of GVCF files into the variant
  store
- **`GvsCreateFilterSet.wdl`**: Variant quality score recalibration (VQSR) and
  filtering
- **`GvsExtractCallset.wdl`**: Cohort extraction and VCF generation

#### Supporting Workflows

- **`GvsAssignIds.wdl`**: Sample ID assignment and validation
- **`GvsImportGenomes.wdl`**: Individual genome import operations
- **`GvsPopulateAltAllele.wdl`**: Alternative allele frequency calculations
- **`GvsCreateVDS.wdl`**: Hail Variant Dataset creation for advanced analytics

#### Utility Workflows

- **`GvsCallsetCost.wdl`**: Cost analysis and tracking
- **`GvsWithdrawSamples.wdl`**: Sample withdrawal for data governance
- **`GvsValidateVDS.wdl`**: Data validation and quality checks

### 4. Data Flow Pipeline

```
GVCF Files → Sample Assignment → Bulk Ingestion → Joint Calling →
Quality Filtering → Cohort Extraction → Output VCFs/Datasets
```

### 5. BigQuery Schema Design

### Core concepts

GVS leverages BigQuery partitioning and clustering to optimize query performance
and cost. The core tables are either partitioned by sample ID or location,
depending on the access pattern for which they are optimized. GVS uses the term
"location" in a very specific way: the chromosome number is multiplied by
1,000,000,000,000 and then the position on the chromosome is added. For example,
chromosome 1 position 1000 becomes `1000000000000 + 1000 = 1000000001000`. GVS
locations use a chromosome number of 23 for the X chromosome and 24 for the Y
chromosome. This location encoding allows for efficient storage and querying in
BigQuery using a single partitionable value.

#### Core Tables

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
  4000 partition limit must be many of these tables (named `vet_001`, `vet_002`,
  etc) in a GVS dataset. Schema information for vet tables can be found in the
  `CreateBQTables.wdl` workflow in the variable `vet_schema_json`.
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

#### Annotation Tables

- **Variant Annotation Tables (VAT)**: External annotation integration
- **Gene annotation tables**: Gene-level metadata and OMIM information

# GVS "Beta" Usage

The GVS Beta represents the general-purpose implementation suitable for research
cohorts and pilot studies.

## Typical Beta Workflow

### 1. Initial Setup

```bash
# Create BigQuery dataset and tables
GvsCreateTables.wdl

# Assign sample IDs
GvsAssignIds.wdl
```

### 2. Data Ingestion

```bash
# Bulk ingest GVCF files
GvsBulkIngestGenomes.wdl
  - Input: Sample manifest with GVCF paths
  - Output: Populated BigQuery tables
```

### 3. Joint Calling Pipeline

```bash
# Complete joint calling workflow
GvsJointVariantCalling.wdl
  - Bulk ingestion
  - Alt allele population
  - VQSR filter creation
  - Cohort preparation
```

### 4. Cohort Extraction

```bash
# Extract specific cohorts
GvsExtractCallset.wdl
  - Input: Sample lists and genomic intervals
  - Output: VCF files with joint-called variants
```

## Beta Configuration

### Resource Requirements

- **Compute**: Standard Google Cloud VMs (n1-standard-4 to n1-highmem-16)
- **Storage**: Regional persistent disks
- **BigQuery**: Standard pricing tier

### Typical Use Cases

- Research cohorts (hundreds to thousands of samples)
- Method development and validation
- Pilot studies for larger production deployments

# GVS AoU (All of Us) Usage

The All of Us Research Program implementation represents a production-scale
deployment optimized for the specific requirements of this national precision
medicine initiative.

## AoU-Specific Features

### 1. Scale Optimization

- **Sample Volume**: Designed for 100,000+ participant genomes
- **Compute Resources**: Enhanced machine types and parallel processing
- **Storage**: Premium storage tiers for performance

### 2. Data Governance

- **Sample Withdrawal**: Automated participant withdrawal capabilities
- **Access Control**: Integration with AoU data access policies
- **Audit Logging**: Comprehensive operation tracking

### 3. Specialized Workflows

#### Precision and Sensitivity Analysis

- **`AoU_PRECISION_SENSITIVITY.md`**: Detailed callset validation procedures
- **GIAB Integration**: Genome in a Bottle truth set comparisons
- **Quality Metrics**: Population-specific QC thresholds

#### Production Pipeline

```bash
# AoU-optimized joint calling
GvsJointVariantCalling.wdl (AoU configuration)
  - Enhanced resource allocation
  - Population-specific filtering
  - Comprehensive quality control
```

### 4. Performance Characteristics

#### Throughput

- **Ingestion**: ~1,000 samples per day capacity
- **Joint Calling**: Full cohort processing in 2-3 days
- **Extraction**: Sub-hour turnaround for most queries

#### Cost Management

- **Resource Optimization**: Preemptible instances where appropriate
- **Storage Tiering**: Automated data lifecycle management
- **Query Optimization**: Partitioned tables and clustering

## AoU Deliverables

### Primary Outputs

1. **Joint-Called VCF**: Population-scale variant calls
2. **Variant Annotation Tables**: Comprehensive functional annotations
3. **Quality Metrics**: Per-sample and per-variant QC statistics
4. **Hail Variant Datasets**: For advanced population genetics analyses

### Integration Points

- **Researcher Workbench**: Direct query access for approved researchers
- **Data Browser**: Web interface for variant exploration
- **Analysis Tools**: Integration with population genetics software

## Configuration Differences

| Aspect          | Beta           | AoU Production         |
| --------------- | -------------- | ---------------------- |
| Sample Scale    | 100s-1000s     | 100,000+               |
| Compute Tier    | Standard       | High-performance       |
| Storage         | Regional       | Multi-regional premium |
| Monitoring      | Basic          | Comprehensive          |
| Data Governance | Research-grade | Production-grade       |
| SLA             | Best effort    | 99.9% availability     |

## Key Scripts and Utilities

### Cost Analysis

- **`workflow_compute_costs.py`**: Workflow cost calculation and tracking
- **`GvsCallsetCost.wdl`**: Automated cost reporting

### Data Validation

- **`vds_validation.py`**: Variant Dataset validation
- **`compare_data.py`**: Cross-platform data comparison

### Monitoring and Observability

- **`summarize_task_monitor_logs.py`**: Task performance analysis
- **Cost observability tables**: Real-time cost tracking in BigQuery

---

_This document provides an overview of the GVS architecture and usage patterns.
For detailed implementation guidance, refer to the individual WDL files and the
Terra quickstart documentation._
