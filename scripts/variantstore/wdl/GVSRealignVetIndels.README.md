# GVS Realign VET Indels Workflow

## Overview

The `GVSRealignVetIndels` workflow performs indel realignment on VET (Variant Evidence Table) data using the `VetIndelRealigner` tool. The workflow is designed to process large cohorts by scattering across individual samples within a VET table range.

## Workflow Purpose

This workflow:
1. **Realigns indels**: Uses GATK's left-alignment algorithms to standardize indel representations
2. **Processes VET data**: Reads from BigQuery VET tables containing variant evidence
3. **Scales efficiently**: Scatters processing across thousands of samples in parallel
4. **Generates statistics**: Produces detailed CSV outputs for analysis and quality control

## Input Parameters

### Required Parameters
- `dataset_name` (String): Name of the BigQuery dataset containing VET tables
- `project_id` (String): Google Cloud project ID where the BigQuery dataset is located
- `vet_table_number` (Int): VET table number to process (e.g., 7 for vet_007)

### Optional Parameters
- `reference_fasta` (String): Path to reference genome FASTA file
  - Default: `"gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"`
  - For canine data: Use UU_Cfam_GSD_1.0_ROSY reference
- `samples_per_vet_table` (Int): Number of samples per VET table (default: 4000)
- `max_records_per_sample` (Int): Maximum records to process per sample (default: 1,000,000)
- `gatk_override` (File): Optional custom GATK JAR file
- `gatk_docker` (String): GATK Docker image to use
- `git_branch_or_tag` (String): Git branch/tag for version tracking
- `git_hash` (String): Git hash for version tracking

### Resource Parameters (per task)
- `memory_gb` (Int): Memory allocation in GB (default: 8)
- `disk_gb` (Int): Disk space in GB (default: 50) 
- `cpu_count` (Int): Number of CPUs (default: 2)
- `preemptible_tries` (Int): Number of preemptible attempts (default: 3)

## Sample ID Calculation

The workflow calculates sample ID ranges using the formula:
```
sample_id_start = (vet_table_number - 1) * samples_per_vet_table + 1
sample_id_end = sample_id_start + samples_per_vet_table - 1
```

For example:
- VET table 1: samples 1-4000
- VET table 2: samples 4001-8000  
- VET table 7: samples 24001-28000

## Output Files

The workflow generates four types of CSV output files per sample:

### 1. Realigned Variants (`*_realigned_variants.csv`)
Contains details of variants that were realigned:
- Sample ID
- Original and realigned genomic locations (encoded format)
- Reference and alternate alleles
- Quality metrics

### 2. Position Bucket Histogram (`*_position_bucket_histogram.csv`) 
Genomic distribution analysis:
- Position bucket ranges
- Count of variants per bucket
- Coverage statistics

### 3. Shifted Indels Per Sample (`*_shifted_indels_per_sample.csv`)
Per-sample realignment statistics:
- Sample ID  
- Number of indels realigned
- Shift distance metrics

### 4. Indel Realignment Summary (`*_indel_realignment_summary.csv`)
Overall summary statistics:
- Total variants processed
- Realignment success rate
- Quality score distributions
- Performance metrics

### 5. Monitoring Logs (`monitoring.log`)
Resource usage and performance monitoring data.

## Usage Examples

### Human Data (hg38)
```json
{
  "GVSRealignVetIndels.dataset_name": "human_gvs_dataset",
  "GVSRealignVetIndels.project_id": "my-human-project", 
  "GVSRealignVetIndels.vet_table_number": 7
}
```

### Canine Data
```json
{
  "GVSRealignVetIndels.dataset_name": "canine_gvs_dataset",
  "GVSRealignVetIndels.project_id": "my-canine-project",
  "GVSRealignVetIndels.vet_table_number": 1,
  "GVSRealignVetIndels.reference_fasta": "gs://my-bucket/UU_Cfam_GSD_1.0_ROSY.fasta"
}
```

### High-Performance Configuration
```json
{
  "GVSRealignVetIndels.dataset_name": "large_dataset",
  "GVSRealignVetIndels.project_id": "production-project",
  "GVSRealignVetIndels.vet_table_number": 15,
  "GVSRealignVetIndels.RealignSingleSample.memory_gb": 16,
  "GVSRealignVetIndels.RealignSingleSample.disk_gb": 100,
  "GVSRealignVetIndels.RealignSingleSample.cpu_count": 4,
  "GVSRealignVetIndels.RealignSingleSample.preemptible_tries": 1
}
```

## Resource Requirements

### Typical Resource Usage
- **Memory**: 8-16 GB per sample task
- **Disk**: 50-100 GB per sample task  
- **CPU**: 2-4 cores per sample task
- **Runtime**: 10-60 minutes per sample (varies by variant density)

### Scaling Considerations
- Processing 4000 samples in parallel requires significant compute quota
- Consider adjusting `preemptible_tries` based on cost vs. reliability needs
- Increase resources for samples with high variant density

## Prerequisites

1. **BigQuery Access**: Read permissions on VET tables
2. **Reference Genome**: Accessible FASTA file with index (.fai) and dictionary (.dict)
3. **GATK Docker**: Access to GATK Docker images
4. **Cloud Storage**: Write permissions for output files
5. **Compute Quota**: Sufficient quota for parallel sample processing

## Integration with GVS Pipeline

This workflow is part of the broader Genome Variant Store (GVS) pipeline:

1. **Upstream**: VET tables populated by variant ingestion tools
2. **This Workflow**: Indel realignment and quality analysis  
3. **Downstream**: Cleaned variants feed into cohort analysis and extraction workflows

## Error Handling

The workflow includes:
- Input validation (documented in comments)
- Retry logic via `preemptible_tries`
- Resource monitoring via background scripts
- Comprehensive error logging

## Performance Optimization

For optimal performance:
1. **Use preemptible instances** for cost savings
2. **Adjust memory** based on variant density per sample
3. **Monitor disk usage** for samples with many variants
4. **Consider regional data locality** for reference files

## Version Tracking

The workflow automatically tracks:
- GATK Docker image versions
- Git hash of the workflow code
- Tool versions used in processing

This ensures reproducibility and enables debugging of results.
