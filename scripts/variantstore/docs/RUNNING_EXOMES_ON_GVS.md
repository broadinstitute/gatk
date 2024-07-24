# Running Exomes on GVS

This document describes the changes necessary to run exome gVCFs through the GVS workflow. The changes needed to run exomes primarily involve using different parameters. Currently this document is written to be at the developer level (that is for experienced developers). For other docs (specifically, for our beta users) see [https://github.com/broadinstitute/gatk/tree/ah_var_store/scripts/variantstore/beta_docs/]([https://github.com/broadinstitute/gatk/tree/ah_var_store/scripts/variantstore/beta_docs/)

## Setup

- Create a Terra workspace and a BigQuery dataset with the necessary corresponding permissions for your PROXY group.
- Populate the workspace with the following workflows:
  - [GvsBulkIngestGenomes](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsBulkIngestGenomes) workflow
  - [GvsAssignIds](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsAssignIds) workflow (only if you want to calculate Precision and Sensitivity)
  - [GvsImportGenomes](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsImportGenomes) workflow (only if you want to calculate Precision and Sensitivity)
  - [GvsPopulateAltAllele](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsPopulateAltAllele) workflow
  - [GvsCreateFilterSet](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsCreateFilterSet) workflow
  - [GvsPrepareRangesCallset](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsPrepareRangesCallset) workflow
  - [GvsExtractCallset](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsExtractCallset) workflow
  - [GvsCalculatePrecisionAndSensitivity](https://dockstore.org/workflows/github.com/broadinstitute/gatk/GvsCalculatePrecisionAndSensitivity) workflow (only if you want to calculate Precision and Sensitivity)

## The Pipeline
1. `GvsBulkIngestGenomes` workflow
   - Run this workflow in order to load all non-control samples into the database tables so that they can be run through the GVS workflow.
   - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
   - Set the `interval_list` input to `"gs://gcp-public-data--broad-references/hg38/v0/bge_exome_calling_regions.v1.1.interval_list"`
 
    To ingest control samples (which you will need to calculate Precision and Sensitivity), you will need to run the`GvsAssignIds` and `GvsImportGenomes` workflows just for them:
   1. `GvsAssignIds` workflow
      - This workflow is set up to be re-run (with the same inputs) if there is a failure, but be sure to check for an existing `sample_id_assignment_lock` table in your dataset; if it exists, delete it before re-running.
      - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
      - The `external_sample_names` input should be the GCS path of a text file that lists all the control sample names (external sample IDs).
      - Set the input `samples_are_controls` to `true`.
   1. `GvsImportGenomes` workflow
      - This will import the re-blocked gVCF files into GVS. The workflow will check whether data for that sample has already been loaded into GVS. It is designed to be re-run (with the same inputs) if there is a failure during one of the workflow tasks (e.g. BigQuery write API interrupts).
      - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
      - You will want to set the `external_sample_names`, `input_vcfs` and `input_vcf_indexes` inputs based on the GCP locations of files that contain lists of these values (in the same order). For NA12878/HG001, the files you will need for ingest are:
          - `input_vcfs`: `"gs://broad-gotc-test-storage/ExomeGermlineSingleSample/truth/scientific/master/D5327.NA12878/NA12878.rb.g.vcf.gz"`
          - `input_vcf_indexes`: `"gs://broad-gotc-test-storage/ExomeGermlineSingleSample/truth/scientific/master/D5327.NA12878/NA12878.rb.g.vcf.gz.tbi"`
      - Set the `interval_list` input to `"gs://gcp-public-data--broad-references/hg38/v0/bge_exome_calling_regions.v1.1.interval_list"`
1. `GvsPopulateAltAllele` workflow
   - This step loads data into the `alt_allele` table from the `vet_*` tables in preparation for running the filtering step.
   - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
1. `GvsCreateFilterSet` workflow
   - This step calculates features from the `alt_allele` table, and trains the VETS model along with site-level QC filters and loads them into BigQuery into a series of `filter_set_*` tables.
   - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
   - Set `GvsCreateFilterSet.interval_list` input to `"gs://gcp-public-data--broad-references/hg38/v0/bge_exome_calling_regions.v1.1.interval_list"` *Make sure you're not setting `interval_list` on any of the called tasks by mistake!*
1. `GvsPrepareRangesCallset` workflow
   - This workflow transforms the data in the vet tables into a schema optimized for VCF extraction.
   - This workflow will need to be run once to extract the callset as a whole and an additional time to create the files used to calculate Precision and Sensitivity (with `control_samples` set to `true`).
   - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
1. `GvsExtractCallset` workflow
   - This workflow is used to extract non-control samples only. See the `GvsCalculatePrecisionAndSensitivity` instructions below if running Precision and Sensitivity on control samples.
   - This workflow takes the tables created in the `Prepare` step to output joint VCF shards.
   - This workflow is run only to extract the callset into VCF shards (`control_samples` set to `false` explicitly or left at its default `false` value).
   - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
   - Set the `output_gcs_dir` to a GCS location to collect all the VCF shards into one place.
1. `GvsCalculatePrecisionAndSensitivity` workflow
   - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
   - Input values with truth inputs for NA12878/HG001:
     - `truth_beds`:  `[ "gs://gvs-internal/truth/HG001.exome_evaluation_regions.v1.1.bed" ]`
     - `truth_vcfs`: `[ "gs://gvs-internal/truth/HG001_exome_filtered.recode.vcf.gz.tbi" ]`
     - `truth_vcf_indices`: `"gs://gvs-internal/truth/HG001_exome_filtered.recode.vcf.gz.tbi"`
     - `vcf_eval_bed_file`: `"gs://gvs-internal/truth/HG001.exome_evaluation_regions.v1.1.bed"`
     - `chromosomes`: `[ "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22" ]`
     - `interval_list`: `"gs://gcp-public-data--broad-references/hg38/v0/bge_exome_calling_regions.v1.1.interval_list"`
     - `target_interval_list`: `"gs://gcp-public-data--broad-references/hg38/v0/bge_TwistAllianceClinicalResearchExome_Covered_Targets_hg38.interval_list"`
     - `is_wgs`: `false`
