# Running Exomes on GVS

This document describes the changes necessary to run exome gVCFs through the GVS workflow. The changes needed to run exomes primarily involve using different parameters.
**NOTE** Currently this document is written to be at the developer level (that is for experienced developers). For other docs (specifically, for our beta users) see https://github.com/broadinstitute/gatk/tree/ah_var_store/scripts/variantstore/beta_docs/

**NOTE** For Exome we want to use the latest BGE exome interval list:
gs://gcp-public-data--broad-references/hg38/v0/bge_exome_calling_regions.v1.1.interval_list

## Setup
- Create a Terra workspace
- Populate the workspace with the following workflows:
    - [GvsBulkIngestGenomes](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsBulkIngestGenomes) workflow
    - [GvsAssignIds](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsAssignIds) workflow
    - [GvsImportGenomes](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsImportGenomes) workflow
    - [GvsPopulateAltAllele](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsPopulateAltAllele) workflow
    - [GvsCreateFilterSet](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsCreateFilterSet) workflow
    - [GvsPrepareRangesCallset](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsPrepareRangesCallset) workflow
    - [GvsExtractCallset](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsExtractCallset) workflow
    - [GvsCalculatePrecisionAndSensitivity](https://dockstore.org/workflows/github.com/broadinstitute/gatk/GvsCalculatePrecisionAndSensitivity) workflow

## The Pipeline
1. `GvsBulkIngestGenomes` workflow
    - Run this workflow in order to load all samples into the database tables so that they can be run through the GVS workflow. This workflow encompasses the tasks described below (in `GvsAssignIds` and `GvsImportGenomes`)
   - Run at the `sample set` level ("Step 1" in workflow submission) with a sample set of all the new samples to be included in the callset.
    - **NOTE** For Exomes, use `gs://gcp-public-data--broad-references/hg38/v0/bge_exome_calling_regions.v1.1.interval_list` for the `interval_list`
- OR you can run the two workflows that `GvsBulkIngestGenomes` calls (for instance if you also need to load control samples) 
1. `GvsAssignIds` workflow
    - To optimize the GVS internal queries, each sample must have a unique and consecutive integer ID assigned. Running the `GvsAssignIds` will create a unique GVS ID for each sample (`sample_id`) and update the BQ `sample_info` table (creating it if it doesn't exist). This workflow takes care of creating the BQ `vet_*`, `ref_ranges_*` and `cost_observability` tables needed for the sample IDs generated.
    - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
    -  The `external_sample_names` input should be the GCS path of a text file that lists all the sample names (external sample IDs).
    - If new controls are being added, they need to be done in a separate run, with the `samples_are_controls` input set to "true" (the referenced Data columns may also be different, e.g. "this.control_samples.control_sample_id" instead of "this.samples.research_id").
2. `GvsImportGenomes` workflow
    - This will import the re-blocked gVCF files into GVS. The workflow will check whether data for that sample has already been loaded into GVS. It is designed to be re-run (with the same inputs) if there is a failure during one of the workflow tasks (e.g. BigQuery write API interrupts).
    - Run at the `sample set` level ("Step 1" in workflow submission).  You can either run this on a sample_set of all the samples and rely on the workflow logic to break it up into batches.
    - You will want to set the `external_sample_names`, `input_vcfs` and `input_vcf_indexes` inputs based on the columns in the workspace Data table, e.g. "this.samples.research_id", "this.samples.reblocked_gvcf_v2" and "this.samples.reblocked_gvcf_index_v2".
   - **NOTE** For Exomes, use `gs://gcp-public-data--broad-references/hg38/v0/bge_exome_calling_regions.v1.1.interval_list` for the `interval_list`

3. `GvsPopulateAltAllele` workflow
    - This step loads data into the `alt_allele` table from the `vet_*` tables in preparation for running the filtering step.
    - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
4. `GvsCreateFilterSet` workflow
    - This step calculates features from the `alt_allele` table, and trains the VQSR filtering model along with site-level QC filters and loads them into BigQuery into a series of `filter_set_*` tables.
    - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
    - **NOTE** For Exomes, use `gs://gcp-public-data--broad-references/hg38/v0/bge_exome_calling_regions.v1.1.interval_list` for the `interval_list`
5. `GvsPrepareRangesCallset` workflow
    - This workflow transforms the data in the vet tables into a schema optimized for callset stats creation and for calculating sensitivity and precision.
    - This workflow may only need to be run once (for controls extract for Precision and Sensitivity). Run it with `control_samples` set to "true" (the default value is `false`).
    - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
    - **NOTE** For Exomes, set the parameter `use_interval_weights` to `false`.  This avoids a bug seen in WeightedSplitIntervals when using exomes, forcing it to use the standard version of SplitIntervals.
6. `GvsCalculatePrecisionAndSensitivity` workflow
    - This workflow needs to be run with the control sample chr20 vcfs from `GvsExtractCallset` step which were placed in the `output_gcs_dir`.
    - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
    - **NOTE** For Exomes, use `gs://gvs-internal/truth/HG001.exome_evaluation_regions.v1.1.bed` as the "truth" bed for NA12878/HG001.

