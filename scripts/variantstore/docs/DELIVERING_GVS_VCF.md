# Genome Variant Store (GVS) Pipeline for VCF Output

## Setup
- Create a Terra workspace and share it with `dsp-variant-team@firecloud.org` so the rest of the team can see it (Reader) and with the individual users who will be running the workflows (Owner).
- Populate the workspace with the following:
  - [Fetch WGS metadata for samples from list](http://app.terra.bio/#workspaces/allofus-drc-wgs-dev/GVS%20AoU%20WGS%20Charlie/notebooks/launch/Fetch%20WGS%20metadata%20for%20samples%20from%20list.ipynb) notebook
  - [GvsAssignIds](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsAssignIds) workflow
  - [GvsImportGenomes](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsImportGenomes) workflow
  - [GvsWithdrawSamples](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsWithdrawSamples) workflow
  - [GvsPopulateAltAllele](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsPopulateAltAllele) workflow
  - [GvsCreateFilterSet](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsCreateFilterSet) workflow
  - [GvsPrepareRangesCallset](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsPrepareRangesCallset) workflow
  - [GvsExtractCallset](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsExtractCallset) workflow
  - [GvsCalculatePrecisionAndSensitivity](https://dockstore.org/workflows/github.com/broadinstitute/gatk/GvsCalculatePrecisionAndSensitivity) workflow (if applicable)
  - [GvsCallsetCost](https://dockstore.org/workflows/github.com/broadinstitute/gatk/GvsCallsetCost) workflow
  - [GvsCallsetStatistics](https://dockstore.org/workflows/github.com/broadinstitute/gatk/GvsCallsetStatistics) workflow (if applicable)
- Run the "Fetch WGS metadata for samples from list" notebook after you have placed the file with the list of the new samples to ingest in a GCS location the notebook (running with your @pmi-ops account) will have access to.  This will grab the samples from the workspace where they were reblocked and bring them into this callset workspace.
  - Set the `sample_list_file_path` variable in that notebook to the path of the file
  - Run the "now that the data have been copied, you can make sample sets if you wish" step if you want to automatically break up the new samples into smaller sample sets.  Set the `SUBSET_SIZE` and `set_name` variables to customize.
- **NOTE** If you want to create a large sample set after you have run the notebook, Terra provides (and recommends you use) this Python [script](https://github.com/broadinstitute/firecloud-tools/tree/master/scripts/import_large_tsv) which allows you to upload a sample set to the workspace.
- You will want to increase the Google quotas for the workspace project (you can find this in the workspace dashboard under Cloud Information > Google Project ID) to these levels (all in the workspace region):
  - Persistent Disk Standard (GB): 1,000,000 (1 PB)
  - CPUs: 6,400
  - In-use IP addresses: 5,000 (this is the most challenging one and will probably require contacting the GCP account team to facilitate)
  - VM instances: 64,000
  - Concurrent connections per project for small regions per region: 2,000
  - AppendRows throughput per project for small regions per minute per region: 2,000
  - Concurrent connections per project for small regions per region: 2,000
- Make a note of the Google project ID, BigQuery dataset name, and callset identifier (e.g. "Bravo") as these will be inputs (`project_id`, `dataset_name` and `call_set_identifier`) to all or most of the GVS workflows. The [naming conventions for other aspects of GVS datasets are outlined here](https://docs.google.com/document/d/1pNtuv7uDoiOFPbwe4zx5sAGH7MyxwKqXkyrpNmBxeow).

## The Pipeline
1. `GvsAssignIds` workflow
   - To optimize the GVS internal queries, each sample must have a unique and consecutive integer ID assigned. Running the `GvsAssignIds` workflow will create a unique GVS ID for each sample (`sample_id`) and update the BQ `sample_info` table (creating it if it doesn't exist). This workflow takes care of creating the BQ `vet_*`, `ref_ranges_*` and `cost_observability` tables needed for the sample IDs generated.
   - Run at the `sample set` level ("Step 1" in workflow submission) with a sample set of all the new samples to be included in the callset (created by the "Fetch WGS metadata for samples from list" notebook mentioned above).
   - You will want to set the `external_sample_names` input based on the column in the workspace Data table, e.g. "this.samples.research_id".
   - If new controls are being added, they need to be done in a separate run, with the `samples_are_controls` input set to "true" (the referenced Data columns may also be different, e.g. "this.control_samples.control_sample_id" instead of "this.samples.research_id").
2. `GvsImportGenomes` workflow
   - This will import the re-blocked gVCF files into GVS. The workflow will check whether data for that sample has already been loaded into GVS. It is designed to be re-run (with the same inputs) if there is a failure during one of the workflow tasks (e.g. BigQuery write API interrupts).
   - Run at the `sample set` level ("Step 1" in workflow submission).  You can either run this on a sample_set of all the samples and rely on the workflow logic to break it up into batches (or manually set the `load_data_batch_size` input) or run it on smaller sample_sets created by the "Fetch WGS metadata for samples from list" notebook mentioned above.  
   - You will want to set the `external_sample_names`, `input_vcfs` and `input_vcf_indexes` inputs based on the columns in the workspace Data table, e.g. "this.samples.research_id", "this.samples.reblocked_gvcf_v2" and "this.samples.reblocked_gvcf_index_v2".
   - **NOTE** It appears that there is a rawls limit on the size of the input (list of gvcf files and indexes) per workflow run. 25K samples in a list worked for the Intermediate call set, 50K did not.
3. `GvsWithdrawSamples` workflow
   - Run if there are any samples to withdraw from the last callset.
4. `GvsPopulateAltAllele` workflow
   - This step loads data into the `alt_allele` table from the `vet_*` tables in preparation for running the filtering step.
   - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
5. `GvsCreateFilterSet` workflow
   - This step calculates features from the `alt_allele` table, and trains the VQSR filtering model along with site-level QC filters and loads them into BigQuery into a series of `filter_set_*` tables.
   - See [naming conventions doc](https://docs.google.com/document/d/1pNtuv7uDoiOFPbwe4zx5sAGH7MyxwKqXkyrpNmBxeow) for guidance on what to use for `filter_set_name`, which you will need to keep track of for the `GvsExtractCallset` WDL. If, for some reason, this step needs to be run multiple times, be sure to use a different `filter_set_name` (the doc has guidance for this, as well).
   - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
6. `GvsPrepareRangesCallset` workflow
   - This workflow transforms the data in the vet and ref_ranges tables into a schema optimized for VCF generation during the Extract step.
   - It will need to be run twice if you are generating Precision and Sensitivity, once with `control_samples` set to "true" (see [naming conventions doc](https://docs.google.com/document/d/1pNtuv7uDoiOFPbwe4zx5sAGH7MyxwKqXkyrpNmBxeow) for guidance on what to use for `extract_table_prefix` or cohort prefix, which you will need to keep track of for the `GvsExtractCallset` WDL); the default value is `false`.
   - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
7. `GvsExtractCallset` workflow
   - It is strongly suggested that you provide the `output_gcs_dir` input which will collect the VCFs, VCF indexes, interval lists, sample name list and manifest in one place.  Make sure it is a path your Terra proxy account has access to (the easiest way to do this is to set it to be within the workspace bucket where you are running these workflows).
   - This workflow extracts the data in BigQuery and transforms it into a sharded joint called VCF incorporating the VQSR filter set data.  We will probably not run this on callsets of more than 100K samples.
   - It also needs to be run twice if you are generating Precision and Sensitivity, once with `control_samples` set to "true", and with the `filter_set_name` and `extract_table_prefix` from step 5 & 6.  Include a valid (and secure) "output_gcs_dir" parameter, which is where the VCF, interval list, manifest, and sample name list files will go.
   - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
8. `GvsCallsetStatistics` workflow
     - Add the "BigQuery Data Viewer" role for your Terra proxy group on the `spec-ops-aou:gvs_public_reference_data.gnomad_v3_sites` table 
     - You will need the `filter_set_name` and `extract_table_prefix` from step 5 & 6.
9. `GvsCalculatePrecisionAndSensitivity` workflow
    - Add the "Storage Object Viewer" access granted for your Terra proxy group on the `gs://broad-dsp-spec-ops/gvs/truth` directory
10. `GvsCallsetCost`
    - The cost from this callset, which represents the total BigQuery cost (which is not represented in the Terra UI total workflow cost) from the GVS pipeline workflows, is used to calculate the cost of the callset as a whole and by sample.
