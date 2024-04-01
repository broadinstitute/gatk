# Running the Genome Variant Store (GVS) Pipeline for AoU

## Setup
- Create a Terra workspace
  - Request a new AoU workspace as described in Section 3.D. of [AoU DRC Protocols](https://docs.google.com/document/d/1ooK0wbHLgSueiepjTTyLLI6Zz7vi1-GhKTjmCd8ZHwU/edit?usp=sharing).
  - As described in the "Getting Started" of [Operational concerns for running Hail in Terra Cromwell/WDL](https://docs.google.com/document/d/1_OY2rKwZ-qKCDldSZrte4jRIZf4eAw2d7Jd-Asi50KE/edit?usp=sharing), this workspace will need permission in Terra to run Hail dataproc clusters within WDL. Contact Emily to request this access as part of setting up the new workspace.
  - There is a quota that needs to be upgraded for the process of Bulk Ingest.
    When we ingest data, we use the Write API, which is part of BQ’s Storage API. Since we are hitting this API with so much data all at once, we want to increase our CreateWriteStream quota. Follow the [Quota Request Template](workspace/CreateWriteStreamRequestIncreasedQuota.md).
    Once that quota has been increased, the `load_data_batch` value needs to be updated based on calculations in the [Quota Request Template](workspace/CreateWriteStreamRequestIncreasedQuota.md) doc. Even if no increased quota is granted, this doc goes over how to choose the value for this param.
  - Create and push a feature branch (e.g. `EchoCallset`) based off the `ah_var_store` branch to the GATK GitHub repo.
    - Update the .dockstore.yml file on that feature branch to add the feature branch for all the WDLs that will be loaded into the workspace in the next step.
- Once the requested workspace has been created and permissioned, populate with the following WDLs:
  - [GvsBulkIngestGenomes](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsBulkIngestGenomes) workflow
  - [GvsWithdrawSamples](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsWithdrawSamples) workflow
  - [GvsPopulateAltAllele](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsPopulateAltAllele) workflow
  - [GvsCreateFilterSet](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsCreateFilterSet) workflow
  - [GvsExtractAvroFilesForHail](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsExtractAvroFilesForHail) workflow
  - [GvsPrepareRangesCallset](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsPrepareRangesCallset) workflow
  - [GvsExtractCallset](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsExtractCallset) workflow
  - [GvsExtractCallsetPgenMerged](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsExtractCallsetPgenMerged) workflow
  - [GvsCallsetStatistics](https://dockstore.org/workflows/github.com/broadinstitute/gatk/GvsCallsetStatistics) workflow
  - [GvsCalculatePrecisionAndSensitivity](https://dockstore.org/workflows/github.com/broadinstitute/gatk/GvsCalculatePrecisionAndSensitivity) workflow
  - [GvsCallsetCost](https://dockstore.org/workflows/github.com/broadinstitute/gatk/GvsCallsetCost) workflow
  - [GvsCreateVDS](https://dockstore.org/workflows/github.com/broadinstitute/gatk/GvsCreateVDS) workflow
  - [GvsCreateVATfromVDS](https://dockstore.org/workflows/github.com/broadinstitute/gatk/GvsCreateVATfromVDS) workflow
  - [GvsValidateVat](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsValidateVat) workflow
- Install the [Fetch WGS metadata for samples from list](./workspace/Fetch%20WGS%20metadata%20for%20samples%20from%20list.ipynb) python notebook in the workspace that has been created.
  - Placed the file with the list of the new samples to ingest in a GCS location the notebook (running with your @pmi-ops account) will have access to.  This will grab the samples from the workspace where they were reblocked and bring them into this callset workspace.
  - Run the steps in the notebook:
      - In section 1.2, Set the `sample_list_file_path` variable in that notebook to the path of the file
      - Run the cells up to and through, section 1.4 ("Perform the copy, in batches")
      - If you want to automatically break up the new samples into smaller sample sets, then run the "now that the data have been copied, you can make sample sets if you wish" step. Set the `SUBSET_SIZE` and `set_name` variables to customize.
- **NOTE** If you want to create a large sample set after you have run the notebook, Terra provides (and recommends you use) this python [script](https://github.com/broadinstitute/firecloud-tools/tree/master/scripts/import_large_tsv) which allows you to upload a sample set to the workspace.
- Create a dataset in the Google project. Make sure that when you are creating the dataset that you set the `location type` to be `Multi-Region`.
- Make a note of the Google project ID ("aou-genomics-curation-prod"), dataset name (e.g. "aou_wgs" — if it does not exist be sure to create one before running any workflows) and callset identifier (e.g. "Bravo") as these will be inputs (`project_id`, `dataset_name` and `call_set_identifier`) to all or most of the GVS workflows. The [naming conventions for other aspects of GVS datasets are outlined here](https://docs.google.com/document/d/1pNtuv7uDoiOFPbwe4zx5sAGH7MyxwKqXkyrpNmBxeow).
- Once the **non-control** samples have been fully ingested into BQ using the `GvsBulkIngestGenomes` workflow, the **control** samples can be manually added to the workspace and loaded in separately

## The Pipeline
1. `GvsBulkIngestGenomes` workflow
   - For use with **non-control** samples only! To ingest control samples (required for running `GvsCalculatePrecisionAndSensitivity`), use the`GvsAssignIds` and `GvsImportGenomes` workflows described below.
   - Set `sample_id_column_name` to "research_id" to use the shorter unique ID from AoU for the `sample_name` column.
   - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
   - Set `process_vcf_headers` to true to include sample header information in your database.
     - This slows down the ingestion process by about 30%, so if you dont care about this addtional data, leave it as false. 
     - Once these headers have been loaded, they can be queried to provide helpful context on the samples pipeline, such as DRAGEN version.
   - **NOTE** Be sure to set the input `drop_state` to ZERO (this will have the effect of dropping GQ0 reference blocks) and set `use_compressed_references` to true (this will further compress the reference data).
   - `GvsBulkIngestGenomes` performs the functions of both `GvsAssignIds` and `GvsImportGenomes` with a much more scalable design. Detailed bulk ingest documentation can be found [here](../gvs-bulk-ingest-details.md).
   - Note: In case of mistakenly ingesting a large number of bad samples, instructions for removing them can be found in [this Jira ticket](https://broadworkbench.atlassian.net/browse/VS-1206)
1. `GvsWithdrawSamples` workflow
   - Run if there are any samples to withdraw from the last callset.
   - When you run the `GvsWithdrawSamples` workflow, you should inspect the output of the workflow.
     - The output `num_samples_withdrawn` indicates the number of samples that have been withdrawn. This number should agree with that which you expect.
     - If the workflow fails, it may have failed if the list of smaples that was supplied to it includes samples that have not yet been ingested. To determine if this is the case, inspect the output (STDOUT) of the workflow and if it includes a list of samples that need to be ingested, then do so (or investigate the discreprancy). Note that there is a boolean variable `allow_uningested_samples` for this workflow that will allow it to pass if this condition occurs.
1. **TBD Workflow to soft delete samples**
1. `GvsPopulateAltAllele` workflow
   - This step loads data into the `alt_allele` table from the `vet_*` tables in preparation for running the filtering step.
   - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
1. `GvsCreateFilterSet` workflow
   - This step calculates features from the `alt_allele` table, and trains the VETS filtering model along with site-level QC filters and loads them into BigQuery into a series of `filter_set_*` tables.
   - See [naming conventions doc](https://docs.google.com/document/d/1pNtuv7uDoiOFPbwe4zx5sAGH7MyxwKqXkyrpNmBxeow) for guidance on what to use for `filter_set_name`, which you will need to keep track of for the `GvsExtractAvroFilesForHail` WDL. If, for some reason, this step needs to be run multiple times, be sure to use a different `filter_set_name` (the doc has guidance for this, as well).
   - **NOTE** For the Echo callset, the `CreateFilteredScoredSNPsVCF` failed twice, needing to be rerun with first more disk and then more memory. Since that time the call to the underlying method has been changed, but it is possible that the future runs may also need to have additional memory or disk - in which case the defaults used by that task may need to be updated.
   - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
1. `GvsExtractAvroFilesForHail` workflow
   - This workflow extracts the data in BigQuery and transforms it into Avro files in a Google bucket, incorporating the VETS filter set data.
   - The extracted Avro files will then be used as an input for `GvsCreateVDS` workflow described below.
   - This workflow needs to be run with the `filter_set_name` from `GvsCreateFilterSet` step.
   - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
1. `GvsCreateVDS` workflow
   - This step creates a VDS based on the Avro files generated from the `GvsExtractAvroFilesForHail` workflow above.
   - You can find what the `avro_path` input should be by going to the `GvsExtractAvroFilesForHail` run in Job Manager; the `GenerateHailScripts` task has an input of `avro_prefix` which is the location of the files created by that workflow.
   - The `vds_path` path input to this workflow represents the output path for the VDS. VDSes should be written under the AoU delivery bucket `gs://prod-drc-broad/`. Ask Lee for the exact path to use for the VDS in `#dsp-variants`.
   - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
   - Once a VDS has been created the Variants team will also generate callset statistics using `GvsCallsetStatistics` as described below. The Variants team then forwards both the path to the VDS and the output callset statistics TSV to Lee to quality check the VDS.
   - If you are debugging a Hail-related issue, you may want to set `leave_hail_cluster_running_at_end` to `true` and refer to [the suggestions for debugging issues with Hail](HAIL_DEBUGGING.md).
1. `GvsCallsetStatistics` workflow
    - You will need to run `GvsPrepareRangesCallset` workflow for callset statistics first, if it has not been run already
       - This workflow transforms the data in the vet tables into a schema optimized for callset stats creation and for calculating sensitivity and precision.
       - The `only_output_vet_tables` input should be set to `true` (the default value is `false`).
       - The `enable_extract_table_ttl` input should be set to `true` (the default value is `false`), which will add a TTL of two weeks to the tables it creates.
       - See [naming conventions doc](https://docs.google.com/document/d/1pNtuv7uDoiOFPbwe4zx5sAGH7MyxwKqXkyrpNmBxeow) for guidance on what to use for `extract_table_prefix` or cohort prefix, which you will need to keep track of for the callset stats.
       - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
    - You will need to have the "BigQuery Data Viewer" role for your @pmi-ops proxy group on the `spec-ops-aou:gvs_public_reference_data.gnomad_v3_sites` table
    - This workflow needs to be run with the `extract_table_prefix` input from `GvsPrepareRangesCallset` step.
    - This workflow needs to be run with the `filter_set_name` input from `GvsCreateFilterSet` step.
    - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
1. `GvsExtractCallset` / `GvsExtractCallsetPgenMerged` workflow
    - You will need to run the `GvsPrepareRangesCallset` workflow for each "[Region](https://support.researchallofus.org/hc/en-us/articles/14929793660948-Smaller-Callsets-for-Analyzing-Short-Read-WGS-SNP-Indel-Data-with-Hail-MT-VCF-and-PLINK)" (interval list) for which a PGEN or VCF deliverable is required for the callset.
        - This workflow transforms the data in the vet, ref_ranges, and samples tables into a schema optimized for extract.
        - The `enable_extract_table_ttl` input should be set to `true` (the default value is `false`), which will add a TTL of two weeks to the tables it creates.
        - `extract_table_prefix` should be set to a name that is unique to the given Region / interval list. See the [naming conventions doc](https://docs.google.com/document/d/1pNtuv7uDoiOFPbwe4zx5sAGH7MyxwKqXkyrpNmBxeow) for guidance on what to use.
        - Specify the `interval_list` appropriate for the PGEN / VCF extraction run you are performing.
        - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
    - Specify the same `call_set_identifier`, `dataset_name`, `project_id`, `extract_table_prefix`, and `interval_list` that were used in the `GvsPrepareRangesCallset` run documented above.
    - Specify the `interval_weights_bed` appropriate for the PGEN / VCF extraction run you are performing. `gs://gvs_quickstart_storage/weights/gvs_full_vet_weights_1kb_padded_orig.bed` is the interval weights BED used for Quickstart.
    - For `GvsExtractCallsetPgenMerged` only, select the workflow option "Retry with more memory" and choose a "Memory retry factor" of 1.5
    - These workflows do not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
1. `GvsCalculatePrecisionAndSensitivity` workflow
    - Please see the detailed instructions for running the Precision and Sensitivity workflow [here](../../tieout/AoU_PRECISION_SENSITIVITY.md).
1. `GvsCallsetCost` workflow
    - This workflow calculates the total BigQuery cost of generating this callset (which is not represented in the Terra UI total workflow cost) using the above GVS workflows; it's used to calculate the cost as a whole and by sample.

## Deliverables (via email to stakeholders once the above steps are complete)

The Callset Stats and S&P files can be simply `gsutil cp`ed to the AoU delivery bucket since they are so much smaller.
1. GCS location of the VDS in the AoU delivery bucket
2. Fully qualified name of the BigQuery dataset (composed of the `project_id` and `dataset_name` inputs from the workflows)
3. GCS location of the CSV output from `GvsCallsetStatistics` workflow in the AoU delivery bucket
4. GCS location of the TSV output from `GvsCalculatePrecisionAndSensitivity` in the AoU delivery bucket



## Running the VAT pipeline
To create a BigQuery table of variant annotations, you may follow the instructions here:
[process to create variant annotations table](../../variant_annotations_table/README.md)

The pipeline takes in the VDS and outputs a variant annotations table in BigQuery.



