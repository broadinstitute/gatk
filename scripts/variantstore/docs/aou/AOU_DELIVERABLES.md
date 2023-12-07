# Running the Genome Variant Store (GVS) Pipeline for AoU

## Setup
- Create a Terra workspace
  - Request a new AoU workspace as described in Section 3.D. of [AoU DRC Protocols](https://docs.google.com/document/d/1ooK0wbHLgSueiepjTTyLLI6Zz7vi1-GhKTjmCd8ZHwU/edit?usp=sharing).
  - As described in the "Getting Started" of [Operational concerns for running Hail in Terra Cromwell/WDL](https://docs.google.com/document/d/1_OY2rKwZ-qKCDldSZrte4jRIZf4eAw2d7Jd-Asi50KE/edit?usp=sharing), this workspace will need permission in Terra to run Hail dataproc clusters within WDL. Contact Emily to request this access as part of setting up the new workspace.
- Create and push a feature branch (e.g. `EchoCallset`) based off the `ah_var_store` branch to the GATK GitHub repo.
    - Update the .dockstore.yml file on that feature branch to add the feature branch for all the WDLs that will be loaded into the workspace in the next step.
- Once the workspace requested above has been created and permissioned, populate with the following WDLs:
  - [GvsBulkIngestGenomes](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsBulkIngestGenomes) workflow
  - [GvsAssignIds](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsAssignIds) workflow
  - [GvsImportGenomes](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsImportGenomes) workflow
  - [GvsWithdrawSamples](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsWithdrawSamples) workflow
  - [GvsPopulateAltAllele](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsPopulateAltAllele) workflow
  - [GvsCreateFilterSet](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsCreateFilterSet) workflow
  - [GvsExtractAvroFilesForHail](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsExtractAvroFilesForHail) workflow
  - [GvsPrepareRangesCallset](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsPrepareRangesCallset) workflow
  - [GvsExtractCallset](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsExtractCallset) workflow
  - [GvsCallsetStatistics](https://dockstore.org/workflows/github.com/broadinstitute/gatk/GvsCallsetStatistics) workflow
  - [GvsCalculatePrecisionAndSensitivity](https://dockstore.org/workflows/github.com/broadinstitute/gatk/GvsCalculatePrecisionAndSensitivity) workflow
  - [GvsCallsetCost](https://dockstore.org/workflows/github.com/broadinstitute/gatk/GvsCallsetCost) workflow
  - [GvsCreateVDS](https://dockstore.org/workflows/github.com/broadinstitute/gatk/GvsCreateVDS) workflow
- Install the [Fetch WGS metadata for samples from list](./workspace/Fetch%20WGS%20metadata%20for%20samples%20from%20list.ipynb) python notebook in the workspace that has been created.
  - Placed the file with the list of the new samples to ingest in a GCS location the notebook (running with your @pmi-ops account) will have access to.  This will grab the samples from the workspace where they were reblocked and bring them into this callset workspace.
  - Run the steps in the notebook:
      - In section 1.2, Set the `sample_list_file_path` variable in that notebook to the path of the file
      - Run the cells up to and through, section 1.4 ("Perform the copy, in batches")
      - If you want to automatically break up the new samples into smaller sample sets, then run the "now that the data have been copied, you can make sample sets if you wish" step. Set the `SUBSET_SIZE` and `set_name` variables to customize.
- **NOTE** If you want to create a large sample set after you have run the notebook, Terra provides (and recommends you use) this python [script](https://github.com/broadinstitute/firecloud-tools/tree/master/scripts/import_large_tsv) which allows you to upload a sample set to the workspace.
- Create a dataset in the Google project. Make sure that when you are creating the dataset that you set the `location type` to be `Multi-Region`
- Make a note of the Google project ID ("aou-genomics-curation-prod"), dataset name (e.g. "aou_wgs" — if it does not exist be sure to create one before running any workflows) and callset identifier (e.g. "Bravo") as these will be inputs (`project_id`, `dataset_name` and `call_set_identifier`) to all or most of the GVS workflows. The [naming conventions for other aspects of GVS datasets are outlined here](https://docs.google.com/document/d/1pNtuv7uDoiOFPbwe4zx5sAGH7MyxwKqXkyrpNmBxeow).

## The Pipeline
1. `GvsBulkIngestGenomes` workflow
   - For use with **non-control** samples only! To ingest control samples (required for running `GvsCalculatePrecisionAndSensitivity`), use the`GvsAssignIds` and `GvsImportGenomes` workflows described below.
   - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
   - **NOTE** Be sure to set the input `drop_state` to ZERO (this will have the effect of dropping GQ0 reference blocks) and `use_compressed_references` to true (this will further compress the reference data).
   - `GvsBulkIngestGenomes` performs the functions of both `GvsAssignIds` and `GvsImportGenomes` with a much more scalable design. Detailed bulk ingest documentation can be found [here](../gvs-bulk-ingest-details.md).
1. `GvsAssignIds` workflow
   - For use with **control** samples only!
   - To optimize the GVS internal queries, each sample must have a unique and consecutive integer ID assigned. Running the `GvsAssignIds` will create a unique GVS ID for each sample (`sample_id`) and update the BQ `sample_info` table (creating it if it doesn't exist). This workflow takes care of creating the BQ `vet_*`, `ref_ranges_*` and `cost_observability` tables needed for the sample IDs generated.
   - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
   - You will want to set the `external_sample_names` input based on the GCS path of a text file that lists all the sample names (external sample IDs).
   - If new controls are being added, they need to be done in a separate run, with the `samples_are_controls` input set to "true" (the referenced Data columns may also be different, e.g. "this.control_samples.control_sample_id" instead of "this.samples.research_id").
   - **NOTE** Set `use_compressed_references` to true.
1. `GvsImportGenomes` workflow
   - For use with **control** samples only!
   - This will import the re-blocked gVCF files into GVS. The workflow will check whether data for that sample has already been loaded into GVS. It is designed to be re-run (with the same inputs) if there is a failure during one of the workflow tasks (e.g. BigQuery write API interrupts).
   - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
   - You will need to set the `external_sample_names` input to be the GCS path of a text file that lists all the sample names (external sample IDs).
   - You will need to set the `input_vcfs` input to be the GCS path of a text file that lists all the re-blocked gVCFs (in the same order as the file used in `external_sample_names`). You can generate the text file by exporting the full paths to the VCFs from the samples table in Terra and formatting appropriately.
   - You will need to set the `input_vcf_indexes` input to be the GCS path of a text file that lists all the index files of the re-blocked gVCFs (in the same order as the file used in `external_sample_names`). You can generate this text file by exporting the full paths to the VCF index files from the samples table in Terra and formatting appropriately.
   - **NOTE** Be sure to set the input `drop_state` to ZERO. This will have the effect of dropping GQ0 reference blocks.
   - **NOTE** Set `use_compressed_references` to true.
   - **NOTE** It appears that there is a rawls limit on the size of the input (list of gvcf files and indexes) per workflow run. 25K samples in a list worked for the Intermediate call set, 50K did not.
1. `GvsWithdrawSamples` workflow
   - Run if there are any samples to withdraw from the last callset.
1. **TBD Workflow to soft delete samples**
1. `GvsPopulateAltAllele` workflow
   - This step loads data into the `alt_allele` table from the `vet_*` tables in preparation for running the filtering step.
   - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
   - This step is cumulative (as are all of the steps prior), and it is the last step that is--so be sure that all samples have been loaded or withdrawn before progressing to the next step
1. `GvsCreateFilterSet` workflow
   - This step calculates features from the `alt_allele` table, and trains the VETS filtering model along with site-level QC filters and loads them into BigQuery into a series of `filter_set_*` tables.
   - See [naming conventions doc](https://docs.google.com/document/d/1pNtuv7uDoiOFPbwe4zx5sAGH7MyxwKqXkyrpNmBxeow) for guidance on what to use for `filter_set_name`, which you will need to keep track of for the `GvsExtractAvroFilesForHail` WDL. If, for some reason, this step needs to be run multiple times, be sure to use a different `filter_set_name` (the doc has guidance for this, as well).
   - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
1. `GvsExtractAvroFilesForHail` workflow
   - This workflow extracts the data in BigQuery and transforms it into Avro files in a Google bucket, incorporating the VETS filter set data.
   - The extracted Avro files will then be used as an input for `GvsCreateVDS` workflow described below.
   - This workflow needs to be run with the `filter_set_name` from `GvsCreateFilterSet` step.
   - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
1. `GvsCreateVDS` workflow
   - This step creates a VDS based on the Avro files generated from the `GvsExtractAvroFilesForHail` workflow above.
   - For the `avro_path` input, Avro files from the above `GvsExtractAvroFilesForHail` step can be found in the Avro directory in the `OutputPath` task: `gs://fc-<workspace-id>/submissions/<submission id>/GvsExtractAvroFilesForHail/<workflow id>/call-OutputPath/avro`
   - The `vds_path` path input to this workflow represents the output path for the VDS. VDSes should be written under the AoU delivery bucket `gs://prod-drc-broad/`. Ask Lee for the exact path to use for the VDS in `#dsp-variants`.
   - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
   - Once a VDS has been created the Variants team will also generate callset statistics using `GvsCallsetStatistics` as described below. The Variants team then forwards both the path to the VDS and the output callset statistics TSV to Lee to quality check the VDS.
1. `GvsCallsetStatistics` workflow
    - You will need to run `GvsPrepareRangesCallset` workflow first, if it has not been run already
       - This workflow transforms the data in the vet tables into a schema optimized for callset stats creation and for calculating sensitivity and precision.
       - The `only_output_vet_tables` input should be set to `true` (the default value is `false`).
       - The `enable_extract_table_ttl` input should be set to `true` (the default value is `false`), which will add a TTL of two weeks to the tables it creates.
       - See [naming conventions doc](https://docs.google.com/document/d/1pNtuv7uDoiOFPbwe4zx5sAGH7MyxwKqXkyrpNmBxeow) for guidance on what to use for `extract_table_prefix` or cohort prefix, which you will need to keep track of for the callset stats.
       - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
    - You will need to have the "BigQuery Data Viewer" role for your @pmi-ops proxy group on the `spec-ops-aou:gvs_public_reference_data.gnomad_v3_sites` table
    - This workflow needs to be run with the `extract_table_prefix` input from `GvsPrepareRangesCallset` step.
    - This workflow needs to be run with the `filter_set_name` input from `GvsCreateFilterSet` step.
    - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
1. `GvsCalculatePrecisionAndSensitivity` workflow
    - Please see the detailed instructions for running the Precision and Sensitivity workflow [here](../../tieout/AoU_PRECISION_SENSITIVITY.md).
1. `GvsCallsetCost` workflow
    - This workflow calculates the total BigQuery cost of generating this callset (which is not represented in the Terra UI total workflow cost) using the above GVS workflows; it's used to calculate the cost as a whole and by sample.

## Deliverables (via email to stakeholders once the above steps are complete)
To get the VDS someplace secure and accessible for delivery, see [Delivering a VDS](vds/delivery/Delivering%20a%20VDS.md)
This includes details of how to get the VDS to the AoU delivery bucket as well as VDS naming conventions.
The Callset Stats and S&P files can be simply gsutil cp-ed to the AoU delivery bucket since they are so much smaller.
1. GCS location of the VDS in the AoU delivery bucket
2. Fully qualified name of the BigQuery dataset (composed of the `project_id` and `dataset_name` inputs from the workflows)
3. GCS location of the CSV output from `GvsCallsetStatistics` workflow in the AoU delivery bucket
4. GCS location of the TSV output from `GvsCalculatePrecisionAndSensitivity` in the AoU delivery bucket



## Running the VAT pipeline
To create a BigQuery table of variant annotations, you may follow the instructions here:
[process to create variant annotations table](../../variant_annotations_table/README.md)

The pipeline takes in the VDS and outputs a variant annotations table in BigQuery.



