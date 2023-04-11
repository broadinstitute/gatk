# Running the Genome Variant Store (GVS) Pipeline for AoU

## Setup
- Create a Terra workspace
  - using your @pmi-ops.org account 
  - with the `allofus-drc-wgs-dev` Terra Billing account 
  - in the `AOU_DRC_WGS_OPERATORS` Authorization Domain
  - Share it with `dsp-variant-team@firecloud.org` so the rest of the team can see it (Reader) and with the individual users who will be running the workflows (Owner).
- Populate the workspace with the following:
  - [Fetch WGS metadata for samples from list](http://app.terra.bio/#workspaces/allofus-drc-wgs-dev/GVS%20AoU%20WGS%20Charlie/notebooks/launch/Fetch%20WGS%20metadata%20for%20samples%20from%20list.ipynb) notebook
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
- Run the "Fetch WGS metadata for samples from list" notebook after you have placed the file with the list of the new samples to ingest in a GCS location the notebook (running with your @pmi-ops account) will have access to.  This will grab the samples from the workspace where they were reblocked and bring them into this callset workspace.
  - Set the `sample_list_file_path` variable in that notebook to the path of the file
  - Run the "now that the data have been copied, you can make sample sets if you wish" step if you want to automatically break up the new samples into smaller sample sets.  Set the `SUBSET_SIZE` and `set_name` variables to customize.
- **NOTE** If you want to create a large sample set after you have run the notebook, Terra provides (and recommends you use) this python [script](https://github.com/broadinstitute/firecloud-tools/tree/master/scripts/import_large_tsv) which allows you to upload a sample set to the workspace.
- Make a note of the Google project ID ("aou-genomics-curation-prod"), dataset name ("aou_wgs") and callset identifier (e.g. "Bravo") as these will be inputs (`project_id`, `dataset_name` and `call_set_identifier`) to all or most of the GVS workflows. The [naming conventions for other aspects of GVS datasets are outlined here](https://docs.google.com/document/d/1pNtuv7uDoiOFPbwe4zx5sAGH7MyxwKqXkyrpNmBxeow).

## The Pipeline
1. `GvsAssignIds` workflow
   - To optimize the GVS internal queries, each sample must have a unique and consecutive integer ID assigned. Running the `GvsAssignIds` will create a unique GVS ID for each sample (`sample_id`) and update the BQ `sample_info` table (creating it if it doesn't exist). This workflow takes care of creating the BQ `vet_*`, `ref_ranges_*` and `cost_observability` tables needed for the sample IDs generated.
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
4. **TBD Workflow to soft delete samples**
5. `GvsPopulateAltAllele` workflow
   - **TODO:** needs to be made cumulative so that it can add data to the existing table instead of creating it from scratch on each run (see [VS-52](https://broadworkbench.atlassian.net/browse/VS-52))
   - This step loads data into the `alt_allele` table from the `vet_*` tables in preparation for running the filtering step.
   - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
6. `GvsCreateFilterSet` workflow
   - This step calculates features from the `alt_allele` table, and trains the VQSR filtering model along with site-level QC filters and loads them into BigQuery into a series of `filter_set_*` tables.
   - See [naming conventions doc](https://docs.google.com/document/d/1pNtuv7uDoiOFPbwe4zx5sAGH7MyxwKqXkyrpNmBxeow) for guidance on what to use for `filter_set_name`, which you will need to keep track of for the `GvsExtractAvroFilesForHail` WDL. If, for some reason, this step needs to be run multiple times, be sure to use a different `filter_set_name` (the doc has guidance for this, as well).
   - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
7. `GvsExtractAvroFilesForHail` workflow
   - This workflow extracts the data in BigQuery and transforms it into Avro files in a Google bucket, incorporating the VQSR filter set data.
   - The extracted Avro files will then be used as an input for the python script `hail_gvs_import.py` which will then be used to produce a Hail VDS.
   - This workflow needs to be run with the `filter_set_name` from `GvsCreateFilterSet` step.
   - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
8. Run the VDS Extract using a notebook terminal and the python script called `hail_gvs_import.py` from the GenerateHailScripts task in the `GvsExtractAvroFilesForHail` workflow. It will look something like: `gs://fc-<workspace-id>/submissions/<submission id>/GvsExtractAvroFilesForHail/<workflow id>/call-GenerateHailScripts/hail_gvs_import.py`
   - This step creates a VDS based on the Avro files 
   - Notebook provisioning suggestions and other necessary set up details (e.g. the Hail wheel to use) can be found in the file: [AoU Delta VDS Cluster Configuration.md](vds/cluster/AoU%20Delta%20VDS%20Cluster%20Configuration.md)
   - Avro files from the above `GvsExtractAvroFilesForHail` step can be found in the avro directory in the OutputPath task: `gs://fc-<workspace-id>/submissions/<submission id>/GvsExtractAvroFilesForHail/<workflow id>/call-OutputPath/avro`
   - We suggest gsutil cp-ing `hail_gvs_import.py` to the notebook and then invoking it directly in the terminal
   - inputs to the python script are:
     1. `--avro-path`: the directory path at which exported GVS Avro files are found in GCP.  This path can be found by navigating to the execution directory of the `OutputPath` task in the `GvsExtractAvroFilesForHail ` workflow.  The files are in a directory there named "avro".
     2. `--vds-path`: the desired output path to which the VDS should be written
     3. `--temp-path`: a convenient path to temporary directory. We suggest a folder under the GCP bucket of the workspace that the notebook is in, e.g. `gs://fc-<workspace-id>/hail_tmp`.
9. `GvsPrepareRangesCallset` workflow
    - This workflow transforms the data in the vet tables into a schema optimized for callset stats creation and for calculating sensitivity and precision.
    - It will need to be run twice. Once with `only_output_vet_tables` set to "true" (the default value is `false`) and once with `control_samples` set to "true" (the default value is `false`).
    - See [naming conventions doc](https://docs.google.com/document/d/1pNtuv7uDoiOFPbwe4zx5sAGH7MyxwKqXkyrpNmBxeow) for guidance on what to use for `extract_table_prefix` or cohort prefix, which you will need to keep track of for the callset stats. Use the same prefix for both runs.
    - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
8. `GvsExtractCallset` workflow
    - This workflow extracts the control data in BigQuery and transforms it into a sharded joint VCF incorporating the VQSR filter set data.
    - It will need to be run once with `control_samples` set to "true", and with the `filter_set_name` from the `GvsCreateFilterSet` step and the `extract_table_prefix` from the `GvsPrepareRangesCallset` step.  
    - Include a valid `output_gcs_dir` parameter, which is where the VCF, interval list, manifest, and sample name list files will go, e.g. `gs://fc-<workspace-id>/extractedcontrol-vcfs/`
    - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
10. `GvsCallsetStatistics` workflow
    - You will need to have the "BigQuery Data Viewer" role for your @pmi-ops proxy group on the `spec-ops-aou:gvs_public_reference_data.gnomad_v3_sites` table
    - This workflow needs to be run with the `extract_table_prefix` input from `GvsPrepareRangesCallset` step.
    - This workflow needs to be run with the `filter_set_name` input from `GvsCreateFilterSet` step.
    - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
11. `GvsCalculatePrecisionAndSensitivity` workflow
    - You will need to have "Storage Object View" access granted for your @pmi-ops proxy group on the `gs://broad-dsp-spec-ops/gvs/truth` directory
    - This workflow needs to be run with the control sample chr20 vcfs from `GvsExtractCallset` step which were placed in the `output_gcs_dir`.
    - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
12. `GvsCallsetCost` workflow
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



