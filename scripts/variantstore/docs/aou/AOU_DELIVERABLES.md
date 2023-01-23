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
  - [GvsExtractAvroFilesForHail](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsExtractAvroFilesForHail) workflow (Hail VDS output)
  - [GvsPrepareRangesCallset](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsPrepareRangesCallset) workflow (Hail VDS output -- used for callset stats)
  - [GvsCalculatePrecisionAndSensitivity](https://dockstore.org/workflows/github.com/broadinstitute/gatk/GvsCalculatePrecisionAndSensitivity) workflow
  - [GvsCallsetCost](https://dockstore.org/workflows/github.com/broadinstitute/gatk/GvsCallsetCost) workflow
- Run the "Fetch WGS metadata for samples from list" notebook after you have placed the file with the list of the new samples to ingest in a GCS location the notebook (running with your @pmi-ops account) will have access to.  This will grab the samples from the workspace where they were reblocked and bring them into this callset workspace.
  - Set the `sample_list_file_path` variable in that notebook to the path of the file
  - Run the "now that the data have been copied, you can make sample sets if you wish" step if you want to automatically break up the new samples into smaller sample sets.  Set the `SUBSET_SIZE` and `set_name` variables to customize.
- **NOTE** If you want to create a large sample set after you have run the notebook, Terra provides (and recommends you use) this python [script](https://github.com/broadinstitute/firecloud-tools/tree/master/scripts/import_large_tsv) which allows you to upload a sample set to the workspace.
- **TBD This is based on VCF output; VDS output might need different increases.** For extracting VCFs as the final output for the callset, you will want to increase the Google quotas for the workspace project (you can find this in the workspace dashboard under Cloud Information > Google Project ID) to these levels (all in the workspace region):
  - Persistent Disk Standard (GB): 1,000,000 GB (1 PB)
  - CPUs: 6,400
  - In-use IP addresses: 5,000 (this is the most challenging one and will probably require contacting the GCP account team to facilitate)
  - VM instances: 64,000
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
   - See [naming conventions doc](https://docs.google.com/document/d/1pNtuv7uDoiOFPbwe4zx5sAGH7MyxwKqXkyrpNmBxeow) for guidance on what to use for `filter_set_name`, which you will need to keep track of for the `GvsExtractCallset` WDL. If, for some reason, this step needs to be run multiple times, be sure to use a different `filter_set_name` (the doc has guidance for this, as well).
   - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
7. `GvsExtractAvroFilesForHail` workflow (Hail VDS output)    
   - This workflow extracts the data in BigQuery and transforms it into Avro files in a Google bucket, incorporating the VQSR filter set data.  We suggest this method for callsets of more than 100k samples.
   - The extracted Avro files will then be used as the input for the below notebook which can then be run to produce a Hail VDS
   - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
8. Run the VDS Extract using a notebook terminal and the python script `hail_gvs_import.py`
   - This step creates a VDS based on the Avro files 
   - Notebook provisioning suggestions can be found in the file: [AoU Delta VDS Cluster Configuration.md](vds/cluster/AoU Delta VDS Cluster Configuration.md)
   - `hail_gvs_import.py` is the python script that needs to be run in the notebook-- we suggest running it directly in the terminal
   - inputs to the python script are:
     1. `--avro-path`: the directory path at which exported GVS Avro files are found in GCP
     2. `--vds-path`: the desired output path to which the VDS should be written
     3. `--temp-path`: a convenient path to temporary directory (we suggest the GCP bucket of the workspace that the notebook is in) 
9. `GvsPrepareRangesCallset` workflow (Hail VDS output -- for callset stats)
   - This workflow transforms the data in the vet tables into a schema optimized for VCF generation during the Extract step.
   - It will need to be run once with `only_output_vet_tables` set to "true" (see [naming conventions doc](https://docs.google.com/document/d/1pNtuv7uDoiOFPbwe4zx5sAGH7MyxwKqXkyrpNmBxeow) for guidance on what to use for `extract_table_prefix` or cohort prefix, which you will need to keep track of for the callset stats); the default value is `false`.
   - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
10. Run the notebook to create callset stats (**TODO**, see [VS-388](https://broadworkbench.atlassian.net/browse/VS-388)), for which you need
   - "BigQuery Data Viewer" role for your @pmi-ops proxy group on the `spec-ops-aou:gvs_public_reference_data.gnomad_v3_sites` table
   - the Google project ID you used for all the GVS WDLs (`project_id` input)
   - the name of the BigQuery dataset you used for all the GVS WDLs (`dataset_name` input)
   - the `extract_table_prefix` input from `GvsExtractCallset` step
   - the `filter_set_name` input from `GvsCreateFilterSet` step
11. [GvsCalculatePrecisionAndSensitivity](../../tieout/AoU_PRECISION_SENSITIVITY.md#generating-callset-precision-and-sensitivity-values) workflow, for which you need
   - "Storage Object View" access granted for your @pmi-ops proxy group on the `gs://broad-dsp-spec-ops/gvs/truth` directory
12. `GvsCallsetCost`
   - The cost from this callset, which represents the total BigQuery cost (which is not represented in the Terra UI total workflow cost) from the GVS pipeline workflows, is used to calculate the cost of the callset as a whole and by sample.


## Deliverables (via email to stakeholders once the above steps are complete)
1. GCS locations of the VCFs, indexes and interval_list files (subpaths of the `output_gcs_dir` input from `GvsExtractCallset`)
2. fully qualified name of the BigQuery dataset (composed of the `project_id` and `dataset_name` inputs from the workflows)
3. callset statistics CSV file (see step #11 above)
4. TSV output from `GvsCalculatePrecisionAndSensitivity` workflow (see step #12 above)

## Running the VAT pipeline
To create a BigQuery table of variant annotations, you may follow the instructions here:
[process to create variant annotations table](../../variant_annotations_table/README.md)

The pipeline takes in jointVCF shards and outputs a variant annotations table in BigQuery.



