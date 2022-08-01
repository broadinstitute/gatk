# Running the Genome Variant Store (GVS) Pipeline for AoU

## Setup
- Create a Terra workspace with your @pmi-ops.org account with the `allofus-drc-wgs-dev` Terra Billing account in the `AOU_DRC_WGS_OPERATORS` Authorization domain. Share it with `dsp-variant-team@firecloud.org` so the rest of the team can see it (Reader) and with the individual users who will be running the workflows (Owner).
- Populate the workspace with the following:
  - [Fetch WGS metadata for samples from list](http://app.terra.bio/#workspaces/allofus-drc-wgs-dev/GVS%20AoU%20WGS%20Charlie/notebooks/launch/Fetch%20WGS%20metadata%20for%20samples%20from%20list.ipynb) notebook
  - GvsAssignIds workflow
  - GvsImportGenomes workflow
  - GvsWithdrawSamples workflow
  - GvsCreateAltAllele workflow
  - GvsCreateFilterSet workflow
  - GvsPrepareRangesCallset workflow
  - GvsExtractCallset workflow
- Run the "Fetch WGS metadata for samples from list" notebook after you have placed the file with the list of the **new** samples to ingest in a place the notebook will have access to.  This will grab the samples from the workspace where they were reblocked and bring them into this callset workspace.
  - Set the `sample_list_file_path` variable in that notebook to the path of the file
  - Run the "now that the data have been copied, you can make sample sets if you wish" step if you want to automatically break up the new samples into smaller sample sets.  Just be sure to set the `SUBSET_SIZE` and `set_name` variables beforehand.
- You will want to increase the Google quotas for the workspace project (you can find this in the workspace dashboard under Cloud Information > Google Project ID) to these levels (all in region `us-central1`):
    1. Persistent Disk Standard (GB): 1,000,000 GB (1 PB)
    2. CPUs: 64,000
    3. In-use IP addresses: 5,000 (this is the most challenging one and will probably require contacting the GCP account team to facilitate)
    4. VM instances: 64,000
- Make a note of the Google project id (`aou-genomics-curation-prod`), dataset name (`aou_wgs`) and callset indentifier (e.g. "Bravo") as these will be inputs to all the GVS workflows. The [naming conventions for other aspects of GVS datasets are outlined here](https://docs.google.com/document/d/1pNtuv7uDoiOFPbwe4zx5sAGH7MyxwKqXkyrpNmBxeow).

## Workflows
1. `GvsAssignIds` at the `sample set` level ("Step 1" in workflow submission) with a sample set of all the new samples to be included in the callset (if new controls are being added, they need to be done in a separate run, with the `samples_are_controls` input set to "true")
2. `GvsImportGenomes` at the `sample set` level ("Step 1" in workflow submission) with sample sets of the new samples (currently we divide them up into sets of max 4000 each, but [we want to raise or remove this limit](https://broadworkbench.atlassian.net/browse/VS-344)).
3. `GvsWithdrawSamples` if there are any samples to withdraw from the last callset.
4. GvsCreateAltAllele
5. GvsCreateFilterSet (see [naming conventions doc](https://docs.google.com/document/d/1pNtuv7uDoiOFPbwe4zx5sAGH7MyxwKqXkyrpNmBxeow) for guidance on what to name the filter set, which you will need to keep track of for the `GvsExtractCallset` WDL).
6. GvsPrepareRangesCallset needs to be run twice, once with `control_samples` set to "true" (see [naming conventions doc](https://docs.google.com/document/d/1pNtuv7uDoiOFPbwe4zx5sAGH7MyxwKqXkyrpNmBxeow) for guidance on what to use for `extract_table_prefix` or cohort prefix, which you will need to keep track of for the `GvsExtractCallset` WDL).
7. GvsExtractCallset needs to be run twice, once with `control_samples` set to "true", and with the `filter_set_name` and `extract_table_prefix` from step 5 & 6.  Include a valid (and secure) "output_gcs_dir" parameter, which is where the VCF and interval list files  will go.
8. Run [process to create callset statistics](callset_QC/README.md), for which you need
    1. permission to query table `spec-ops-aou:gvs_public_reference_data.gnomad_v3_sites`
    2. the data_project you used for all the GVS WDLs
    3. the default_dataset you used for all the GVS WDLs
    4. the `extract_table_prefix` input from `GvsExtractCallset` step
    5. the `filter_set_name` input from `GvsCreateFilterSet` step

## Deliverables (via email once the above steps are complete)
1. location of the VCFs and interval_list files (`output_gcs_dir` input from GvsExtractCallset)
2. fully qualified name of the BigQuery dataset (input `dataset_name` in the workflows)
3. [callset statistics](callset_QC/README.md) CSV
4. [precision and sensitivity results](tieout/AoU_PRECISION_SENSITIVITY.md) TSV

## Running the VAT pipeline
To create a BigQuery table of variant annotations, you may follow the instructions here:
[process to create variant annotations table](variant_annotations_table/README.md)

The pipeline takes in a jointVCF and outputs a variant annotations table in BigQuery.



