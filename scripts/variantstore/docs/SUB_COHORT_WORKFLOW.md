# How to Run the GvsExtractCohortFromSampleNames Workflow

The purpose of `GvsExtractCohortFromSampleNames.wdl` is to take advantage of an existing GVS (Genomic Variant Store) in BigQuery, complete with filter model, to generate a callset with a subset of specified samples' data.  It calls existing WDLs for the "Prepare" and "Extract" steps and allows for the data required to create the subset callset to be stored and queried (and, therefore, paid for) by a different Google project than the "parent" GVS (most probably paid for by AoU).

Required inputs:

| Parameter              | Description                                                                                                                                                                                                                              |
|------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| call_set_identifier    | a unique name for this cohort (used for cost tracking)                                                                                                                                                                                   |
 | cohort_table_prefix    | a unique name for this cohort (can be the same as `call_set_identifier`) should only contain letters and underscores                                                                                                                     |
 | destination_dataset_name | name of the BigQuery dataset that's in the destination Google project ID                                                                                                                                                                 |
 | destination_project_id | Google project ID where the `destination_dataset_name` lives                                                                                                                                                                             |
 | extraction_uuid | unique name used as a query label for the "Prepare" step, can be the same as call_set_identifier and/or cohort_table_prefix (requested as an input by Verily for cost-trackings; if this isn't Verily-related, the value doesn't matter) |
 | filter_set_name | the same input that was used for `GvsCreateFilterSet` (or `GvsJointVariantCalling`)                                                                                                                                                      |
 | gvs_dataset | the same input value that was used for `dataset_name` in the GVS step WDLs (or `GvsJointVariantCalling`)                                                                                                                                 |
 | gvs_project | the same input value that was used for `project_id` in the GVS step WDLs (or `GvsJointVariantCalling`)                                                                                                                                   |
 | output_file_base_name | the base file name for the VCFs that will be created in the "Extract" step                                                                                                                                                               |
 | query_project | Google project ID for the permissions/billing of the "Prepare" and "Extract" steps for this sub-cohort (can be the same as `destination_project_id`)                                                                                     |
 | scatter_count | Number of shards to scatter "Extract" step across                                                                                                                                                                                        |
 | control_samples | GCS path to a file that contains a list of the samples (`sample_name` field in GVS) to include in the sub-cohort; if not set, `cohort_sample_names_array` must be set                                                                    |
 | cohort_sample_names_array | Array of sample identifiers (`sample_name` field in GVS) to include in the sub-cohort; if not set, `control_samples` must be set                                                                                                         |

Optional inputs of interest:

| Parameter              | Description                              |
|------------------------|--------------------------------------------|
 | drop_state     | This should correspond to the same value set in `GvsImportGenomes` (or `GvsJointVariantCalling`)     |
 | output_gcs_dir | GCS path to a directory to copy the interval list files, the extract VCFs and a sample manifest into |
