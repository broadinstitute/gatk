# Running the Genome Variant Store (GVS) Pipeline for AoU

## Steps
1. If you haven't already, clone the [QuickStart Workspace](https://app.terra.bio/#workspaces/broad-dsde-firecloud-billing/Genomic%20Variant%20Store%20-%20GVS%20Quickstart)
2. In that cloned workspace, run the workflows outlined in the Workspace description; using the same “data_project” and “default_dataset” throughout
   1. GvsAssignIds
   2. GvsImportGenomes
   3. GvsCreateAltAllele
   4. GvsCreateFilterSet
      - make note of the “filter_set_name” input used
   5. GvsPrepareCallset
   6. GvsExtractCallset
       - include a valid (and secure) “output_gcs_dir” parameter, this is where the VCF files will go
       - after `SplitIntervals` task has completed, make note of location of interval_list files (e.g. “fc-secure-…./…/GvsExtractCallset/…/call-SplitIntervals/.../*.interval_list”
       - make note of the “cohort_extract_table_prefix”

## Deliverables
1. location of the VCF files (“output_gcs_dir” input from GvsExtractCallset)
2. location of the interval_list files (output of `SplitIntervals` task in `GvsExtractCallset`)
3. results of [running QC](callset_QC/README.md), for which you need
   1. permission to query table `spec-ops-aou:gvs_public_reference_data.gnomad_v3_sites`
   2. the data_project you used for all the Quickstart WDLs
   3. the default_dataset you used for all the Quickstart WDLs
   4. the “cohort_extract_table_prefix” input from `GvsExtractCallset` step
   5. the “filter_set_name” input from `GvsCreateFilterSet` step

## Cost estimation for BQ queries
To get cost information for the AoU queries, only the AoU Service account has access. 
After using the json file to auth as that service account, you can query the 
`region-us`.INFORMATION_SCHEMA.JOBS_BY_USER table to narrow down the query you are interested in.

    gcloud auth activate-service-account --key-file=<keyfile>
    # use queries like this to narrow down the time window of the queries you are interested in
    bq --project_id <project> query --use_legacy_sql=false 'select total_bytes_billed, job_id, creation_time, start_time, query from `region-us`.INFORMATION_SCHEMA.JOBS_BY_USER where creation_time > "2021-04-16 00:00:00" and creation_time < "2021-04-18 00:00:00" and total_bytes_billed > 1 order by creation_time asc limit 100'
    
    # use a query like this to get the total number of bytes billed
    bq --project_id <project> query --use_legacy_sql=false 'select sum(total_bytes_billed) as total from `region-us`.INFORMATION_SCHEMA.JOBS_BY_USER where creation_time > "2021-04-16 00:00:00" and creation_time < "2021-04-18 00:00:00" and total_bytes_billed > 1 limit 100'


## Running the VAT pipeline
The VAT pipeline is a set of WDLs
 - GvsSitesOnlyVCF.wdl
 - GvsValidateVAT.wdl

The pipeline takes in a jointVCF and outputs a table in BigQuery.
