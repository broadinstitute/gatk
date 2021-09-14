# Running the Genome Variant Store (GVS) Pipeline for AoU

## Steps
1. If this is the first time running the GVS pipeline in a particular Google billing project, submit a ticket requesting the following changes to enable read API clustering:
   - TBD
2. This list assumes you are starting off with a Terra workspace that contains the re-blocked gVCFs for your callset as "sample sets" in the "Data" tab.  For workflows, you will be following the steps outlined in the [GVS Quickstart workspace](https://app.terra.bio/#workspaces/broad-dsde-firecloud-billing/Genomic%20Variant%20Store%20-%20GVS%20Quickstart).  Make special note of the "Prerequisites" section.
3. Run the workflows outlined in the Quickstart description in the workspace that contains the data by using the "Copy to Another Workspace" function.  The first two steps (GvsAssignIds and GvsImportGenomes) can be run multiple times to load all necessary sample sets into one "instance" of GVS.  The rest of the steps are only run once, in order, and only after all the samples for the callset are loaded.
   1. GvsAssignIds at the sample set level ("Step 2")
   2. GvsImportGenomes at the sample set level ("Step 2")
   3. GvsCreateAltAllele
   4. GvsCreateFilterSet
      - make note of the “filter_set_name” input used
   5. GvsPrepareCallset
      - if the callset has more than 20K genomes, you will want to reserve GCP flex slots for this step and then release them once it has completed
   6. GvsExtractCallset
       - include a valid (and secure) “output_gcs_dir” parameter, this is where the VCF and interval list files  will go
       - make note of the “cohort_extract_table_prefix” input that you enter
   7. Run [callset QC](callset_QC/README.md), for which you need
      1. permission to query table `spec-ops-aou:gvs_public_reference_data.gnomad_v3_sites`
      2. the data_project you used for all the GVS WDLs
      3. the default_dataset you used for all the GVS WDLs
      4. the “cohort_extract_table_prefix” input from `GvsExtractCallset` step
      5. the “filter_set_name” input from `GvsCreateFilterSet` step

## Deliverables
1. location of the VCF and interval_list files (“output_gcs_dir” input from GvsExtractCallset)
2. results of running callset QC

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
