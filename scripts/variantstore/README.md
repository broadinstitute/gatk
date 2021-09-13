# Running the Genome Variant Store (GVS) Pipeline for AoU

## Steps
1. If this is the first time running the pipeline in a particular Google billing project, submit a ticket requesting the following changes to enable read API clustering:
   - TBD
2. If you haven't already, clone the [QuickStart Workspace](https://app.terra.bio/#workspaces/broad-dsde-firecloud-billing/Genomic%20Variant%20Store%20-%20GVS%20Quickstart), selecting the AoU authorization domain.

[comment]: <> (The genome center will deliver the manifests and load the samples into the datamodel in a Terra workspace and create sample sets.  Anyhing up to a certain date goes into a release AoU_DRCV_SEQ_2021_08_06)

[comment]: <> (Run the re-blocking workflow ReblockGVCF on the sample &#40;step 1&#41; and then select the sample sets &#40;step 2&#41;)

[comment]: <> (allofus-drc-wgs-dev-prod/AoU_DRC_WGS is the place where AoU samples are getting loaded/re-blocked)

[comment]: <> (do we need to support case where path to VCF file is null after re-blocking?)

[comment]: <> (why does re-blocking fail?)

[comment]: <> (- the path to the file&#40;s&#41; is wrong)

[comment]: <> (- because index that's provided with the cram needs to be re-generated &#40;we add the GC to do this&#41;)

[comment]: <> (- there may be an error in the re-blocking WF &#40;Laura&#41;)

3. Re-block gVCFs and put into cloned Quickstart workspace.
4. In that cloned workspace, run the workflows outlined in the Workspace description; using the same “data_project” and “default_dataset” throughout
   1. GvsAssignIds at the sample set level
   2. GvsImportGenomes at the sample set level
   3. GvsCreateAltAllele
   4. GvsCreateFilterSet
      - make note of the “filter_set_name” input used
   5. GvsPrepareCallset
      - if the callset has more than 20K genomes, you will want to reserve flex slots for this step and then release them once it has completed
   6. GvsExtractCallset
       - include a valid (and secure) “output_gcs_dir” parameter, this is where the VCF files will go
       - after `SplitIntervals` task has completed, make note of location of interval_list files (e.g. “fc-secure-…./…/GvsExtractCallset/…/call-SplitIntervals/.../*.interval_list”
       - make note of the “cohort_extract_table_prefix”
5. Run [callset QC](callset_QC/README.md), for which you need
   1. permission to query table `spec-ops-aou:gvs_public_reference_data.gnomad_v3_sites`
   2. the data_project you used for all the Quickstart WDLs
   3. the default_dataset you used for all the Quickstart WDLs
   4. the “cohort_extract_table_prefix” input from `GvsExtractCallset` step
   5. the “filter_set_name” input from `GvsCreateFilterSet` step

## Deliverables
1. location of the VCF files (“output_gcs_dir” input from GvsExtractCallset)
2. location of the interval_list files (output of `SplitIntervals` task in `GvsExtractCallset`)
3. results of running callset QC

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
