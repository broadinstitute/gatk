


## Cost estimation for BQ queries

### Using the AoU Service account
To get cost information for the AoU queries, only the AoU Service account has access. 
After using the json file to auth as that service account, you can query the 
`region-us`.INFORMATION_SCHEMA.JOBS_BY_USER table to narrow down the query you are interested in.

    gcloud auth activate-service-account --key-file=<keyfile>
    # use queries like this to narrow down the time window of the queries you are interested in
    bq --project_id <project> query --use_legacy_sql=false 'select total_bytes_billed, job_id, creation_time, start_time, query from `region-us`.INFORMATION_SCHEMA.JOBS_BY_USER where creation_time > "2021-04-16 00:00:00" and creation_time < "2021-04-18 00:00:00" and total_bytes_billed > 1 order by creation_time asc limit 100'
    
    # use a query like this to get the total number of bytes billed
    bq --project_id <project> query --use_legacy_sql=false 'select sum(total_bytes_billed) as total from `region-us`.INFORMATION_SCHEMA.JOBS_BY_USER where creation_time > "2021-04-16 00:00:00" and creation_time < "2021-04-18 00:00:00" and total_bytes_billed > 1 limit 100'
