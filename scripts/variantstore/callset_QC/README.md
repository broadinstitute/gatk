# Callset Stats

The generation of GVS callset stats is performed after a callset has been extracted to validate the quality of the callset with respect to the VQSR filter set used. Most of these steps are manual at this point and are executed in the BQ console.

**NOTE:** You will need permission to query table `spec-ops-aou:gvs_public_reference_data.gnomad_v3_sites`.

## Use the prepared VET_DATA table
We need a table that contains all the vet data. Use:

	<dataset>.<prefix>__VET_DATA

and `prefix` is the "cohort_extract_table_prefix" from `GvsExtractCallset` 
Note that this table is partitioned by position / location

## Create the sample_metrics table

- Copy the contents of `filter_set_samples_create_metrics.example.sql` to the BQ console
- Replace `$FQ_PREFIX` with `<project>.<dataset>.<prefix>`
- Replace `$NAME_OF_FILTER_SET` with the name of the filter set you used to filter the callset
- Execute the query

This will currently create/overwrite the sample_metrics table in the dataset. (TODO we may want this to just update the table so we can support multiple filter sets)

## Evaluate and extract the data

- Copy the contents of `filter_set_samples_calculate_threshold.example.sql` to the BQ console
- Replace `$FQ_PREFIX` with `<project>.<dataset>.<prefix>`
- Replace `$NAME_OF_FILTER_SET` with the name of the filter set you used to filter the callset
- Execute the query
- Click the SAVE RESULTS button to export the results, or you can save the query as an SQL file and then run this command locally to get a CSV file in your local environment:

	`bq --project_id=<project> --use_legacy_sql=false query --format=csv --use_legacy_sql=false "$(cat <path to SQL file>)" > <path to destination CSV file>`

	**NOTE:** using the size of the callset (e.g. 22K) as an identifier will probably be least confusing way to distinguish it from other callset data (as opposed to the BigQuery prefix and/or filter set name referenced above).
