# Sample set QC

Sample set qc is performed after a callset has been extracted to validate the quality of the callset with respect to the VQSR filter set used. Most of these steps are manual at this point and are executed in the BQ console. 

## Create a `vet_all` view
We need to create a view that contains all the vet data. Run

	bq mk --project_id=<project> --use_legacy_sql=false --view <query> <dataset>.<prefix>_vet_all

where `<query>` is of the form:

	"SELECT *
     FROM (
        SELECT * FROM <dataset>.vet_001 UNION ALL
        SELECT * FROM <dataset>.vet_002 UNION ALL
        SELECT * FROM <dataset>.vet_003) unioned_vets
     WHERE
        unioned_vets.sample_id IN
        (SELECT sample_id FROM <dataset>.<prefix>__SAMPLES)"

and `prefix` is the "cohort_extract_table_prefix" from `GvsExtractCallset` and there is a `SELECT * FROM <dataset>.vet_` in the `<query>` for every `vet_` table from the sample set.

## Create the sample_metrics table

- Copy the contents of `filter_set_samples_create_metrics.example.sql` to the BQ console
- Replace `$FQ_DATASET` with `<project>.<dataset>` 
- Replace `$FQ_VET_ALL` with `<project>.<dataset>.<prefix>_vet_all`
- Replace `$NAME_OF_FILTER_SET` with the name of the filter set you used to filter the callset
- Execute the query

This will currently create/overwrite the sample_metrics table in the dataset. (TODO we may want this to just update the table so we can support multiple filter sets)

## Evaluate and extract the qc data

- Copy the contents of `filter_set_samples_calculate_threshold.example.sql` to the BQ console
- Replace `$FQ_DATASET` with `<project>.<dataset>`
- Replace `$NAME_OF_FILTER_SET` with the name of the filter set you used to filter the callset
- Execute the query
- Click save results to export the results

