## Testing WDLs using Service Account authentication

Some data (including All of Us data) will be processed using a service account (SA) that has access to the relevant data (GCP storage) and BigQuery (BQ) resources.

Our WLDs are written so that they can be run EITHER (1) as the user (no `service_account_json` file provided as input) OR (2) as the SA (by providing a `service_account_json` input file to the WDL).

## Manual steps to create testing resources
To test running WDLs using a SA, we have created a project, SA, buckets, and BQ datasets that mirror the permissions in the AoU setup. These include:
* Project `specops-variantstore-sa-tests`
    * Owners: `marymorg@broadinstitute.org`, `kcibul@broadinstitute.org`
* SA `variant-store-test-sa@specops-variantstore-sa-tests.iam.gserviceaccount.com`
* SA key file at `gs://specops-variantstore-sa-tests-bucket/sa/variant-store-test-sa_key.json`
* Bucket to house data `gs://variantstore-sa-tests-data-bucket`
    * note that data can live here but sample map needs to be in the workspace bucket so it can be localized to cromwell
    * gvcfs and index files for 100 samples from the 1000G dataset are loaded into a directory `1000g_gvcfs`
    * data table for these files is in the workspace: https://app.terra.bio/#workspaces/broad-dsp-spec-ops-fc/gvs_sa_testing
* BQ datasets:
    * `variantstore_sa_tests` - data copied from `spec-ops-aou:anvil_100_for_testing` - (to test extract processes)
    * `variantstore_sa_tests_filtering` - same as previous, but does not have filtering tables (to test their generation)
    * `temp_tables` (part of extract processes)

## Permissions
* Permissions for SA:
    * Project-level permissions:
        * `BigQuery User`
        * `BigQuery Data Owner`
    * Bucket-level permissions:
        * `Storage Object Admin` on `gs://variantstore-sa-tests-data-bucket`
    * BQ dataset-level permissions are NOT required with project-level permissions described above.
* Permissions for user (and Terra pet, or use Terra proxy group):
    * `Storage Object Viewer` permissions to bucket containing SA key file (`gs://specops-variantstore-sa-tests-bucket`)
    
