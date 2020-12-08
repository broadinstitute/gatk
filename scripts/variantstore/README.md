
# NOTE: this is a work in progress, please send feedback to dsp-spec-ops@broadinstitute.org

# Before you start

Make sure you have the following ready:

- A collection of HG38 reblocked GVCFs
- A google project containing
  - A BigQuery dataset to hold the genomes
  - A BigQuery dataset for temporary storage (can be same as above, or can be a separate dataset with different permissions or table expiration defaults)
  - A Google Cloud Storage bucket
- A cromwell instance running as a user with access to the BigQuery dataset and the google buckets
- A sample mapping file (see below)
- Python 3.6+ installed with the bigquery client library installed (`google-cloud-bigquery`)

Note:  Where a parameter is prefixed with `fq` it is a fully qualified name (e.g. `project.dataset` or `project.dataset.table`)

# QuickStart

## Overview

There are 3 primary steps use cases in the Variant Store:

1. Importing Genomes
2. Training a Filtering Model
3. Extracting a Cohort with filtering applied

The workflow for each of these is encapsulated in a Cromwell WDL workflow.  These workflows can be found at

 `https://github.com/broadinstitute/gatk/tree/ah_var_store/scripts/variantstore/wdl`

## Importing Genomes

gVCFs from genomes are imported through the `ImportGenomes.wdl` WDL.  An example inputs json for this workflow can be found at `ImportGenomes.example.inputs.json`

Most of the parameters in the WDL are constant, such as the reference genome.  However, several should be changed for each run.  These are:

| Parameter Name| Description |
| ------------- | ------------- |
| ImportGenomes.project_id  | The name of the google project containing the BigQuery Dataset  |
| ImportGenomes.dataset_name  | The name of the BigQuery dataset  |
| ImportGenomes.sample_map | The GCS Path to the sample map file |
| ImportGenomes.input_vcfs| An array of the reblocked gVCF files to be loaded |
| ImportGenomes.output_directory | A GCS path to be used for temporary storage of data while loading |

### Application Notes
- After a successful import the temporary directory specified by `output_directory` can be deleted.  It is not deleted by default currently by the WDL



## Training Filtering Model

A powerful variant filtering model is trained and uploaded into BigQuery through the the `ngs_filter_extract.wdl` WDL.  An example inputs json for this workflow can be found at `ngs_filter_extract.inputs.json`

Most of the parameters in the WDL are constant, such as the reference genome.  However, several should be changed for each run.  These are:

| Parameter Name| Description |
| ------------- | ------------- |
| NgsFilterExtract.fq_sample_table | Specifies which samples to use for training, to use all samples use the `metadata` table |
| NgsFilterExtract.fq_alt_allele_table | FQ name of the alt allele table, typically `<project>.<dataset>.ALT_ALLELE` (DEVNOTE: should make optional) |
| NgsFilterExtract.fq_filter_set_info_table | FQ name of the alt allele table, typically `<project>.<dataset>.filter_set_info` (DEVNOTE: should make optional) |
| NgsFilterExtract.fq_filter_set_tranches_table | FQ name of the alt allele table, typically `<project>.<dataset>.filter_set_tranches` (DEVNOTE: should make optional) |
| NgsFilterExtract.query_project           | Google project where queries are executed, for billing purposes |
| NgsFilterExtract.filter_set_name | Dataset-wide unique name used to identify this training run (e.g. `wdl_testing-1234`) |
  NgsFilterExtract.output_file_base_name | Base name for temporary files in the dataset (DEVNOTE: should make optional) |

  
## Extracting Cohort

Extracting a cohort currently consists of two steps

1. Prepare - A python script to transform the data in BigQuery into a staging table for extract by the WDL
2. Extract - A WDL that scatters and extracts the cohort into sharded VCFs

Eventually this python script will be turned into it's own WDL or incorporated into the extract WDL.

### Prepare

The python script is located in `https://github.com/broadinstitute/gatk/tree/ah_var_store/scripts/variantstore/wdl/extract/ngs_cohort_extract.py`

This script requires the following parameters:

| Parameter Name| Description |
| ----------------------- | ------------- |
| fq_petvet_dataset       | Dataset containing the data to be extracted|
| fq_temp_table_dataset   | Dataset for temp tables, described in prereqs  |
| fq_destination_dataset  | Dataset for prepared cohort, can be same as temp table dataset |
| destination_table       | Unique name for the destination extract table |
| fq_cohort_sample_names  | Table containing list of samples to be extracted.  Has a single column `sample_id` with the numeric id of the each sample from the sample map |
| min_variant_samples     | (Optional) minimum number of variant samples at each site to be extracted |
| query_project           | Google project where queries are executed, for billing purposes |
| fq_sample_mapping_table | Complete sample mapping table, typically the `metadata` table in the dataset |


An example invocation is as follows
```
PROJECT="my-example-project"
DATASET="bq_wgs_example"

python ngs_cohort_extract.py \
  --fq_petvet_dataset ${PROJECT}.${DATASET} \
  --fq_temp_table_dataset ${PROJECT}.temp_tables \
  --fq_destination_dataset ${PROJECT}.${DATASET} \
  --destination_table exported_cohort_100_test \
  --fq_cohort_sample_names ${PROJECT}.${DATASET}.cohort_10_of_35 \
  --min_variant_samples 0 \
  --query_project ${PROJECT} \
  --fq_sample_mapping_table ${PROJECT}.${DATASET}.metadata
```



### Extract

Once the prepare step has been completed, the `ngs_cohort_extract.wdl` WDL.  An example inputs json for this workflow can be found at `ngs_cohort_extract.inputs.json`

Most of the parameters in the WDL are constant, such as the reference genome.  However, several should be changed for each run.  These are:

| Parameter Name| Description |
| ------------- | ------------- |
| NgsCohortExtract.fq_sample_table | Table of samples to be extracted, same as `fq_cohort_sample_names` from prepare
| NgsCohortExtract.fq_cohort_extract_table | Fully qualified name of the extract table, comprised of the `fq_destination_dataset` and `destination_table` from the prepare script|
| NgsCohortExtract.fq_filter_set_table | Name of the filter set table to use, from the filter model training |
| NgsCohortExtract.filter_set_name | Name of the filter set |
| NgsCohortExtract.query_project | Google project used to execute queries, for billing purposes |
| NgsCohortExtract.output_file_base_name | base name for generated output VCF shards|

### Application Notes

- As these exported VCF can be very large, and downstream tools like Hail can imported sharded VCFs, they are left in a sharded state.
