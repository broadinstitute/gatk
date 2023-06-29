# Bulk Ingest Details


To enable the ingestion of >10k samples at a time, we now support bulk ingestion. It is most straightforward to use the default instructions in [Run Your Own Samples](https://github.com/broadinstitute/gatk/blob/ah_var_store/scripts/variantstore/beta_docs/run-your-own-samples.md). However, if your data doesn't follow the recommended table structure, this document outlines the optional parameters you can use to use GVS.

---

#### Input descriptions

The table below describes the GVS ingest variables:

| Input variable name         | Description                                                                                                                                                                | Type |
|-----------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------| --- |
| samples_table_name          | (Optional) The name of the data table; This table holds the GVCFs to be ingested; `sample` is the default.                                                                 | String |
| sample_id_column_name       | (Optional) User defined sample id column name; `sample_id` is the recommended name;                                                                                        | String |
| vcf_files_column_name       | (Optional) Column that contains the GCS paths to the sample GVCF files. Determined by file inspection if not supplied.                                                     | String |
| vcf_index_files_column_name | (Optional) Column that contains the GCS paths to the sample GVCF index files. Determined by file inspection if not supplied.                                               | String |
| sample_set_name             | (Optional) The name of the set of samples. Only required when doing subsets of data in the `sample` table.                                                                 | String |
| dataset_name                | Name of the BigQuery dataset used to hold input samples, filtering model data, and other tables created during the workflow.                                               | String |
| project_id                  | Name of the Google project that contains the BigQuery dataset.                                                                                                             | String |
| call_set_identifier         | Used to name the filter model, BigQuery extract tables, and final joint VCF shards. Should begin with a letter, valid characters include A-z, 0-9, “.”, “,”, “-“, and “_”. | String |


If the user does not supply all 4 optional inputs, GVS uses a mix of default values and heuristics to determine these parameters.
The first parameter, `samples_table_name`, is the name of the main data table and is usually "sample," which is the default value for the input.
The second parameter, `sample_id_column_name`, corresponds to the `sample_name` value in the GVS database. Some advanced users prefer to use a custom column. This defaults to the value of `data_table_name` + "_id" which is most likely "sample_id".

The next two parameters, `vcf_files_column_name` and `vcf_index_files_column_name`, name the two columns that track the GCP locations of the GVCFs and their corresponding index files.
These two columns have no explicit defaults. If they are not input by the user, the workflow uses heuristics to determine what they might be
The workflow does this by looking at all the possible columns and a subset of the data in each to predict which columns have the sample GVCF and index file paths. The workflow looks for paths ending in `vcf` or similar.
The workflow will also validate that any user input for these parameters fits with the expected relationships between them.

The fifth parameter (`sample_set_name`) is required if a sample set is being used for ingest. We only recommend using a sample set when joint calling a subset of the data in the sample table.

##### Running on a sample set
1. On the data tab of your workspace, select the samples that you want to include.
2. Click edit and select "Save selection as set". 
3. Name your sample set. The name you choose is what you will use as the `sample_set_name` parameter.
4. Return to the input configuration and put the sample set name into the input `sample_set_name`. We validate that the requested sample_set exists when it is input as a parameter.

#### Sample data

The [GVS beta workspace](https://app.terra.bio/#workspaces/gvs-prod/Genomic_Variant_Store_Beta) has been configured with example reblocked GVCF files that you can use to test the workflow.

The GVS workflow takes in reblocked single sample GVCF files and their corresponding index files as `input_vcfs` and `input_vcf_indexes`, respectively. 

In the Beta user workflow, the value examples are as follows:
the samples_table_name value is `sample`
the entity_id_column_name value is `sample_id`
the vcf_files_column_name value is `gvcf`
the vcf_index_files_column_name value is `gvcf_index`
the sample_set_name value is `sample_set`