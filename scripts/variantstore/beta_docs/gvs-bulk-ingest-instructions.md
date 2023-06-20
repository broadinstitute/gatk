# Bulk Ingest Instructions


The bulk ingest workflow has one distinct difference from its non-bulk counterpart and requires several more input parameters.
It can load in more than 10k samples at a time by directly mapping sample ids to their GCS storage path locations.
In order to do this, the bulk ingest workflow needs to collect some additional information from the user.

Since we ask for so many additional inputs, we have tried to streamline the process in two ways.
The first is that we have defaults of the most common values that are input automatically. These are listed below. 
The second is that we use the information we have to the best of our abilities to make educated guesses about what columns you may want to use.

The logic for these guesses is also below.

Since you are planning to use a bulk ingest workflow, we recommend that you use a sample set to load your data.

---

#### Input descriptions

The table below describes the GVS bulk ingest variables:

| Input variable name      | Description                                                                                                                                                                     | Type |
|--------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------| --- |
| data_table_name          | (Optional) The name of the data table; This table hold the GVCFs to be ingested; `sample` is most likely the name of the data table in Terra, so `sample` is set as the default | String |
| entity_id_column_name    | (Optional) User defined column id name; `sample_id` is the recommended name;`research_id` is also commonly used.                                                                | String |                                                         | String |
| vcf_files_column_name    | (Optional) Column that contains the GCS paths to the sample GVCF files.                                                                                                         | String |
| vcf_index_files_column_name | (Optional) Column that contains the GCS paths to the sample GVCF index files.                                                                                                   | String |
| sample_set_name          | (Optional) The name of the set of samples / entities.                                                                                                                           | String |
| dataset_name             | Name of the BigQuery dataset used to hold input samples, filtering model data, and other tables created during the workflow.                                                    | String |
| project_id               | Name of the Google project that contains the BigQuery dataset.                                                                                                                  | String |
| call_set_identifier      | Used to name the filter model, BigQuery extract tables, and final joint VCF shards. Should begin with a letter, valid characters include A-z, 0-9, “.”, “,”, “-“, and “_”.      | String |


The first four parameters are labelled as optional, but the pipeline needs a value for each of them to process the data. This is because the system can likely determine the four values it needs.
If you know the values that you want for these 4, feel free to fill them in, but it is not necessary.
If the user does not supply all 4 values, the system uses a mix of default values and heuristics to determine these parameters.
The first parameter, `data_table_name`, is the name of the data entity and is usually "sample," which is the default value for the input.
The second parameter, `entity_id_column_name`, corresponds to the `sample_name` value in the GVS database. Some advanced users prefer to use a custom column. This defaults to the value of `data_table_name` + "_id" which is most likely "sample_id".

The next two parameters, `vcf_files_column_name` and `vcf_index_files_column_name`, name the two columns that track the GCP locations of the GVCFs and their corresponding index files.
These two columns have no explicit defaults. If they are not input by the user, the workflow uses heuristics to determine what they might be
The workflow does this by looking at all the possible columns and a subset of the data in each to predict which columns have the sample GVCF and index file paths.

The workflow looks for paths ending in `vcf` or similar.
The workflow will also validate that any user input for these parameters fits with the expected relationships between them.

The fifth parameter (`sample_set_name`) is required if a sample set is being used for ingest, which we do recommend for large batches. 
You can create a sample set from the data tab of your workspace, by selecting the samples that you desire to include, clicking edit and selecting "Save selection as set". 
You will then be given the option to name your sample set. That name is what you will use as the `sample_set_name` parameter.
We validate that the requested sample_set exists when it is input as a parameter. There is a `<data_table_name>_id` column that is used by the sample_set, but a different sample_name column can still be specified independently to differentiate the samples.




## Setup

### Workflow requirements

#### Input GVCF files

The GVS workflow takes in reblocked GVCF files as input. If your files are not already reblocked, you can reblock them using the [WARP reblocking workflow](https://github.com/broadinstitute/warp/blob/master/pipelines/broad/dna_seq/germline/joint_genotyping/reblocking/ReblockGVCF.wdl), which is configured in the [ReblockGVCF Terra workspace](https://app.terra.bio/#workspaces/warp-pipelines/ReblockGVCF). For more information about reblocking, check out [WARP Whole Genome and Exome Pipelines Produce Reblocked GVCFs](https://broadinstitute.github.io/warp/blog/tags/reblock/).

The workflow also requires specific annotations in input GVCF files, which are described in the tutorial, [Upload data to Terra and run the GVS workflow](./run-your-own-samples.md).

### Inputs

The GVS workflow inputs are described in the sections below and are specified in the Terra GVS beta workspace [workflow configuration](https://app.terra.bio/#workspaces/gvs-prod/Genomic_Variant_Store_Beta/workflows/help-terra/GvsJointVariantCalling).

The workflow is configured to use hg38 (aka GRCh38) as the reference genome.

#### Sample data

The [GVS beta workspace](https://app.terra.bio/#workspaces/gvs-prod/Genomic_Variant_Store_Beta) has been configured with example reblocked GVCF files that you can use to test the workflow.

The GVS workflow takes in reblocked single sample GVCF files and their corresponding index files as `input_vcfs` and `input_vcf_indexes`, respectively. 

In the Beta user workflow, the value examples are as follows:
the samples_table_name value is `sample`
the entity_id_column_name value is `sample_id`
the vcf_files_column_name value is `gvcf`
the vcf_index_files_column_name value is `gvcf_index`
the sample_set_name value is `sample_set_name`