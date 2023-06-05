# Bulk Ingest Instructions


The bulk ingest workflow has one distinct difference from it's non-bulk counterpart and required several more input parameters
It can load in more than 10k samples at a time
We do this by directly mapping the sample ids and their gcp storage path locations and loading them that way.
In order to do this, we need some additional information from you the user.

Since we ask for so many additional inputs, we have tried to streamline the process in two ways.
The first is that we have defaults of the most common values that are input automatically. They are listed below. 
The second is that we use the information we have to the best of our abilities to make educated guesses about what columns you may want to use

The logic for these guesses is also below.

Since you are planning to use a bulk ingest workflow, we recommend that you use a sample set to load your data

## TODO: fellow variant team members, should this be required?!?!? Why would they bulk ingest and NOT use a sample set?

---

#### Input descriptions

The table below describes the GVS bulk ingest variables:

| Input variable name         | Description                                                                                                                                                                | Type |
|-----------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------| --- |
| samples_table_name          | (Optional) The entity that is used for the data table; `sample` is most likely the name of the data table in Terra.                                                        | String |
| sample_id_column_name       | (Optional) TODO: okay this is probably just based on the above value if there's no value below; `sample_id` in the sample data table in Terra.                             | String |
| entity_id_column_name       | (Optional) User defined column id name; `sample_id` is the recommended name--`research_id` is also commonly used.                                                          | String |
| vcf_files_column_name       | (Optional) Column that contains the cloud paths to the sample GVCF files                                                                                                   | String |
| vcf_index_files_column_name | (Optional) Column that contains the cloud paths to the sample GVCF index files                                                                                             | String |
| sample_set_name             | (Optional) TODO: should this be optional!??!      The name of the set of entities                                                                                          | String |
| dataset_name                | Name of the BigQuery dataset used to hold input samples, filtering model data, and other tables created during the workflow.                                               | String |
| bq_project_id               | Name of the Google project that contains the BigQuery dataset.                                                                                                             | String |
| call_set_identifier         | Used to name the filter model, BigQuery extract tables, and final joint VCF shards; should begin with a letter; valid characters include A-z, 0-9, “.”, “,”, “-“, and “_”. | String |

## TODO: do we care about documenting the additional optional input values?

        # File? gatk_override

        # File interval_list = "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list"

        # String drop_state = "NONE"

        # The larger the `load_data_batch_size` the greater the probability of preemptions and non-retryable BigQuery errors,
        # so if specifying `load_data_batch_size`, adjust preemptible and maxretries accordingly. Or just take the defaults, as those should work fine in most cases.
        Int? load_data_batch_size
        Int? load_data_preemptible_override
        Int? load_data_maxretries_override


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

(TODO: as an example, should I choose quickstart instead, since the Beta workspace doesn't even has this as an option?)

The GVS workflow takes in reblocked single sample GVCF files and their corresponding index files as `input_vcfs` and `input_vcf_indexes`, respectively. 

Note:
the samples_table_name value is `samples`
the sample_id_column_name value is `sample_id`
the entity_id_column_name value does not need to be specified as it is `sample_id`
the vcf_files_column_name value is `gvcf`
the vcf_index_files_column_name value is `gvcf_index`
the sample_set_name value is `sample_set_name`




### Parameter estimates and hueristics

step 1:
What is our entity type?
This is how we determine the table name, the entity_set and the entity_id info

A user can easily pass this information as the param "sample_table_name"
The default value for this is "sample"

If we have a user defined entity type, we check that it exists as a table.

If we have a user defined entity_set name, we check that it exists as a table and that the entity that is used has a corresponding id column

step 2:
What columns do you care about?
The ingest needs to know about the GVCF locations in GCP
The GVCF paths and their index files can be discerned based on the limited logic that they will likely have files named appropriately
The workflow looks for paths ending in vcf as well as columns with familiar names such as "vcf" or "reblocked"

### Import sample GVCF files

This step validates that sample GVCF files contain required annotations and loads samples into BigQuery tables. 
Once the samples are loaded into BQ, the next step can be run, which creates a table of Alternate Alleles. This table can be used to create a filter or to enable search

The next recommended wdl to run is:
[GvsPopulateAltAllele](https://github.com/broadinstitute/gatk/blob/ah_var_store/scripts/variantstore/wdl/GvsPopulateAltAllele.wdl) 