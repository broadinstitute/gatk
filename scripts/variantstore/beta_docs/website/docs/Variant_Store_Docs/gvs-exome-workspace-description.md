---
sidebar_position: 5
---
# Genomic Variant Store Exomes workflow (beta) for variant discovery

The [Genomic Variant Store (GVS)](https://github.com/broadinstitute/gatk/blob/ah_var_store/scripts/variantstore/gvs-product-sheet.pdf) was developed as a solution for variant discovery on a large scale. The GVS is powered by BigQuery and creates large joint callsets more reliably with decreased time and cost compared to previous solutions.

This workspace contains a fully reproducible example workflow for Exomes and Blended Genome Exome variant discovery using the GVS workflow. If you have genome data, see the [GVS Genome workspace](https://app.terra.bio/#workspaces/gvs-prod/Genomic_Variant_Store_Beta).

Scroll down for an overview of the workflow, example data, cost estimates, and additional resources.

The materials in this workspace were developed by the Data Sciences Platform at the Broad Institute.

## Workflow Overview

![A diagram depicting the Genomic Variant Store workflow. Sample GVCF files are imported into the core data model. A filtering model is trained using the VETS toolchain, and then applied while the samples are extracted as cohorts in sharded joint VCF files. Each step integrates BigQuery and GATK tools.](https://storage.googleapis.com/terra-featured-workspaces/Genomic-Variant-Store/genomic-variant-store_diagram.png)

### What does it do?

The [GVS workflow](https://github.com/broadinstitute/gatk/blob/ah_var_store/scripts/variantstore/wdl/GvsJointVariantCalling.wdl) is an open-source, cloud-optimized workflow for joint calling at a large scale using the GVS. The workflow takes in single sample GVCF files, loads them into [BigQuery](https://cloud.google.com/bigquery/docs) tables, and combines them into a variant filtering model driven by machine learning. The model is uploaded back into BigQuery and applied to the data. The workflow produces sharded joint VCF files with indices, a manifest file, and metrics.

For workflow documentation, see the [Genomic Variant Store workflow overview](https://github.com/broadinstitute/gatk/blob/ah_var_store/scripts/variantstore/beta_docs/gvs-overview.md).

---

**Important Update:**

Starting September 1 2023, variants in GVS are filtered using the GATK Variant Extract-Train-Score (VETS) toolchain with an isolation-forest model by default. GVS maintains the ability to run VQSR on callsets up to 10,000 genomes to reproduce the same results from past analyses. See the [release notes](https://github.com/broadinstitute/gatk/blob/ah_var_store/scripts/variantstore/docs/release_notes/VETS_Release.pdf) for more information.

---

### What data does it require as input?

- Reblocked single sample GVCF files
- GVCF index files

To see details about input requirements, see [Run Your Own Samples](https://github.com/broadinstitute/gatk/blob/ah_var_store/scripts/variantstore/beta_docs/run-your-own-samples.md). Example GVCF and index files in the Data tab of this workspace are hosted in a public Google bucket and links are provided in the sample data table.

While the GVS has been tested with almost 190,000 single sample exome GVCF files as input, only datasets of up to 100,000 exomes are currently supported by the beta workflow.

Note that, by default, the Exomes Beta workflow uses two interval lists:
1. The calling region interval list defaults to the Blended Genome Exome Interval List: `gs://gcp-public-data--broad-references/hg38/v0/bge_exome_calling_regions.v1.1.interval_list` . For information on how that interval list was generated, see `gs://gcp-public-data--broad-references/hg38/v0/bge_exome_calling_regions.v1.1.interval_list.README.md`.
2. The target region interval list defaults to the [Twist Alliance Clinical Research Exome targets](https://www.twistbioscience.com/resources/safety-data-sheet/twist-alliance-clinical-research-exome-349-mb-bed-files). This interval list is publicly available here: `gs://gcp-public-data--broad-references/hg38/v0/bge_TwistAllianceClinicalResearchExome_Covered_Targets_hg38.interval_list`. 
   1. Any variants called outside these interval will be *genotype* filtered with `OUTSIDE_OF_TARGETS`. See [GVS output documentation](https://github.com/broadinstitute/gatk/blob/ah_var_store/scripts/variantstore/beta_docs/gvs-outputs.md) for more details on this filter.
   2. If your exome data was run using different targets, please obtain the target interval list from your library construction method to utilize this filter.
   3. This filter is a nice to have for standard Exomes. However, it is highly recommended to use this filter for Blended Genome Exome data because in BGE data there are low coverage reads that cause more off target calls. 
   4. If you do not have your target interval list, this is optional and can be omitted from the inputs.

### What does it return as output?

The following files are stored in the Google Cloud Storage path specified in the `extract_output_gcs_dir` workflow input.

- Sharded joint VCF files, index files, the interval lists for each sharded VCF, and a list of the sample names included in the callset.
- A list of the sample names included in the callset called `sample-name-list.txt`
- Size of output VCF files in megabytes
- Manifest file containing the destinations and sizes in bytes of the output sharded joint VCF and index files

There are example outputs from the sample data in the workspace bucket. For more information on GVS outputs, including the number of shards to expect, VCF content, and filters, see the [GVS outputs documentation](https://github.com/broadinstitute/gatk/blob/ah_var_store/scripts/variantstore/beta_docs/gvs-outputs.md).

## Setup

You will need to set up several accounts and projects before you can begin testing the GVS workflow in this workspace. For detailed step-by-step instructions, see the [GVS Beta Quickstart](https://github.com/broadinstitute/gatk/blob/ah_var_store/scripts/variantstore/beta_docs/gvs-quickstart.md) where you will learn how to:

1. Register for Terra
1. Create a Terra billing project
1. Create a GCP project
1. Create a BigQuery dataset
1. Find your Terra proxy group
1. Configure permissions in GCP
1. Clone the GVS beta workspace

For troubleshooting or questions, contact the [Broad Variants team](mailto:variants@broadinstitute.org).

## Running the workflow

The `GvsBeta` workflow in the GVS beta workspace is pre-configured to use 10 sample reblocked GVCF files in the workspace Data tab. To run the example data:

1. **Select the workflow** from the Workflows tab.
1. Select the 'Run workflow with inputs defined by file paths' radio button.
1. Configure the workflow inputs.
    1. Enter a **name for the callset** as a string with the format “*CALLSET_NAME*” for the `call_set_identifier` variable. This string is used as to name several variables and files and should begin with a letter. Valid characters include A-z, 0-9, “.”, “,”, “-“, and “_”.
    1. Enter the name of your **BigQuery dataset** as a string with the format “*DATASET_NAME*” for the `dataset_name` variable. Valid characters include A-z, 0-9, “,”, “-”, and “_”.
    1. Enter the name of the **GCP project** that holds the BigQuery dataset as a string with the format “*PROJECT_NAME*” for the `project_id` variable.
    2. Update the name of the **column of reblocked GVCFs** in your samples table in the workspace Data tab for the `vcf_files_column_name` variable with the format "*COLUMN_NAME*" if it is different from the default. If you used the ReblockGVCF workflow in the workspace without modification, this will be the default string _"reblocked_gvcf"_.
    3. Update the name of the **column of reblocked GVCF index files** in your samples table in the workspace Data tab for the `vcf_index_files_column_name` variable with the format "*COLUMN_NAME*" if it is different from the default. If you used the ReblockGVCF workflow in the workspace without modification, this will be the default string _"reblocked_gvcf_index"_.
    4. Update the name of the **column with your sample IDs** that will be used to identify samples in the callset for the `sample_id_column_name` variable as a string with the format "*COLUMN_NAME*" if it is different from the default. Note that the supplied IDs **MUST** be unique.
    5. Enter the name of the **output gcs bucket** where all outputs listed above will go in the variable `extract_output_gcs_dir` in the format `gs://bucket_name/my_run`. We recommend using the workspace google bucket, which you can find on the Dashboard tab under "Cloud Information">Bucket Name, and making a subdirectory under it for each run.
    1. Ensure that the (optional) `is_wgs` parameter is set to false.
1. **Save** the workflow configuration.
1. **Run** the workflow.

To run the GVS workflow on your own sample data, follow the instructions in the tutorial, [Upload data to Terra and run the GVS workflow](https://github.com/broadinstitute/gatk/blob/ah_var_store/scripts/variantstore/beta_docs/run-your-own-samples.md).

See the "Job History" tab in the Genomic_Variant_Store_Exomes_Beta workspace for a recent example configuration.

### Time and cost
Below are several examples of the time and cost of running the workflow. Generally the cost is around a half a cent per sample, including Terra compute and BigQuery compute cost. This does not include storage costs.

| Number of Samples | Wall Clock Time (hh:mm) | Cost $    | Cost per Sample |
|-------------------|-------------------------|-----------|-----------------|
| 10                | 3:08:00                 | $0.76     | $0.07562        |
| 5000              | 5:21:00                 | $20.41    | $0.00408        |
| 20000             | 10:38:00                | $94.60    | $0.00473        |
| 50000             | 20:29:00                | $250.07   | $0.00500        |
| 100000            | 45:11:00                | $274.15   | $0.00274        |

**Note:** The time and cost listed above represent a single run of the GVS workflow. Actual time and cost may vary depending on BigQuery and Terra load at the time of the callset creation. For example, in practice we've seen 10,000 genome runs take from 13-20 hours.

For more information about controlling Cloud costs, see [this article](https://support.terra.bio/hc/en-us/articles/360029748111).


#### Storage cost

The GVS workflow produces several intermediate files in your BigQuery dataset, and storing these files in the cloud will increase the storage cost associated with your callset. To reduce cloud storage costs, you can delete some of the intermediate files after your callset has been created successfully.

If you plan to create subcohorts of your data, you can delete the tables with `_REF_DATA`, `_SAMPLES`, and `_VET_DATA` at the end of the table name in your BigQuery dataset by following the instructions in the Google Cloud article, [Managing tables](https://cloud.google.com/bigquery/docs/managing-tables#deleting_a_table).

If you don’t plan to create subcohorts of your data, you can delete your BigQuery dataset by following the instructions in the Google Cloud article, [Managing datasets](https://cloud.google.com/bigquery/docs/managing-datasets#deleting_a_dataset). Note that the data will be deleted permanently from this location, but output files can still be found in the workspace bucket.

---

### Additional Resources
* For questions regarding GATK-related tools and Best Practices, see the [GATK website](https://gatk.broadinstitute.org/hc/en-us).
* For Terra-specific documentation and support, see the [Terra Support](https://support.terra.bio/hc/en-us).
* To learn more about the GATK Variant Extract-Train-Score (VETS) toolchain, see the [release notes](https://github.com/broadinstitute/gatk/blob/ah_var_store/scripts/variantstore/docs/release_notes/VETS_Release.pdf).

### Contact Information
* If you have questions or issues while running the GVS workflow in this workspace, contact the [Broad Variants team](mailto:variants@broadinstitute.org).

* You can also Contact the Terra team for questions from the Terra main menu. When submitting a request, it can be helpful to include:
    * Your Project ID
    * Your workspace name
    * Your Bucket ID, Submission ID, and Workflow ID
    * Any relevant log information

### Workspace Citation
If you use plan to publish data analyzed using the GVS workflow, please cite the [GVS Exomes beta workspace](https://app.terra.bio/#workspaces/gvs-prod/Genomic_Variant_Store_Exomes_Beta).

Details on citing Terra workspaces can be found in [How to cite Terra](https://support.terra.bio/hc/en-us/articles/360035343652).

Data Sciences Platform, Broad Institute (*Year, Month Day that this workspace was last modified*) gvs-prod/Genomic_Variant_Store_Exomes_Beta [workspace] Retrieved *Month Day, Year that workspace was retrieved*, https://app.terra.bio/#workspaces/gvs-prod/Genomic_Variant_Store_Exomes_Beta

### License
**Copyright Broad Institute, 2024 | Apache**  
The workflow script is released under the Apache License, Version 2.0 (full license text at https://github.com/broadinstitute/gatk/blob/master/LICENSE.TXT). Note however that the programs called by the scripts may be subject to different licenses. Users are responsible for checking that they are authorized to run all programs before running these tools.

### Workspace Change Log
| Date       | Change                                   | Author             |
|------------|------------------------------------------|--------------------|
| 07/25/2024 | Add new input to support BGE data        | Miguel Covarrubias |
| 06/18/2024 | Add ReblockGVCF and update documentation | Kylee Degatano     |
| 05/08/2024 | Third release of the exomes workspace.   | Bec Asch           |
| 01/09/2024 | Second release of the exomes workspace.  | Bec Asch           |
| 09/08/2023 | First release of the exomes workspace.   | George Grant       |