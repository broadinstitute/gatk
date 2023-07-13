### Genomic Variant Store workflow (beta) for variant discovery

The [Genomic Variant Store (GVS)](https://github.com/broadinstitute/gatk/blob/ah_var_store/scripts/variantstore/gvs-product-sheet.pdf) was developed as a solution for variant discovery on a large scale. The GVS is powered by BigQuery and creates large joint callsets more reliably with decreased time and cost compared to previous solutions.

This workspace contains a fully reproducible example workflow for variant discovery using the GVS workflow.

Scroll down for an overview of the workflow, example data, cost estimates, and additional resources.

The materials in this workspace were developed by the Data Sciences Platform at the Broad Institute.

## Workflow Overview

![A diagram depicting the Genomic Variant Store workflow. Sample GVCF files are imported into the core data model. A filtering model is trained using Variant Quality Score Recalibration, or VQSR, and then applied while the samples are extracted as cohorts in sharded joint VCF files. Each step integrates BigQuery and GATK tools.](https://storage.googleapis.com/terra-featured-workspaces/Genomic-Variant-Store/genomic-variant-store_diagram.png)

### What does it do?

The [GVS workflow](https://github.com/broadinstitute/gatk/blob/ah_var_store/scripts/variantstore/wdl/GvsJointVariantCalling.wdl) is an open-source, cloud-optimized workflow for joint calling at a large scale using the GVS. The workflow takes in single sample GVCF files, loads them into [BigQuery](https://cloud.google.com/bigquery/docs) tables, and combines them into a variant filtering model driven by machine learning. The model is uploaded back into BigQuery and applied to the data. The workflow produces sharded joint VCF files with indices, a manifest file, and metrics.

For workflow documentation, see the [Genomic Variant Store workflow overview](https://github.com/broadinstitute/gatk/blob/ah_var_store/scripts/variantstore/beta_docs/gvs-overview.md).

### What data does it require as input?

- Reblocked single sample GVCF files (`input_vcfs`)
- GVCF index files (`input_vcf_indexes`)

To see details about input requirements, see [Run Your Own Samples](https://github.com/broadinstitute/gatk/blob/ah_var_store/scripts/variantstore/beta_docs/run-your-own-samples.md). Example GVCF and index files in the Data tab of this workspace are hosted in a public Google bucket and links are provided in the sample data table.

While the GVS has been tested with 250,000 single sample genome GVCF files as input, only datasets of up to 25,000 genomes are being used for beta testing.

### What does it return as output?

The following files are stored in the workspace workflow execution bucket under Data>Files, or the google bucket specified in the inputs.

- Sharded joint VCF files, index files, the interval lists for each sharded VCF, and a list of the sample names included in the callset.
- Size of output VCF files in megabytes
- Manifest file containing the destinations and sizes in bytes of the output sharded joint VCF and index files

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

The workflow in the GVS beta workspace is pre-configured to use 10 sample GVCF files in the workspace Data tab.

The workflow is configured to call this input from the data table. To run:

1. **Select the workflow** from the Workflows tab.
1. Configure the workflow inputs.
    1. Enter a **name for the callset** as a string with the format “*CALLSET_NAME*” for the `call_set_identifier` variable. This string is used as to name several variables and files and should begin with a letter. Valid characters include A-z, 0-9, “.”, “,”, “-“, and “_”.
    1. Enter the name of your **BigQuery dataset** as a string with the format “*DATASET_NAME*” for the `dataset_name` variable.
    1. Enter the name of the **GCP project** that holds the BigQuery dataset as a string with the format “*PROJECT_NAME*” for the `project_id` variable.
    2. Enter the name of a **directory for writing outputs**. Enter a string in the format "*gs://your_bucket/here*" in `extract_output_gcs_dir`. If you want the data in your Terra workspace bucket you can find that on the Workspace Dashboard, right panel, under Cloud Information. Copy the "Bucket Name" and use it to the inputs to create a gs path as a string like this "*gs://fc-338fe040-3522-484c-ba48-14b48f9950c2*". If you do not enter a bucket here, the outputs will be in the execution directory under *Files*.
1. **Save** the workflow configuration.
1. **Run** the workflow.

To run the GVS workflow on your own sample data, follow the instructions in the tutorial, [Upload data to Terra and run the GVS workflow](https://github.com/broadinstitute/gatk/blob/ah_var_store/scripts/variantstore/beta_docs/run-your-own-samples.md).

### Important configuration notes

By default, the workflow is set up to write outputs to the workspace Google bucket. If you want to write the outputs to a different cloud storage location, you can specify the cloud path in the `extract_output_gcs_dir` optional input in the workflow configuration.

### Time and cost
Below are several examples of the time and cost of running the workflow. Generally the cost is around 6 cents per sample, including Terra compute and BigQuery compute cost. This does not include storage costs.

| Number of Samples | Wall Clock Time (hh:mm) | Cost $  | Cost per Sample |
|-------------------|-----------------|---------|-----------------|
| 10                | 04:30        | $1.35   |  $0.14          |
| 1000              | 07:24        | $59.64     |  $0.06          |
| 2500              | 08:45        | $141.28     |  $0.06          |
| 5000              | 12:00        | $286.71  |  $0.06          |
| 10000             | 13:41        | $604.97 |  $0.06          |

**Note:** For more information about controlling Cloud costs, see [this article](https://support.terra.bio/hc/en-us/articles/360029748111).

#### Storage cost

The GVS workflow produces several intermediate files in your BigQuery dataset, and storing these files in the cloud will increase the storage cost associated with your callset. To reduce cloud storage costs, you can delete some of the intermediate files after your callset has been created successfully.

If you plan to create subcohorts of your data, you can delete the tables with `_REF_DATA`, `_SAMPLES`, and `_VET_DATA` at the end of the table name in your BigQuery dataset by following the instructions in the Google Cloud article, [Managing tables](https://cloud.google.com/bigquery/docs/managing-tables#deleting_a_table).

If you don’t plan to create subcohorts of your data, you can delete your BigQuery dataset by following the instructions in the Google Cloud article, [Managing datasets](https://cloud.google.com/bigquery/docs/managing-datasets#deleting_a_dataset). Note that the data will be deleted permanently from this location, but output files can still be found in the workspace bucket.

---

### Additional Resources
* For questions regarding GATK-related tools and Best Practices, see the [GATK website](https://gatk.broadinstitute.org/hc/en-us).
* For Terra-specific documentation and support, see the [Terra Support](https://support.terra.bio/hc/en-us).
* To learn more about Variant Quality Score Recalibration (VQSR), see the [GATK tool index](https://gatk.broadinstitute.org/hc/en-us/articles/5257893583259).

### Contact Information
* If you have questions or issues while running the GVS workflow in this workspace, contact the [Broad Variants team](mailto:variants@broadinstitute.org).

* You can also Contact the Terra team for questions from the Terra main menu. When submitting a request, it can be helpful to include:
    * Your Project ID
    * Your workspace name
    * Your Bucket ID, Submission ID, and Workflow ID
    * Any relevant log information

### Workspace Citation
If you use plan to publish data analyzed using the GVS workflow, please cite the [GVS beta workspace](https://app.terra.bio/#workspaces/gvs-prod/Genomic_Variant_Store_Beta).

Details on citing Terra workspaces can be found in [How to cite Terra](https://support.terra.bio/hc/en-us/articles/360035343652).

Data Sciences Platform, Broad Institute (*Year, Month Day that this workspace was last modified*) gvs-prod/Genomic_Variant_Store_Beta [workspace] Retrieved *Month Day, Year that workspace was retrieved*, https://app.terra.bio/#workspaces/gvs-prod/Genomic_Variant_Store_Beta

### License
**Copyright Broad Institute, 2023 | Apache**  
The workflow script is released under the Apache License, Version 2.0 (full license text at https://github.com/broadinstitute/gatk/blob/master/LICENSE.TXT). Note however that the programs called by the scripts may be subject to different licenses. Users are responsible for checking that they are authorized to run all programs before running these tools.

### Workspace Change Log
| Date | Change                                   | Author |
| --- |------------------------------------------| --- |
| 06/29/2023 | Updated to support larger sample sizes.  | Kylee Degatano | 
| 09/08/2022 | Updated information on workflow outputs. | Aaron Hatcher |
| 09/07/2022 | Added licensing information.             | Kaylee Mathews |
| 09/06/2022 | Added note about storage cost.           | Kaylee Mathews |
| 06/24/2022 | First release of the workspace.          | Kaylee Mathews |