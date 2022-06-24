# Genomic Variant Store Beta Quickstart

In this Quickstart, you will learn how to use the Genomic Variant Store (GVS) in a [Terra workspace](https://app.terra.bio/#workspaces/gvs-prod/Genomic_Variant_Store_Beta) with provided example data.

The [GVS](../gvs-product-sheet.pdf) is a solution for variant discovery on a large scale developed by the Data Sciences Platform at the Broad Institute of MIT and Harvard.

The [GVS beta workspace](https://app.terra.bio/#workspaces/gvs-prod/Genomic_Variant_Store_Beta) contains a fully reproducible example workflow for variant discovery using the GVS workflow.

## Workflow Overview

![Diagram depicting the Genomic Variant Store workflow. Sample GVCF files are imported into the core data model. A filtering model is trained using Variant Quality Score Recalibration, or VQSR, and then applied while the samples are extracted as cohorts in sharded joint VCF files. Each step integrates BigQuery and GATK tools.](/scripts/variantstore/beta_docs/genomic-variant-store_diagram.png)

The [GVS workflow](https://github.com/broadinstitute/gatk/blob/rc-vs-483-beta-user-wdl/scripts/variantstore/wdl/GvsJointVariantCalling.wdl) is an open-source, cloud-optimized workflow for joint calling at using the GVS. The workflow takes in single sample GVCF files with indices and produces sharded joint VCF files with indices, a manifest file, and metrics.

To learn more about the GVS workflow, see the [Genomic Variant Store workflow overview](./gvs-overview.md).

### What data does it require as input?

- unblocked single sample GVCF files (`input_vcfs`)
- GVCF index files (`input_vcf_indexes`)

Example GVCF and index files in the Data tab of the [GVS beta workspace](https://app.terra.bio/#workspaces/gvs-prod/Genomic_Variant_Store_Beta) are hosted in a public Google bucket and links are provided in the sample data table.

While the GVS workflow has been tested with 100,000 single sample GVCF files as input, only datasets of up to 10,000 files are being used for beta testing.

### What does it return as output?

The following files are stored in the workspace Google bucket and links to the files are written to the `sample_set` data table:

- sharded joint VCF files and index files
- size of output VCF files in MB
- manifest file containing the destinations and sizes in B of the output sharded joint VCF and index files

## Setup

You will need to set up several accounts and projects before you can begin testing the GVS workflow in the [GVS beta workspace](https://app.terra.bio/#workspaces/gvs-prod/Genomic_Variant_Store_Beta). To configure the prerequisites of the workflow, follow the instructions below.

For troubleshooting or questions, contact the [Broad Variants team](mailto:variants@broadinstitute.org).

### Step 1. Register for Terra

The GVS workflow requires a Terra data table to run properly.  If you are new to Terra, you’ll need to [register for a Terra account](https://support.terra.bio/hc/en-us/articles/360028235911).

If you already have a Terra account, you can skip this step.

### Step 2. Create a Terra billing project

While Terra is an open-access platform, it runs on the Google Cloud Platform (GCP), which means you’ll need to pay for storage and analyses you perform in the cloud using a billing project in Terra. For information on the cost of running the GVS workflow in the workspace, see the Time and cost estimates section below.

If you’ve never used GCP before, you are eligible to receive $300 in GCP credits. You can use these credits to create a billing project in Terra by following the step-by-step instructions in the article [Set up billing with $300 Google credits to explore Terra](https://support.terra.bio/hc/en-us/articles/360046295092).

If you have used GCP before but need to create a billing project in Terra, follow steps 2 and 3 in the article [Set up billing with $300 Google credits to explore Terra](https://support.terra.bio/hc/en-us/articles/360046295092).

If you already have a Terra billing account that you would like to use to test the GVS workflow, you can skip this step.

### Step 3. Create a GCP project

To create a new GCP project for testing the GVS workflow, follow the instructions in the Google Cloud documentation article, [Creating and managing projects](https://cloud.google.com/resource-manager/docs/creating-managing-projects#creating_a_project).

When the project has been created successfully, you will see it appear in your list of projects.

If you already have a GCP project that you would like to use to test the GVS workflow, you can skip this step.

### Step 4. Create a BigQuery dataset

BigQuery datasets store the tables created during the execution of the GVS workflow. A new dataset should be created for each callset you want to create using the workflow. Samples can be added to BiqQuery datasets cumulatively and any samples in the dataset will be included in the extracted callset. 

Create a dataset in BigQuery inside the GCP project you created in Step 3 (above) by following the instructions in the Google Cloud documentation article, [Creating datasets](https://cloud.google.com/bigquery/docs/datasets#create-dataset).

Click on the arrow next to the name of your GCP project. When the BigQuery dataset has been created successfully, it will appear under the name of the GCP project.

### Step 5. Find your Terra proxy group

Proxy groups permit Terra to interface with GCP on your behalf, allowing you to run workflows on data in the cloud. Follow these steps to find your Terra proxy group.

1. Select the main menu (three horizontal line icon) at the top left of any Terra page. 
1. Click on your account name to open a dropdown menu.
1. Select Profile to open your profile in Terra. 
1. Find the Proxy Group field, which lists your proxy group underneath. 

For more information on proxy groups in Terra, see [Pet service accounts and proxy groups](https://support.terra.bio/hc/en-us/articles/360031023592). 

You will need your proxy group to grant Terra access to your GCP project and BigQuery dataset in the next step.

### Step 6. Configure permissions in GCP

Your Terra proxy group needs to be granted specific roles in GCP so Terra can create, read, and write to BigQuery tables as part of the GVS workflow.

Grant your proxy group the BigQuery Data Editor, BigQuery Job User, and BigQuery Read Session User roles on the Google Project that you created in Step 3 (above) that contains your BigQuery dataset. For step-by-step instructions, see the Google Cloud documentation article, [Manage access to projects, folders, and organizations](https://cloud.google.com/iam/docs/granting-changing-revoking-access#grant-single-role).

If you’ve done this correctly, you should see your Terra proxy group listed in the table of principals along with the roles you granted it.

### Step 7. Clone the GVS beta workspace

The GVS beta workspace in Terra is read-only, so you’ll need to clone the workspace to create a copy where you can upload your own data and run the workflow. Clone the workspace using the billing project you created in Step 2 (above) by following the instructions in the article [Make your own project workspace](https://support.terra.bio/hc/en-us/articles/360026130851).

## Running the workflow

The workflow in the GVS beta workspace is pre-configured to use 10 sample GVCF files as a sample set in the sample_set data table in the workspace Data tab.

The workflow is configured to call this input from the data table. To run:

1. Select the workflow from the Workflows tab.
1. In the configuration page, select the sample_set entity in Step 1.
1. Select the sample_set dataset in Step 2.
1. Configure the workflow inputs.
    1. Enter a name for the callset as a string with the format “*CALLSET_NAME*” for the `callset_identifier` variable. This string is used as to name several variables and files and should begin with a letter. Valid characters include A-z, 0-9, “.”, “,”, “-“, and “_”.
    1. Enter the name of the BigQuery dataset you created in Step 4 of Setup as a string with the format “*DATASET_NAME*” for the `dataset_name` variable.
    1. Enter the name of the GCP project you created in Step 3 of Setup as a string with the format “*PROJECT_NAME*” for the `project_id` variable.
1. Save the workflow configuration.
1. Run the workflow.

To run the GVS workflow on your own sample data, follow the instructions in the tutorial, [Upload data to Terra and run the GVS workflow](./run-your-own-samples.md).
   
### Important configuration notes

By default, the workflow is set up to write outputs to the workspace Google bucket. If you want to write the outputs to a different cloud storage location, you can specify the cloud path in the `extract_output_gcs_dir` optional input in the workflow configuration. 

### Time and cost estimates
Below is an example of the time and cost of running the workflow.

| Number of Samples | Wall Clock Time | Cost $ | Cost per Sample |
| ---  | --- | --- | --- |
| 10 | 04:30:00 | $0.84 | ~$0.08 |
| 1528 |
| 10379 | 05:38:00 | $683.33 | $0.0658 |

**Note:** For more information about controlling Cloud costs, see [this article](https://support.terra.bio/hc/en-us/articles/360029748111).

---

### Additional Resources
* For questions regarding GATK-related tools and Best Practices, see the [GATK website](https://gatk.broadinstitute.org/hc/en-us).
* For Terra-specific documentation and support, see the [Terra Support](https://support.terra.bio/hc/en-us).
* To learn more about Variant Quality Score Recalibration (VQSR), see the [GATK tool index](https://gatk.broadinstitute.org/hc/en-us/articles/5257893583259).

### Contact Information
* If you have questions or issues while running the GVS workflow in the GVS beta workspace, contact the [Broad Variants team](mailto:variants@broadinstitute.org).

* You can also Contact the Terra team for questions from the Terra main menu. When submitting a request, it can be helpful to include:
    * Your Project ID
    * Your workspace name
    * Your Bucket ID, Submission ID, and Workflow ID
    * Any relevant log information

### Citing the GVS workflow
If you use plan to publish data analyzed using the GVS workflow, please cite the [GVS beta workspace](https://app.terra.bio/#workspaces/gvs-prod/Genomic_Variant_Store_Beta).

Details on citing Terra workspaces can be found here: [How to cite Terra](https://support.terra.bio/hc/en-us/articles/360035343652)

Data Sciences Platform, Broad Institute (*Year, Month Day that the workspace was last modified*) gvs-prod/Genomic_Variant_Store_Beta [workspace] Retrieved *Month Day, Year that workspace was retrieved*, https://app.terra.bio/#workspaces/gvs-prod/Genomic_Variant_Store_Beta

### License
**Copyright Broad Institute, 2020 | BSD-3**  
All code provided in the workspace is released under the WDL open source code license (BSD-3) (full license text at https://github.com/broadinstitute/warp/blob/develop/LICENSE). Note however that the programs called by the scripts may be subject to different licenses. Users are responsible for checking that they are authorized to run all programs before running these tools.