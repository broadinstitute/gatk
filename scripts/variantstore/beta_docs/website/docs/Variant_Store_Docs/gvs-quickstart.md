---
sidebar_position: 2
---
# Genomic Variant Store Beta Quickstart

In this Quickstart, you will learn how to use the Genomic Variant Store (GVS) in a [Terra workspace](https://app.terra.bio/#workspaces/gvs-prod/Genomic_Variant_Store_Beta) with provided example data. Using the example data in the [Exome GVS workspace](https://app.terra.bio/#workspaces/gvs-prod/Genomic_Variant_Store_Exomes_Beta) is a similar process.

The [GVS](../gvs-product-sheet.pdf) is a solution for variant discovery on a large scale developed by the Data Sciences Platform at the Broad Institute of MIT and Harvard.

The [GVS beta workspace](https://app.terra.bio/#workspaces/gvs-prod/Genomic_Variant_Store_Beta) contains a fully reproducible example workflow for variant discovery using the GVS workflow.

## Workflow Overview

![Diagram depicting the Genomic Variant Store workflow. Sample GVCF files are imported into the core data model. A filtering model is trained using Variant Extract-Train-Score, or VETS, and then applied while the samples are extracted as cohorts in sharded joint VCF files. Each step integrates BigQuery and GATK tools.](./genomic-variant-store_diagram.png)

The [GVS workflow](https://github.com/broadinstitute/gatk/blob/ah_var_store/scripts/variantstore/wdl/GvsJointVariantCalling.wdl) is an open-source, cloud-optimized workflow for joint calling at a large scale using the Terra platform. The workflow takes in single sample GVCF files with indices and produces sharded joint VCF files with indices, a manifest file, and metrics.

To learn more about the GVS workflow, see the [Genomic Variant Store workflow overview](./gvs-overview.md).

### What data does it require as input?

- Reblocked single sample GVCF files
- GVCF index files

Example GVCF and index files in the Data tab of the [GVS beta workspace](https://app.terra.bio/#workspaces/gvs-prod/Genomic_Variant_Store_Beta) are hosted in a public Google bucket and links are provided in the sample data table.

While the GVS has been tested with 410,000 single sample whole genome GVCF files as input, only datasets of up to 25,000 whole genomes and 100,000 whole exomes are currently supported by the beta workflow.

### What does it return as output?

The following files are stored in the Google Cloud Storage path specified in the `extract_output_gcs_dir` workflow input or in the workspace workflow execution bucket under Data>Files (within the left-hand menu on the "Data" workspace tab , under "Other Data", there is a "Files" link that allows you to navigate the files in the workspace bucket).

- Sharded joint VCF files, index files, the interval lists for each sharded VCF, and a list of the sample names included in the callset.
- A list of the sample names included in the callset called `sample-name-list.txt`
- Size of output VCF files in MB
- Manifest file containing the destinations and sizes in B of the output sharded joint VCF and index files

## Setup

You will need to set up several accounts and projects before you can begin testing the GVS workflow in the [GVS beta workspace](https://app.terra.bio/#workspaces/gvs-prod/Genomic_Variant_Store_Beta). To configure the prerequisites of the workflow, follow the instructions below.

For troubleshooting or questions, contact the [Broad Variants team](mailto:variants@broadinstitute.org).

### Step 1. Register for Terra

The GVS workflow requires Terra to run properly.  If you are new to Terra, you’ll need to [register for a Terra account](https://support.terra.bio/hc/en-us/articles/360028235911).

If you already have a Terra account, you can skip this step.

### Step 2. Create a Terra billing project

While Terra is an open-access platform, it runs on the Google Cloud Platform (GCP), which means you’ll need to pay for storage and analyses you perform in the cloud using a billing project in Terra. For information on the cost of running the GVS workflow in the workspace, see the Time and cost estimates section below.

If you’ve never used GCP before, you are eligible to receive $300 in GCP credits, but unfortunately, you **can not** use these credits to run the GVS pipeline. Create a billing project in Terra by following the step-by-step instructions in the article [Set up billing with $300 Google credits to explore Terra](https://support.terra.bio/hc/en-us/articles/360046295092) but do not accept the invitation to use the free credits.

If you have used GCP before but need to create a billing project in Terra, follow steps 2 and 3 in the article [Set up billing with $300 Google credits to explore Terra](https://support.terra.bio/hc/en-us/articles/360046295092).

If you already have a Terra billing account that you would like to use to test the GVS workflow, you can skip this step.

### Step 3. Create a GCP project

To create a new GCP project for testing the GVS workflow, follow the instructions in the Google Cloud documentation article, [Creating and managing projects](https://cloud.google.com/resource-manager/docs/creating-managing-projects#creating_a_project).

When the project has been created successfully, you will see it appear in your list of projects.

If you already have a GCP project that you would like to use to test the GVS workflow, you can skip this step.

### Step 4. Create a BigQuery dataset

BigQuery datasets store the tables created during the execution of the GVS workflow. A new dataset should be created for each callset you want to create using the workflow. 

Create a dataset in BigQuery inside the GCP project you created in Step 3 (above) by following the instructions in the Google Cloud documentation article, [Creating datasets](https://cloud.google.com/bigquery/docs/datasets#create-dataset). If you already have a google project but are not an owner, you will need at least BigQuery Create Dataset ("bigquery.datasets.create") permissions on the project to do this.

We recommend that you select "Multi-region" for the dataset's "Location type" in order to avoid more restrictive Google quotas.

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

The GVS beta workspaces in Terra are read-only, so you’ll need to clone the workspace to create a copy where you can upload your own data and run the workflow. Clone the workspace using the billing project you created in Step 2 (above) by following the instructions in the article [Make your own project workspace](https://support.terra.bio/hc/en-us/articles/360026130851).

The GVS beta workspace for genomes is [here](https://app.terra.bio/#workspaces/gvs-prod/Genomic_Variant_Store_Beta) and the workspace for exomes is [here](https://app.terra.bio/#workspaces/gvs-prod/Genomic_Variant_Store_Exomes_Beta). 

## Running the workflow

The workflow in the GVS beta workspace is pre-configured to use the 10 sample reblocked GVCF files in the workspace Data tab. See the "Job History" tab in the Genomic_Variant_Store_Beta workspace for a recent example configuration.

To run the GVS workflow on your own sample data, follow the instructions in the tutorial, [Upload data to Terra and run the GVS workflow](https://github.com/broadinstitute/gatk/blob/ah_var_store/scripts/variantstore/beta_docs/run-your-own-samples.md).

### Time and cost
Generally, once over ~100 samples, the GVS costs $0.06 USD in compute per sample to run on genomes and $0.005 USD in compute per sample to run on exomes. Cost to run the example data is below and for more information on cost see the [run your own samples](./run-your-own-samples.md) document. 

Cost and runtime of a callset of the **whole genome** sample data pre-loaded in the genome workspace.

| Number of Genome Samples | Elapsed Time (hh:mm) | Terra Cost $ | BigQuery Cost | Total Cost | Measured Cost per Sample |
|-------------------|----------------------|------------|---------------|------------|-----------------------------|
| 10                | 04:30                | $0.84      | $0.51         | $1.35      | $0.14                       |

Cost and runtime of a callset of the **exome** sample data pre-loaded in the exome workspace.

| Number of Samples | Elapsed Time (hh:mm) |	Cost $ |Cost per Sample | 
|-------------------|----------------------|--------|----------------| 
| 10                | 	03:08               | 	$0.76	 | $0.07562 |

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
**Copyright Broad Institute, 2023 | Apache**
The workflow script is released under the Apache License, Version 2.0 (full license text at https://github.com/broadinstitute/gatk/blob/master/LICENSE.TXT). Note however that the programs called by the scripts may be subject to different licenses. Users are responsible for checking that they are authorized to run all programs before running these tools.
