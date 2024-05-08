# Upload data to Terra and run the GVS workflow

In this tutorial, you will learn how to upload your own sample GVCF and index files to a Terra workspace and run the Genomic Variant Store (GVS) workflow using that data.

## What is the GVS?

The [GVS](../gvs-product-sheet.pdf) is a solution for variant discovery on a large scale developed by the Data Sciences Platform at the Broad Institute of MIT and Harvard.

The [GVS workflow](https://github.com/broadinstitute/gatk/blob/ah_var_store/scripts/variantstore/wdl/GvsJointVariantCalling.wdl) is an open-source, cloud-optimized workflow for joint calling at a large scale using the Terra platform. The workflow takes in single sample GVCF files with indices and produces sharded joint VCF files with indices, a manifest file, and metrics.

To learn more about the GVS workflow, see the [Genomic Variant Store workflow overview](./gvs-overview.md).

### What does it require as input?

- Reblocked single sample GVCF files with specific annotations described below
- GVCF index files 

While the GVS has been tested with 250,000 single sample whole genome GVCF files as input, only datasets of up to 25,000 whole genomes are being used for beta testing.

#### Reblocked GVCF files

If your GVCF files have not been reblocked, you can reblock them using the [WARP reblocking workflow](https://github.com/broadinstitute/warp/blob/master/pipelines/broad/dna_seq/germline/joint_genotyping/reblocking/ReblockGVCF.wdl), which is configured in the [ReblockGVCF Terra workspace](https://app.terra.bio/#workspaces/warp-pipelines/ReblockGVCF). 

For more information about reblocking, check out [WARP Whole Genome and Exome Pipelines Produce Reblocked GVCFs](https://broadinstitute.github.io/warp/blog/tags/reblock/).

#### GVCF annotations

Input GVCF files for the GVS workflow must include the annotations described in the table below:

| Annotation                                | Description                                                                               | Notes                      |
|-------------------------------------------|-------------------------------------------------------------------------------------------|----------------------------|
| Ref                                       | Reference allele.                                                                         | ---                        |
| Alt                                       | Alternate allele.                                                                         | ---                        |
| AS_RAW_MQ, RAW_MQandDP, or RAW_MQ         | RMS mapping quality (‘AS’: allele-specific).                                              | Required for VETS training |
| AS_RAW_MQRankSum or Map_QUAL_RANK_SUM_KEY | Z-score from Wilcoxon rank sum test of alternate versus reference read mapping qualities. | Required for VETS training |
| QUALapprox                                | Sum of PL[0] values; used to approximate the QUAL score.                                  | Required for VETS training |
| AS_QUALapprox                             | Allele-specific sum of PL[0] values; used to approximate the QUAL score.                  | Required for VETS training |
| AS_SB_TABLE or STRAND_BIAS_BY_SAMPLE      | Allele-specific forward/reverse read counts for strand bias tests.                        | Required for VETS training |
| AS_VarDP, VarDP, or DP                    | Depth over variant genotypes, or read depth  (‘AS’: allele-specific).                     | Required for VETS training |
| call_GT                                   | Genotype.                                                                                 | ---                        |
| call_GQ                                   | Genotype quality.                                                                         | ---                        |

### What does it return as output?

The following files are stored in the Google Cloud Storage path specified in the `extract_output_gcs_dir` workflow input.

- Sharded joint VCF files, index files, the interval lists for each sharded VCF, and a list of the sample names included in the callset.
- Size of output VCF files in MB
- Manifest file containing the output destination of additional files and other metadata

Note:
The interval lists are named consistently with the vcfs: 00000000.vcf.gz.interval-list will go with 00000000.vcf.gz and 00000000.vcf.gz.tbi

## Setup

Before you can begin uploading your data to Terra, you’ll need to setup some accounts and permissions that will allow Terra to access your data and use BigQuery to run the workflow. Follow the step-by-step instructions in [GVS Beta Quickstart](./gvs-quickstart.md).

## Upload data to your workspace

To run the GVS workflow, your single sample GVCF files need to be stored in the cloud and loaded into a data table in your clone of the [GVS workspace](https://app.terra.bio/#workspaces/gvs-prod/Genomic_Variant_Store_Beta). The procedure is a little different, depending on whether your samples are already stored in the cloud. Follow the step-by-step instructions below to load your sample files into the workspace based on where your files are stored.

### Download the sample table template

The GVS workflow relies on the table structure and names in the Beta workspace. Download a TSV of the `sample` table to use as a template for uploading your own data.

1. Navigate to the **Data tab** in your clone of the GVS workspace.
2. Click on the **three vertical dots icon** next to the **sample** data table inside the TABLES sidebar on the left side of the page.
3. Select **Download TSV** to download the data table to your local machine.

### Delete the example data

You will need to delete the example sample data from your workspace so that it is not included in your callset.

1. Navigate to the **Data tab** in your clone of the GVS workspace.
2. Click on the **sample** data table.
3. Click the **top check mark** to select all samples.
4. Click **Edit** and click **Delete selected rows** to delete the data.

### Edit TSV to upload data stored in the cloud

If your data is already stored in the cloud, you’ll need to upload a TSV file to Terra containing the cloud paths to your files and update the cloud permissions on the data to allow your Terra proxy group to access it.

1. **Open the TSV file** with a spreadsheet editor of your choice.
2. **Replace the cloud paths** to the example GVCF and index files in the second and third columns with the cloud paths to your GVCF and index files.
3. **Update the `sample_id` field** for each sample to anything you’d like. These will be used to name the samples in the output joint VCF files.

---

**Warning:**
The workflow in the GVS beta workspace is configured based on the format of the TSV file. To avoid reconfiguring the workflow, do **not** rearrange or rename the columns in the TSV file. If you need to customize the column names, read more about the parameters to use in [GVS Bulk Ingest Details](https://github.com/broadinstitute/gatk/blob/ah_var_store/scripts/variantstore/docs/gvs-bulk-ingest-details.md).

---

4. Follow steps 2 and 3 in [How to make a data table from scratch or a template](https://support.terra.bio/hc/en-us/articles/6197368140955) to **save and upload the TSV file** to Terra.
5. Grant your Terra proxy group the Storage Object Creator and Storage Object Viewer roles on the Google Cloud Storage (GCS) bucket that holds your sample data by following the **Add a principal to a bucket-level policy** instructions in the Google Cloud documentation article, [Use IAM permissions](https://cloud.google.com/storage/docs/access-control/using-iam-permissions).

### Edit TSV to upload data NOT stored in the cloud

If your data is not stored in the cloud, you’ll need to upload it to your workspace storage bucket along with a TSV file containing each of the data file names.

1. **Open the TSV file** with a spreadsheet editor of your choice.
2. **Replace the cloud paths** to the example GVCF and index files in the second and third columns with the names of your GVCF and index files for each sample in your dataset.
3. **Update the `sample_id` field** to anything you’d like for each sample in your dataset. These will be used to name the samples in the output joint VCF file.

---

**Warning:**        
The workflow in the GVS beta workspace is configured based on the format of the TSV file. To avoid reconfiguring the workflow, do **not** rearrange or rename the columns in the TSV file. If you need to customize the column names, read more about the parameters to use in [GVS Bulk Ingest Details](https://github.com/broadinstitute/gatk/blob/ah_var_store/scripts/variantstore/docs/gvs-bulk-ingest-details.md).

---
4. Follow the steps in [How to use the Data Uploader](https://support.terra.bio/hc/en-us/articles/4419428208411) to **upload your data and TSV file** to Terra.

## Run the workflow

Now that your samples are loaded into data table in Terra, it’s time to setup and run the GVS workflow!

1. **Select the workflow** from the Workflows tab.
1. Configure the workflow inputs.
    1. Enter a **name for the callset** as a string with the format “*CALLSET_NAME*” for the `call_set_identifier` variable. This string is used as to name several variables and files and should begin with a letter. Valid characters include A-z, 0-9, “.”, “,”, “-“, and “_”.
    1. Enter the name of your **BigQuery dataset** as a string with the format “*DATASET_NAME*” for the `dataset_name` variable. Valid characters include A-z, 0-9, “,”, “-”, and “_”.
    1. Enter a GCS path for `extract_output_gcs_dir`, which is where the callset VCFs, VCF indexes, interval lists, and manifest will be copied.
    1. Enter the name of the **GCP project** that holds the BigQuery dataset as a string with the format “*PROJECT_NAME*” for the `project_id` variable.
    1. Ensure that the (optional) `is_wgs` parameter is set to false.
    1. Check the workspace Data tab, and fill out the below inputs according to the data column names that contain the values you want to use:
        - `sample_id_column_name`: a unique identifier for each sample (the workflow will fail if there are any duplicates)
        - `vcf_files_column_name`: the path to the sample VCF file
        - `vcf_index_files_column_name`: the path to the sample VCF index file
1. **Save** the workflow configuration.
1. **Run** the workflow.

### Time and cost
Below are several examples of the time and cost of running the workflow.

| Number of Samples | Elapsed Time (hh:mm) | Terra Cost | BigQuery Cost | Total Cost | Approximate Cost per Sample |
|-------------------|----------------------|------------|---------------|------------|-----------------------------|
| 10                | 04:30                | $0.84      | $0.51         | $1.35      | $0.14                       |
| 1000              | 07:24                | $13.02     | $46.62        | $59.64     | $0.06                       |
| 2500              | 08:45                | $25.10     | $116.18       | $141.28    | $0.06                       |
| 5000              | 12:00                | $54.00     | $232.71       | $286.71    | $0.06                       |
| 10000             | 13:41                | $138.1     | $466.87       | $604.97    | $0.06                       |


**Note:** The time and cost listed above each represent a single run of the GVS workflow. Actual time and cost may vary depending on BigQuery and Terra load at the time of the callset creation.

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

### Citing the GVS workflow
If you use plan to publish data analyzed using the GVS workflow, please cite the [GVS beta workspace](https://app.terra.bio/#workspaces/gvs-prod/Genomic_Variant_Store_Beta).

Details on citing Terra workspaces can be found here: [How to cite Terra](https://support.terra.bio/hc/en-us/articles/360035343652)

Data Sciences Platform, Broad Institute (*Year, Month Day that the workspace was last modified*) gvs-prod/Genomic_Variant_Store_Beta [workspace] Retrieved *Month Day, Year that workspace was retrieved*, https://app.terra.bio/#workspaces/gvs-prod/Genomic_Variant_Store_Beta

### License
**Copyright Broad Institute, 2023 | Apache**
The workflow script is released under the Apache License, Version 2.0 (full license text at https://github.com/broadinstitute/gatk/blob/master/LICENSE.TXT). Note however that the programs called by the scripts may be subject to different licenses. Users are responsible for checking that they are authorized to run all programs before running these tools.
