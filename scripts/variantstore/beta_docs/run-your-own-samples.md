# Upload data to Terra and run the GVS workflow

In this tutorial, you will learn how to upload your own sample GVCF and index files to a Terra workspace and run the Genomic Variant Store (GVS) workflow using that data.

## What is the GVS?

The [GVS](../gvs-product-sheet.pdf) is a solution for variant discovery on a large scale developed by the Data Sciences Platform at the Broad Institute of MIT and Harvard.

The [GVS workflow](https://github.com/broadinstitute/gatk/blob/ah_var_store/scripts/variantstore/wdl/GvsJointVariantCalling.wdl) is an open-source, cloud-optimized workflow for joint calling at a large scale using the Terra platform. The workflow takes in single sample GVCF files with indices and produces sharded joint VCF files with indices, a manifest file, and metrics.

To learn more about the GVS workflow, see the [Genomic Variant Store workflow overview](./gvs-overview.md).

### What does it require as input?

- reblocked single sample GVCF files (`input_vcfs`) with specific annotations described below
- GVCF index files (`input_vcf_indexes`)

While the GVS workflow has been tested with 100,000 single sample GVCF files as input, only datasets of up to 10,000 files are being used for beta testing.


#### Reblocked GVCF files

If your GVCF files have not been reblocked, you can reblock them using the [WARP reblocking workflow](https://github.com/broadinstitute/warp/blob/master/pipelines/broad/dna_seq/germline/joint_genotyping/reblocking/ReblockGVCF.wdl), which is configured in the [ReblockGVCF Terra workspace](https://app.terra.bio/#workspaces/warp-pipelines/ReblockGVCF). 

For more information about reblocking, check out [WARP Whole Genome and Exome Pipelines Produce Reblocked GVCFs](https://broadinstitute.github.io/warp/blog/tags/reblock/).

#### GVCF annotations

Input GVCF files for the GVS workflow must include the annotations described in the table below:

| Annotation | Description | Notes |
| --- | --- | --- |
| Ref | Reference allele. | --- |
| Alt | Alternate allele. | --- |
| AS_RAW_MQ, RAW_MQandDP, or RAW_MQ | RMS mapping quality (‘AS’: allele-specific). | Required for VQSR Data |
| AS_RAW_MQRankSum or Map_QUAL_RANK_SUM_KEY | Z-score from Wilcoxon rank sum test of alternate versus reference read mapping qualities. | Required for VQSR Data |
| QUALapprox | Sum of PL[0] values; used to approximate the QUAL score. | Required for VQSR Data |
| AS_QUALapprox | Allele-specific sum of PL[0] values; used to approximate the QUAL score. | Required for VQSR Data |
| AS_SB_TABLE or STRAND_BIAS_BY_SAMPLE | Allele-specific forward/reverse read counts for strand bias tests. | Required for VQSR Data |
| AS_VarDP, VarDP, or DP | Depth over variant genotypes, or read depth  (‘AS’: allele-specific). | Required for VQSR Data |
| call_GT | Genotype. | --- |
| call_GQ | Genotype quality. | --- |

### What does it return as output?

The following files are stored in the workspace Google bucket and links to the files are written to the sample_set data table:

- sharded joint VCF files and index files
- size of output VCF files in MB
- manifest file containing the output destination of additional files and other metadata

## Setup

Before you can begin uploading your data to Terra, you’ll need to setup some accounts and permissions that will allow Terra to access your data and use BigQuery to run the workflow. Follow the step-by-step instructions in [GVS Beta Quickstart](./gvs-quickstart.md).

## Upload data to your workspace

To run the GVS workflow, your single sample GVCF files need to be stored in the cloud and loaded into a data table in your clone of the [GVS workspace](https://app.terra.bio/#workspaces/gvs-prod/Genomic_Variant_Store_Beta). The procedure is a little different, depending on whether your samples are already stored in the cloud. Follow the step-by-step instructions below to load your sample files into the workspace based on where your files are stored.

### Data stored in the cloud

If your data is already stored in the cloud, you’ll need to upload a TSV file to Terra containing the cloud paths to your files and update permissions on the data to allow Terra to access it. 

1. Navigate to the **Data tab** in your clone of the GVS workspace.

2. Click on the **three vertical dots icon** next to the **sample** data table inside the TABLES sidebar on the left side of the page.

3. Select **Download TSV** to download the data table to your local machine. 

4. **Open the TSV file** with a spreadsheet editor of your choice.

5. **Replace the cloud paths** to the example GVCF and index files in the second and third columns with the cloud paths to your GVCF and index files.

6. **Update the `sample_id` field** for each sample to anything you’d like. These will be used to name the samples in the output joint VCF file.

---

**Warning:**      
The workflow in the GVS beta workspace is configured based on the format of the TSV file. To avoid reconfiguring the workflow, do **not** rearrange or rename the columns in the TSV file.

---

7. Follow steps 2 and 3 in [How to make a data table from scratch or a template](https://support.terra.bio/hc/en-us/articles/6197368140955) to **save and upload the TSV file** to Terra.

8. Grant your Terra proxy group the Storage Object Creator and Storage Object Viewer roles on the Google Cloud Storage (GCS) bucket that holds your sample data by following the **Add a principal to a bucket-level policy** instructions in the Google Cloud documentation article, [Use IAM permissions](https://cloud.google.com/storage/docs/access-control/using-iam-permissions).

### Data NOT stored in the cloud

If your data is not stored in the cloud, you’ll need to upload it to your workspace storage bucket along with a TSV file containing each of the data file names.

1. Navigate to the **Data tab** in your clone of the GVS workspace.

2. Click on the **three vertical dots icon** next to the **sample** data table inside the TABLES sidebar on the left side of the page.

3. Select **Download TSV** to download the data table to your local machine. 

4. **Open the TSV file** with a spreadsheet editor of your choice.

5. **Replace the cloud paths** to the example GVCF and index files in the second and third columns with the names of your GVCF and index files for each sample in your dataset.

6. **Update the `sample_id` field** to anything you’d like for each sample in your dataset. These will be used to name the samples in the output joint VCF file.

---

**Warning:**        
The workflow in the GVS beta workspace is configured based on the format of the TSV file. To avoid reconfiguring the workflow, do **not** rearrange or rename the columns in the TSV file.

---

7. Follow the steps in [How to use the Data Uploader](https://support.terra.bio/hc/en-us/articles/4419428208411) to **upload your data and TSV file** to Terra.

## Run the workflow

Now that your samples are loaded into data table in Terra, it’s time to setup and run the GVS workflow! The workflow is configured to call inputs from and write outputs to the sample_set table.

1. **Select the workflow** from the Workflows tab.
1. In the configuration page, select the **sample_set entity** in Step 1.
1. Create a new sample_set dataset in Step 2.
    1. Select **Create a new sample_set from selected samples**.
    1. Click on the **arrow next to the tickbox** in the header row of the table.
    1. Select All in the dropdown menu that appears to select all of the samples.
    1. If you do **not** want to include the example data in your callset, **click on the tickbox** next to each of the 10 samples to deselect them.
    1. At the bottom of the page, **name your new sample_set**. This will appear in the sample_set data table.
    1. Click **OK** to create the sample_set.
1. Configure the workflow inputs.
    1. Enter a **name for the callset** as a string with the format “*CALLSET_NAME*” for the `call_set_identifier` variable. This string is used as to name several variables and files and should begin with a letter. Valid characters include A-z, 0-9, “.”, “,”, “-“, and “_”.
    1. Enter the name of your **BigQuery dataset** as a string with the format “*DATASET_NAME*” for the `dataset_name` variable.
    1. Enter the name of the **GCP project** that holds the BigQuery dataset as a string with the format “*PROJECT_NAME*” for the `project_id` variable.
1. **Save** the workflow configuration.
1. **Run** the workflow.

### Adding additional data to the BigQuery dataset
The workflow imports GVCFs to the BigQuery dataset. The example GVCFs for the workspace are listed as a set in the workspace sample_set table. Multiple sample sets may be added to the same BigQuery dataset to appear in the same callset, but the workflow only runs on one sample set at a time. 

If you run the workflow with the example data pre-loaded into the workspace, and then load additional data into that same dataset, the example data will be in your callset. Be sure to make a new dataset for your own data if you run the example data (unless you want the example data in there too–the more the merrier!).

### Important configuration notes

By default, the workflow is set up to write outputs to the workspace Google bucket. If you want to write the outputs to a different cloud storage location, you can specify the cloud path in the `extract_output_gcs_dir` optional input in the workflow configuration. 

Accessing the interval lists for all of the output VCF shards:
Start by clicking into your job from the job history tab. This page for your job will have a table with metadata from your job in it (such as workflow ID and the Run Cost) and the table will have a Links column.
The icons/links in that column have additional job information.
* Job Manager
* Workflow Dashboard
* Execution directory


There are two ways to get the interval lists from here. You can either use the UI in the Job Manager, or navigate for the correct path from the Execution directory.

#### Using the Job Manager:
If you prefer to use the job manager, click on the job manager link from the links column. From the List View tab, click into the GvsUnified page. This page's List View tab will show a list of all the sub-workflows in your job. From there, click into the GvsExtractCallset sub-workflow page. This page's List View tab will show a list of all tasks in the extract sub-workflow. You are interested in the SplitIntervals task which should be the 8th task. In the SplitIntervals row, find the icon in the outputs column and click on it. A popover will load, which will contain an array of the interval_files. These are the paths for the interval lists.
If the popover does not load, you can click into the SplitIntervals Execution directory (which is in the Links column of this table).

#### Using the Execution directory:
Your jobs Execution directory contains artifacts from your job, including logs and the interval lists.
You can access the interval lists by navigating to the directory for the SplitIntervals task inside the GvsExtractCallset sub-workflow.
Once you have clicked on the link to the Execution directory, it will drop you into your workflow bucket.
From there, click into the call-GvsUnified --> GvsUnified --> <unified job id (this will be the only option, so no lookup is needed)>
Then select the extract sub-workflow: call-GvsExtractCallset --> GvsExtractCallset --> <extract job id (only option--no lookup needed)>
Finally, select the SplitIntervals task: call-SplitIntervals and find the glob-<> directory. Inside that directory are the interval lists for the call set.

The interval list paths will look something like this:
`gs://<workspace bucket id>/<submission id>/GvsJointVariantCalling/<workflow id>/call-GvsUnified/GvsUnified/<unified job id (only option--no lookup needed)>/call-GvsExtractCallset/GvsExtractCallset/<extract job id (only option--no lookup needed)>/call-SplitIntervals/glob-<task id>/0000000000-<callset id>.vcf.gz.interval_list`

There will also be a `glob-<task id>` file with a list of interval lists, but not their paths.

Note:
The interval lists are named consistently with the vcfs: 00000000.vcf.gz.interval-list will go with 00000000.vcf.gz and 00000000.vcf.gz.tbi 

### Time and cost estimates
Below are several examples of the time and cost of running the workflow.

| Number of Samples | Wall Clock Time | Cost $  | Cost per Sample |
|-------------------|-----------------|---------|-----------------|
| 10                | 04:30:00        | $0.84   | ~$0.08          |
| 1000              | 07:24:00        | $$$     | ~$0.07          |
| 2500              | 08:45:00        | $$$     | ~$0.07          |
| 5000              | 12:00:00        | $$$     | ~$0.07          |
| 10000             | 00:00:00        | $683.33 | ~$0.07          |

**Note:** For more information about controlling Cloud costs, see [this article](https://support.terra.bio/hc/en-us/articles/360029748111).

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

### Citing the GVS workflow
If you use plan to publish data analyzed using the GVS workflow, please cite the [GVS beta workspace](https://app.terra.bio/#workspaces/gvs-prod/Genomic_Variant_Store_Beta).

Details on citing Terra workspaces can be found here: [How to cite Terra](https://support.terra.bio/hc/en-us/articles/360035343652)

Data Sciences Platform, Broad Institute (*Year, Month Day that the workspace was last modified*) gvs-prod/Genomic_Variant_Store_Beta [workspace] Retrieved *Month Day, Year that workspace was retrieved*, https://app.terra.bio/#workspaces/gvs-prod/Genomic_Variant_Store_Beta

### License
**Copyright Broad Institute, 2020 | BSD-3**  
All code provided in this workspace is released under the WDL open source code license (BSD-3) (full license text at https://github.com/broadinstitute/warp/blob/develop/LICENSE). Note however that the programs called by the scripts may be subject to different licenses. Users are responsible for checking that they are authorized to run all programs before running these tools.