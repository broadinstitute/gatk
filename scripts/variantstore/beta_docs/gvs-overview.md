# Genomic Variant Store pipeline overview
 
| Pipeline Version | Date Updated | Documentation Authors | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [GvsJointVariantCalling](https://github.com/broadinstitute/gatk/blob/rc-vs-483-beta-user-wdl/scripts/variantstore/wdl/GvsJointVariantCalling.wdl) | June, 2022 | [Kaylee Mathews](mailto:kmathews@broadinstitute.org) and [Aurora Cremer](mailto:aurora@broadinstitute.org) | If you have questions or feedback, contact the [Broad Variants team](mailto:variants@broadinstitute.org) |
 
![Diagram depicting the Genomic Variant Store workflow. Sample GVCF files are imported into the core data model. A filtering model is trained using Variant Quality Score Recalibration, or VQSR, and then used to extract cohorts and produce sharded joint VCF files. Each step integrates BigQuery and GATK tools.](/scripts/variantstore/genomic-variant-store_diagram.png)


## Introduction to the Genomic Variant Store pipeline

The [Genomic Variant Store (GVS)](https://github.com/broadinstitute/gatk/blob/ah_var_store/scripts/variantstore/gvs-product-sheet.pdf) was developed by the Data Sciences Platform at the Broad Institute of MIT and Harvard as a solution for variant discovery on a large scale. The GVS is powered by BigQuery and creates joint callsets of up to 100,000 human genomes more reliably with decreased time and cost compared to previous solutions.

The [GVS pipeline](https://github.com/broadinstitute/gatk/blob/rc-vs-483-beta-user-wdl/scripts/variantstore/wdl/GvsJointVariantCalling.wdl) is an open-source, cloud-optimized workflow for joint calling at a large scale using the GVS. The workflow takes in up to 100,000 single sample gVCF files and combines them into a variant filtering model driven by machine learning. The model is applied to the data, and a sharded joint VCF with variant calls is output.

The filtering model is based on the [WARP Joint Genotyping workflow](https://github.com/broadinstitute/warp/blob/master/pipelines/broad/dna_seq/germline/joint_genotyping/JointGenotyping.wdl) and uses the [Variant Quality Score Recalibration (VQSR)](https://gatk.broadinstitute.org/hc/en-us/articles/360035531612) technique which uses machine learning to model the technical profile of variants to flag probable artifacts.

## Quickstart table

The following table provides a quick overview of the GVS pipeline features:

| Pipeline features | Description | Source | 
| --- | --- | --- |
| Overall workflow | End-to-end joint calling pipeline that imports samples, trains the filtering model, and extracts VCF files | Code available from [GitHub](https://github.com/broadinstitute/gatk/blob/rc-vs-483-beta-user-wdl/scripts/variantstore/wdl/GvsJointVariantCalling.wdl) |
| Filtering model | Machine-learning powered; uses VQSR and sample annotations | [VQSR](https://gatk.broadinstitute.org/hc/en-us/articles/360035531612) |
| Workflow language | WDL 1.0 | [openWDL](https://github.com/openwdl/wdl) |
| Genomic reference sequence | GRCh38 (hg38) human genome primary sequence | Genome Reference Consortium [GRCh38](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39) |
| Data input file format | File format in which input data is provided | [gVCF](https://gatk.broadinstitute.org/hc/en-us/articles/360035531812) |
| Data output file formats | File formats in which outputs are provided | [VCF](https://gatk.broadinstitute.org/hc/en-us/articles/360035531692) and associated index files; TXT (manifest) |

## Setup

### Getting started with Terra

The GVS workflow is available on Terra, a cloud-based platform for biomedical research. After registering, you can access the [GVS workspace](LINK), which contains instructions for performing joint calling with the GVS using example data uploaded to the workspace. For more information on using the Terra platform, see the [Terra Support Center](https://support.terra.bio/hc/en-us).

If you are new to Terra, you’ll need to [register for a Terra account](https://support.terra.bio/hc/en-us/articles/360028235911). Beta testers of the GVS will be given access to a Broad-managed billing project, but to use Terra for other projects, you’ll need to [set up a billing project](https://support.terra.bio/hc/en-us/articles/360046295092).

The GVS beta workspace in Terra is read-only, so you’ll need to clone the workspace to create a copy where you can upload your own data and run the workflow. For instructions on how to clone workspaces in Terra, see the article [Make your own project workspace](https://support.terra.bio/hc/en-us/articles/360026130851). For a tour of workspaces in Terra, see [Working with workspaces](https://support.terra.bio/hc/en-us/articles/360024743371).

### Workflow requirements

#### Set up BigQuery

To start using the GVS workflow in Terra, you need to create a BigQuery dataset with permissions that will allow Terra to access the dataset by following the instructions below:

**1.** First, **create a dataset in BigQuery** by following the instructions in [Creating datasets](https://cloud.google.com/bigquery/docs/datasets) from the BigQuery documentation.

**2.** Find your Terra proxy group by **selecting the main menu** (three horizontal line icon) at the top left of any Terra page. Next, **click on your account name**, followed by **Profile**. Here, you should see a field called **Proxy Group**, which lists your proxy group underneath. For more information on proxy groups in Terra, see [Pet service accounts and proxy groups](https://support.terra.bio/hc/en-us/articles/360031023592). 

**3.** Grant your proxy group the **BigQuery Data Editor** role on the dataset you created in Step 1 by following the instructions in [Controlling access to datasets](https://cloud.google.com/bigquery/docs/dataset-access-controls) from the BigQuery documentation.

**4.** Grant your proxy group the **BigQuery Data Editor**, **BigQuery Job User**, and **BigQuery Read Session User** roles on the Google Project that contains the dataset you created in Step 1 by following the instructions in [Manage access to projects, folders, and organizations](https://cloud.google.com/iam/docs/granting-changing-revoking-access) from the Google Cloud documentation.

#### Reblocked gVCF files

The GVS workflow takes in reblocked gVCF files as input. If your files are not already reblocked, you can reblock them using the [WARP reblocking workflow](https://github.com/broadinstitute/warp/blob/master/pipelines/broad/dna_seq/germline/joint_genotyping/reblocking/ReblockGVCF.wdl), which is configured in the [ReblockGVCF Terra workspace](https://app.terra.bio/#workspaces/warp-pipelines/ReblockGVCF). For more information about reblocking, check out [WARP Whole Genome and Exome Pipelines Produce Reblocked GVCFs](https://broadinstitute.github.io/warp/blog/tags/reblock/).

#### gVCF annotations

Input gVCF files for the GVS workflow must include the annotations described in the table below:

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

### Inputs

The GVS workflow inputs are described in the sections below and are specified in the Terra GVS beta workspace [workflow configuration](LINK).

#### Sample data input

The GVS workflow takes in reblocked single sample gVCF files and their corresponding index files as `input_vcfs` and `input_vcf_indexes`, respectively. While the GVS workflow takes in up to 100,000 single sample gVCF files as input, only datasets of up to 10,000 files are being used for beta testing.

The [Terra GVS beta workspace](LINK) has been configured with example reblocked gVCF files that you can use to test the workflow. To test the workflow on your own sample data, you’ll need to create an additional data table in your workspace that specifies the locations of the files in the cloud. For step-by-step instructions on creating data tables in Terra, see [How to make a data table from scratch or a template](https://support.terra.bio/hc/en-us/articles/6197368140955). If your data is stored locally, you can upload it to your workspace bucket using the Data Uploader in Terra by following the instructions in the article, [How to use the Data Uploader](https://support.terra.bio/hc/en-us/articles/4419428208411). You’ll also need to ensure that you have a BigQuery dataset and Google project with the correct permissions and that your gVCF files are reblocked, as described in the [Workflow requirements](#workflow-requirements) section above.

#### Input descriptions

The table below describes the GVS workflow input variables:

| Input variable name | Description | Type |
| --- | --- | --- |
| dataset_name | String describing the name of the BigQuery dataset. | String |
| project_id | String describing the name of the Google project that contains the sample data. | String |
| external_sample_names | Strings describing the unique sample IDs in the Terra data model. | Array of strings |
| input_vcfs | Cloud paths to the sample gVCF files. | Array of files |
| input_vcf_indexes | Cloud paths to the sample gVCF index files. | Array of files |
| callset_identifier | String used as the name of the filter model and as the prefix to the names of the BigQuery extract tables and final joint VCF shards; should begin with a letter; valid characters include A-z, 0-9, “.”, “,”, “-“, and “_”. | String |
| extract_output_gcs_dir | Optional string describing the desired cloud path for output files. | String |

## Tasks and tools

Overall, the GVS workflow:

1. Imports sample gVCF files.
1. Trains the filtering model.
1. Extracts VCF files.

The steps of the workflow, in addition to the tools used in each step, are described below.

### 1. Import sample gVCF files

Custom tools used: CreateVariantIngestFiles

This step validates that sample gVCF files contain required annotations and loads samples into BigQuery tables.

### 2. Train the filtering model

GATK tools used: [SplitIntervals (GATK)](https://gatk.broadinstitute.org/hc/en-us/articles/5358914364699), [GatherVcfsCloud (GATK)](https://gatk.broadinstitute.org/hc/en-us/articles/5358884598555), [VariantRecalibrator](https://gatk.broadinstitute.org/hc/en-us/articles/5358906115227), [GatherTranches](https://gatk.broadinstitute.org/hc/en-us/articles/5358889613339)

Custom tools used: ExtractFeatures, CreateFilteringFiles, CreateSiteFilteringFiles

This step splits alternate alleles, calculates annotations to be used for filtering, and creates a non-deterministic filtering model using VQSR and annotations from a random subset of the input samples. SNPs are recalibrated using the annotations `AS_FS`, `AS_ReadPosRankSum`, `AS_MQRankSum`, `AS_QD`, and `AS_SOR`. Indels are recalibrated using the annotations `AS_FS`, `AS_ReadPosRankSum`, `AS_MQRankSum`, `AS_QD`, `AS_SOR`, and `AS_MQ`.

### 3. Extract VCF files

GATK tools used: [SplitIntervals (GATK)](https://gatk.broadinstitute.org/hc/en-us/articles/5358914364699)

Custom tools used: ExtractCohort

This step calculates annotations including allele count (`AC`), allele number (`AN`), and allele frequency (`AF`), and creates and outputs a sharded joint VCF file. The output VCF file includes desired samples, calculated annotations, and flagged probable artifacts. The output VCF file is split so shards do not span multiple chromosomes.

## Outputs

The final outputs of the GVS workflow are described in the table below: 

| Output variable name | Description | Type | 
| ------ | ------ | ------ |
| output_vcfs | Sharded VCF files with variant calls. | Array of VCF files |
| output_vcf_indexes | Sharded VCF index files for the output VCF files. | Array of VCF index files |
| total_vcfs_size_mb | Float describing the total size of the output VCF files in MB. | Float |
| manifest | TXT file listing output file destinations and other metadata. | TXT | 

The GVS workflow outputs a sharded joint VCF file containing filter sites and genotypes flagged as probable artifacts and annotations calculated during the `GvsExtractCallset` subworkflow, including allele count (`AC`), allele number (`AN`), and allele frequency (`AF`). The output VCF file is sharded so that no shards span multiple chromosomes. 

The [Terra GVS beta workspace](LINK) is configured to write the outputs of the workflow back to the `sample_set` data table.

## Feedback

Please help us improve our tools by contacting the [Broad Variants team](mailto:variants@broadinstitute.org) for pipeline-related suggestions or questions.