# Genomic Variant Store workflow overview
 
| Workflow Version | Date Updated | Documentation Authors | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [GvsJointVariantCalling](https://github.com/broadinstitute/gatk/blob/ah_var_store/scripts/variantstore/wdl/GvsJointVariantCalling.wdl) | June, 2022 | [Kaylee Mathews](mailto:kmathews@broadinstitute.org) and [Aurora Cremer](mailto:aurora@broadinstitute.org) | If you have questions or feedback, contact the [Broad Variants team](mailto:variants@broadinstitute.org) |
 
![Diagram depicting the Genomic Variant Store workflow. Sample GVCF files are imported into the core data model. A filtering model is trained using Variant Extract-Train-Score, or VETS, and then used to extract cohorts and produce sharded joint VCF files. Each step integrates BigQuery and GATK tools.](./genomic-variant-store_diagram.png)

## Introduction to the Genomic Variant Store workflow

The [Genomic Variant Store (GVS)](../gvs-product-sheet.pdf) was developed by the Data Sciences Platform at the Broad Institute of MIT and Harvard as a solution for variant discovery on a large scale. The GVS is powered by BigQuery and creates large joint callsets more reliably with decreased time and cost compared to previous solutions.

The [GVS workflow](https://github.com/broadinstitute/gatk/blob/ah_var_store/scripts/variantstore/wdl/GvsJointVariantCalling.wdl) is an open-source, cloud-optimized workflow for joint calling at a large scale using the Terra platform. The workflow takes in single sample GVCF files, loads them into [BigQuery](https://cloud.google.com/bigquery/docs) tables, and combines them into a variant filtering model driven by machine learning. The model is uploaded back into BigQuery and applied to the data. The workflow produces sharded joint VCF files with indices, a manifest file, and metrics.

The filtering approach follows the new best practice to use the GATK Variant Extract-Train-Score (VETS) toolchain with an isolation-forest model to flag probable artifacts.

---

**Want to try out the GVS workflow?**     
To get started using the GVS workflow in Terra with example data, follow the instructions in the [GVS Beta Quickstart](./gvs-quickstart.md). To run the GVS workflow on your own sample data, follow the instructions in the tutorial, [Upload data to Terra and run the GVS workflow](./run-your-own-samples.md).

---

## Quickstart table

The following table provides a quick overview of the GVS workflow features:

| Workflow features          | Description                                                                                                | Source                                                                                                                                     | 
|----------------------------|------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------|
| Overall workflow           | End-to-end joint calling workflow that imports samples, trains the filtering model, and extracts VCF files | Code available from [GitHub](https://github.com/broadinstitute/gatk/blob/ah_var_store/scripts/variantstore/wdl/GvsJointVariantCalling.wdl) |
| Filtering model            | Isolation-forest machine learning model; uses GATK VETS toolchain and sample annotations                   | [VETS](https://github.com/broadinstitute/gatk/blob/ah_var_store/scripts/variantstore/docs/release_notes/VETS_Release.pdf)                  |
| Workflow language          | WDL 1.0                                                                                                    | [openWDL](https://github.com/openwdl/wdl)                                                                                                  |
| Genomic reference sequence | GRCh38 (hg38) human genome primary sequence                                                                | Genome Reference Consortium [GRCh38](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39)                                               |
| Data input file format     | File format in which input data is provided                                                                | [Reblocked](https://broadinstitute.github.io/warp/blog/Nov21_ReblockedGVCF) [GVCF](https://gatk.broadinstitute.org/hc/en-us/articles/360035531812)                                                       |
| Data output file formats   | File formats in which outputs are provided                                                                 | [VCF](https://gatk.broadinstitute.org/hc/en-us/articles/360035531692) and associated index files; TXT (manifest)                           |

## Setup

### Workflow requirements

#### Terra, Google Cloud, and BigQuery

The GVS workflow needs to be run in [Terra](https://app.terra.bio/), a cloud-based platform for biomedical research. The workflow also requires that you have a Terra account, billing project, and BigQuery dataset with permissions that allow Terra to access it. For step-by-step instructions for setting up these requirements, see the [GVS Beta Quickstart](./gvs-quickstart.md).

#### Input GVCF files

The GVS workflow takes in reblocked GVCF files as input. If your GVCF files have not been reblocked, you can reblock them using the [WARP reblocking workflow](https://github.com/broadinstitute/warp/blob/ReblockGVCF_v2.1.12/pipelines/broad/dna_seq/germline/joint_genotyping/reblocking/ReblockGVCF.wdl), which is configured in the GVS workspace for your convenience. You can read instructions on running Reblocking here: [ReblockGVCF Terra workspace](https://app.terra.bio/#workspaces/warp-pipelines/ReblockGVCF). For more information about reblocking, check out [WARP Whole Genome and Exome Pipelines Produce Reblocked GVCFs](https://broadinstitute.github.io/warp/blog/tags/reblock/).

The workflow also requires specific annotations in input GVCF files, which are described in the tutorial, [Upload data to Terra and run the GVS workflow](./run-your-own-samples.md).

### Inputs

The GVS workflow inputs are described in the sections below and are specified in the Terra GVS beta workspace [workflow configuration](https://app.terra.bio/#workspaces/gvs-prod/Genomic_Variant_Store_Beta/workflows/gvs-prod/GvsBeta).

The workflow is configured to use hg38 (aka GRCh38) as the reference genome.

#### Sample data

The GVS workflow takes in reblocked single sample GVCF files and their corresponding index files from the user-indicated columns of the `sample` data table. While the GVS workflow has been tested with 410,000 single sample GVCF files as input, only datasets of up to 25,000 samples are being used for beta testing.

The [GVS beta workspace](https://app.terra.bio/#workspaces/gvs-prod/Genomic_Variant_Store_Beta) has been configured with example reblocked GVCF files that you can use to test the workflow.

#### Input descriptions

The table below describes the GVS workflow input variables:

| Input variable name         | Description                                                                                                                                                                                                                                                         | Type    |
|-----------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------|
| dataset_name                | Name of the BigQuery dataset used to hold input samples, filtering model data, and other tables created during the workflow.                                                                                                                                        | String  |
| project_id                  | Name of the Google project that contains the BigQuery dataset.                                                                                                                                                                                                      | String  |
| call_set_identifier         | Used to name the filter model, BigQuery extract tables, and final joint VCF shards; should begin with a letter; valid characters include A-z, 0-9, “.”, “,”, “-”, and “_”.                                                                                          | String  |
| vcf_files_column_name       | Name of the column in the `sample` workspace data table with the reblocked GVCF file gcs paths.                                                                                                                                                                     | String  |
| vcf_index_files_column_name | Name of the column in the `sample` workspace data table with the reblocked GVCF index file gcs paths.                                                                                                                                                               | String  |
| sample_id_column_name       | Name of the column in the `sample` workspace data table with the desired sample identifiers. These IDs will be used to name the samples in the VCFs. They **MUST** be unique.                                                                                       | String  |
| extract_output_gcs_dir      | Desired cloud path for output files. If unspecified, the workflow will write outputs to the location `gs://<workspace bucket>/output_vcfs/by_submission_id/<submission id>`. The actual output location will be specified by the workflow output `output_gcs_path`. | String  |
| use_classic_VQSR            | Optional; defaults to false since September 1, 2023.                                                                                                                                                                                                                | Boolean |
| billing_project_id          | Optional; Google project ID to charge for the egress of the GVCFs and index files, useful if the bucket the GVCFs are in has "requester pays" enabled                                                                                                               | String  |

## Tasks and tools

Overall, the GVS workflow:

1. Imports the GVCFs into BigQuery.
2. Trains the filtering model.
3. Extracts VCF files.

The steps of the workflow, in addition to the GATK tools used in each step, are described below. Custom tools are used throughout the workflow to read from and write to the BigQuery dataset.

### 1. Import sample GVCF files

This step validates that sample GVCF files contain required annotations and loads samples into BigQuery tables. 

### 2. Train the filtering model

GATK tools used: [SplitIntervals](https://gatk.broadinstitute.org/hc/en-us/articles/5358914364699), [GatherVcfsCloud](https://gatk.broadinstitute.org/hc/en-us/articles/5358884598555), [ExtractVariantAnnotations](https://gatk.broadinstitute.org/hc/en-us/articles/9570430059803-ExtractVariantAnnotations-BETA-), [TrainVariantAnnotationsModel](https://gatk.broadinstitute.org/hc/en-us/articles/9570371970587-TrainVariantAnnotationsModel-BETA-), [ScoreVariantAnnotations](https://gatk.broadinstitute.org/hc/en-us/articles/9570326304155-ScoreVariantAnnotations-BETA-), and [GatherTranches](https://gatk.broadinstitute.org/hc/en-us/articles/5358889613339) (only used for VQSR).

This step splits alternate alleles, calculates annotations to be used for filtering, and creates a filtering model using VETS and annotations from the input samples. To learn more about the VETS release in GVS see this [release note](https://github.com/broadinstitute/gatk/blob/ah_var_store/scripts/variantstore/docs/release_notes/VETS_Release.pdf). 

SNPs and Indels are filtered using the annotations `AS_FS`, `AS_ReadPosRankSum`, `AS_MQRankSum`, `AS_QD`, `AS_SOR` and `AS_MQ`. The annotations are not output into the final callset.

### 3. Extract VCF files

GATK tools used: [SplitIntervals](https://gatk.broadinstitute.org/hc/en-us/articles/5358914364699)

This step calculates annotations including allele count (`AC`), allele number (`AN`), and allele frequency (`AF`), and creates and outputs a sharded joint VCF file. The output VCF file includes desired samples, calculated annotations, and flagged probable artifacts. The output VCF file is split so shards do not span multiple chromosomes.

## Outputs

The final outputs of the GVS workflow are described in the table below: 

| Output variable name | Description | Type | 
| ------ | ------ | ------ |
| output_vcfs | Sharded VCF files with variant calls. | Array of VCF files |
| output_vcf_indexes | Sharded VCF index files for the output VCF files. | Array of VCF index files |
| total_vcfs_size_mb | Float describing the total size of the output VCF files in MB. | Float |
| manifest | TXT file listing output file destinations and other metadata. | TXT | 

The GVS workflow outputs a sharded joint VCF file containing filter sites and genotypes flagged as probable artifacts and annotations calculated during the `GvsExtractCallset` subworkflow, including allele count (`AC`), allele number (`AN`), and allele frequency (`AF`). The output VCF file is sharded so that no shards span multiple chromosomes. See [GVS outputs documentation](./gvs-outputs.md) for more information on the expected number of shards.

The [GVS beta workspace](https://app.terra.bio/#workspaces/gvs-prod/Genomic_Variant_Store_Beta) is configured to write the outputs of the workflow to the Google bucket specified in the inputs at runtime.

## Citing the GVS workflow

If you use plan to publish data analyzed using the GVS workflow, please cite the [GVS beta workspace](https://app.terra.bio/#workspaces/gvs-prod/Genomic_Variant_Store_Beta).

Details on citing Terra workspaces can be found here: [How to cite Terra](https://support.terra.bio/hc/en-us/articles/360035343652)

Data Sciences Platform, Broad Institute (*Year, Month Day that the workspace was last modified*) gvs-prod/Genomic_Variant_Store_Beta [workspace] Retrieved *Month Day, Year that workspace was retrieved*, https://app.terra.bio/#workspaces/gvs-prod/Genomic_Variant_Store_Beta

### License
**Copyright Broad Institute, 2024 | Apache**
The workflow script is released under the Apache License, Version 2.0 (full license text at https://github.com/broadinstitute/gatk/blob/master/LICENSE.TXT). Note however that the programs called by the scripts may be subject to different licenses. Users are responsible for checking that they are authorized to run all programs before running these tools.  

## Feedback

Please help us improve our tools by contacting the [Broad Variants team](mailto:variants@broadinstitute.org) for workflow-related suggestions or questions.
