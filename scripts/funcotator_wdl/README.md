# Running the Funcotator WDL

## Background Information
Funcotator (**FUNC**tional ann**OTATOR**) is a functional annotation tool in the core GATK toolset and was designed to handle both somatic and germline use cases. It analyzes given variants for their function (as retrieved from a set of data sources) and produces the analysis in a specified output file.  Funcotator reads in a VCF file, labels each variant with one of twenty-three distinct variant classifications, produces gene information (e.g. affected gene, predicted variant amino acid sequence, etc.), and associations to information in datasources. Default supported datasources include GENCODE (gene information and protein change prediction), dbSNP, gnomAD, and COSMIC (among others). The corpus of datasources is extensible and user-configurable and includes cloud-based datasources supported with Google Cloud Storage. Funcotator produces either a Variant Call Format (VCF) file (with annotations in the INFO field) or a Mutation Annotation Format (MAF) file.

Funcotator allows the user to add their own annotations to variants based on a set of data sources.  Each data source can be customized to annotate a variant based on several matching criteria.  This allows a user to create their own custom annotations easily, without modifying any Java code.

## Setup 

To run the Funcotator WDL you must have access to a cromwell server that can run your job.

Once your cromwell instance is active, you will need to generate input arguments to pass to funcotator.wdl.  These arguments re passed in as a JSON file (see below for a non-working example).

Once a JSON file has been created you can submit your job to a cromwell server directly (i.e. using a tool such as Cromshell) or through Terra/FireCloud.

## WDL Input Parameters

The input parameters to the Funcotator WDL are as follows:

### Required Inputs:
String gatk_docker                  - GATK Docker image in which to run

File ref_fasta                      - Reference FASTA file.

File ref_fasta_index                - Reference FASTA file index.

File ref_fasta_dict                 - Reference FASTA file sequence dictionary.

File variant_vcf_to_funcotate       - Variant Context File (VCF) containing the variants to annotate.

File variant_vcf_to_funcotate_index - Index file corresponding to the input Variant Context File (VCF) containing the variants to annotate.

String reference_version            - Version of the reference being used.  Either `hg19` or `hg38`.

String output_file_name             - Path to desired output file.

String output_format                - Output file format (either VCF or MAF).

Boolean compress				    - Whether to compress the resulting output file.

Boolean use_gnomad                  - If true, will enable the gnomAD data sources in the data source tar.gz, if they exist.


### Optional Inputs:
File? interval_list                      - Intervals to be used for traversal.  If specified will only traverse the given intervals.

File? data_sources_tar_gz                - Path to tar.gz containing the data sources for Funcotator to create annotations.

String? transcript_selection_mode        - Method of detailed transcript selection.  This will select the transcript for detailed annotation (either `CANONICAL` or `BEST_EFFECT`).

Array[String]? transcript_selection_list - Set of transcript IDs to use for annotation to override selected transcript.

Array[String]? annotation_defaults       - Annotations to include in all annotated variants if the annotation is not specified in the data sources (in the format <ANNOTATION>:<VALUE>).  This will add the specified annotation to every annotated variant if it is not already present.

Array[String]? annotation_overrides      - Override values for annotations (in the format <ANNOTATION>:<VALUE>).  Replaces existing annotations of the given name with given values.

File? gatk4_jar_override                 - Override Jar file containing GATK 4.0.  Use this when overriding the docker JAR or when using a backend without docker.

String? funcotator_extra_args            - Extra command-line arguments to pass through to Funcotator.  (e.g. " --exclude-field foo_field --exclude-field bar_field ")

## Example JSON File (Non-Working)

The follwing is an example of a JSON input file.  It will not work as-is but is provided as a starting point for you to create your own input file:

```
{
  "Funcotator.gatk_docker": "broadinstitute/gatk:latest",
  
  "Funcotator.ref_fasta": "gs://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta",
  "Funcotator.ref_fasta_index": "gs://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
  "Funcotator.ref_dict": "gs://broad-references/hg38/v0/Homo_sapiens_assembly38.dict",
  
  "Funcotator.reference_version": "hg38",
  "Funcotator.output_format": "VCF",

  "Funcotator.compress": "false",
  "Funcotator.use_gnomad": "false",
  "Funcotator.data_sources_tar_gz": "gs://broad-public-datasets/funcotator/funcotator_dataSources.v1.6.20190124s.tar.gz",

  "Funcotator.variant_vcf_to_funcotate": "variants.vcf",
  "Funcotator.variant_vcf_to_funcotate_index": "variants.vcf.idx",
  
  "Funcotator.output_file_base_name": "variants.funcotated"
}
```

## Further Information
 - https://software.broadinstitute.org/gatk/documentation/article?id=11193
 - https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.4.0/org_broadinstitute_hellbender_tools_funcotator_Funcotator.php
 
 