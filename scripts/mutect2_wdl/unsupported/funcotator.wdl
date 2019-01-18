version 1.0
# Run Funcotator on a set of called variants from Mutect 2.
#
# Description of inputs:
#
#   Required:
#     gatk_docker                    -  GATK Docker image in which to run
#     ref_fasta                      -  Reference FASTA file.
#     ref_fasta_index                -  Reference FASTA file index.
#     ref_fasta_dict                 -  Reference FASTA file sequence dictionary.
#     variant_vcf_to_funcotate       -  Variant Context File (VCF) containing the variants to annotate.
#     variant_vcf_to_funcotate_index -  Index for the Variant Context File (VCF) containing the variants to annotate.
#     reference_version              -  Version of the reference being used.  Either `hg19` or `hg38`.
#     output_file_name               -  Path to desired output file.
#
#   Optional:
#     interval_list                  -  Intervals to be used for traversal.  If specified will only traverse the given intervals.
#     data_sources_tar_gz            -  Path to tar.gz containing the data sources for Funcotator to create annotations.
#     transcript_selection_mode      -  Method of detailed transcript selection.  This will select the transcript for detailed annotation (either `CANONICAL` or `BEST_EFFECT`).
#     transcript_selection_list      -  Set of transcript IDs to use for annotation to override selected transcript.
#     annotation_defaults            -  Annotations to include in all annotated variants if the annotation is not specified in the data sources (in the format <ANNOTATION>:<VALUE>).  This will add the specified annotation to every annotated variant if it is not already present.
#     annotation_overrides           -  Override values for annotations (in the format <ANNOTATION>:<VALUE>).  Replaces existing annotations of the given name with given values.
#     gatk4_jar_override             -  Override Jar file containing GATK 4.0.  Use this when overriding the docker JAR or when using a backend without docker.
#     funcotator_extra_args          -  Extra command-line arguments to pass through to Funcotator.  (e.g. " --exclude-field foo_field --exclude-field bar_field ")
#
# This WDL needs to decide whether to use the ``gatk_jar`` or ``gatk_jar_override`` for the jar location.  As of cromwell-0.24,
# this logic *must* go into each task.  Therefore, there is a lot of duplicated code.  This allows users to specify a jar file
# independent of what is in the docker file.  See the README.md for more info.
#
# NOTE: This only does VCF output right now!
#
workflow Funcotator {

  input {
    String gatk_docker
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File variant_vcf_to_funcotate
    File variant_vcf_to_funcotate_index
    String reference_version
    String output_file_base_name

    File? interval_list
    File? data_sources_tar_gz
    String? transcript_selection_mode
    Array[String]? transcript_selection_list
    Array[String]? annotation_defaults
    Array[String]? annotation_overrides
    File? gatk4_jar_override

    String? funcotator_extra_args
  }

  call Funcotate {
    input:
      ref_fasta                 = ref_fasta,
      ref_fasta_index           = ref_fasta_index,
      ref_dict                  = ref_dict,
      input_vcf                 = variant_vcf_to_funcotate,
      input_vcf_idx             = variant_vcf_to_funcotate_index,
      reference_version         = reference_version,
      interval_list             = interval_list,
      output_file_base_name     = output_file_base_name,
      output_format             = "VCF",
      data_sources_tar_gz       = data_sources_tar_gz,
      transcript_selection_mode = transcript_selection_mode,
      transcript_selection_list = transcript_selection_list,
      annotation_defaults       = annotation_defaults,
      annotation_overrides      = annotation_overrides,
      gatk_override             = gatk4_jar_override,
      gatk_docker               = gatk_docker,
      extra_args                = funcotator_extra_args
  }

  output {
    File funcotated_vcf_out = Funcotate.funcotated_vcf
    File funcotated_vcf_out_idx = Funcotate.funcotated_vcf_index
  }
}


task Funcotate {
  input {
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File input_vcf
    File input_vcf_idx
    String reference_version
    String output_file_base_name
    String output_format
    Boolean compress = true

    File? data_sources_tar_gz
    String? transcript_selection_mode
    Array[String]? transcript_selection_list
    Array[String]? annotation_defaults
    Array[String]? annotation_overrides
    Boolean filter_funcotations = false
    File? interval_list

    String? extra_args

    # runtime args
    String gatk_docker
    Int machine_memory = 3000
    Int preemptible_attempts = 3
    Int additional_disk = 0
    Int cpu_threads = 1
    Boolean use_ssd = false
    File? gatk_override
  }

  String output_vcf = output_file_base_name + if compress then ".vcf.gz" else".vcf"
  String output_vcf_index = output_vcf + if compress then ".tbi" else ".idx"

  Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")
  # We need to compute and multiply data source size because uncompressing it can get really big.
  # The default data source is 14 GB, so we multiply that by 10 for uncompressed size
  Int compressed_data_size = if defined(data_sources_tar_gz) then ceil(size(data_sources_tar_gz, "GB") * 10) else 140
  Int disk_size = ceil(size(input_vcf, "GB") + ref_size) + compressed_data_size + 20 + additional_disk
  String disk_type = if use_ssd then "SSD" else "HDD"

  # Mem is in units of GB but our command and memory runtime values are in MB
  Int command_memory = machine_memory - 1000

  command <<<
    set -e
    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

    DATA_SOURCES_FOLDER="$PWD/datasources_dir"
    DATA_SOURCES_TAR_GZ=~{default="" data_sources_tar_gz}
    if [[ ! -e $DATA_SOURCES_TAR_GZ ]] ; then
      DOWNLOADED_DATASOURCES_NAME="downloaded_datasources.tar.gz"
      gatk FuncotatorDataSourceDownloader --germline --output $DOWNLOADED_DATASOURCES_NAME
      DATA_SOURCES_TAR_GZ=$DOWNLOADED_DATASOURCES_NAME
    fi
    # Extract provided the tar.gz:
    mkdir $DATA_SOURCES_FOLDER
    tar zxvf $DATA_SOURCES_TAR_GZ -C $DATA_SOURCES_FOLDER --strip-components 1

    gatk --java-options "-Xmx~{command_memory}m" Funcotator \
      --data-sources-path $DATA_SOURCES_FOLDER \
      --ref-version ~{reference_version} \
      --output-file-format ~{output_format} \
      -R ~{ref_fasta} \
      -V ~{input_vcf} \
      -O ~{output_vcf} \
      ~{true="-L" false="" defined(interval_list)} ~{default="" interval_list} \
      ~{true="--transcript-selection-mode" false="" defined(transcript_selection_mode)} ~{default="" transcript_selection_mode} \
      ~{true="--transcript-list" false="" defined(transcript_selection_list)} ~{default="" sep=" --transcript-list " transcript_selection_list} \
      ~{true="--annotation-default" false="" defined(annotation_defaults)} ~{default="" sep=" --annotation-default " annotation_defaults} \
      ~{true="--annotation-override" false="" defined(annotation_overrides)} ~{default="" sep=" --annotation-override " annotation_overrides} \
      ~{true="--remove-filtered-variants" false="" filter_funcotations} \
      ~{default="" extra_args}
  >>>

  runtime {
    docker: gatk_docker
    memory: "~{machine_memory} MB"
    disks: "local-disk ~{disk_size} ~{disk_type}"
    preemptible: preemptible_attempts
    cpu: cpu_threads
  }

  output {
    File funcotated_vcf = "~{output_vcf}"
    File funcotated_vcf_index = "~{output_vcf_index}"
  }
 }