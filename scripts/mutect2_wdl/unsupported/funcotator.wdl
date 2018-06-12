# Run Funcotator on a set of called variants in a VCF file.
# Variants are assumed to have been called by Mutect 2.
#
# =============================================================================
# Description of Inputs:
#
#   Required:
#     gatk_docker                                         - GATK Docker image in which to run
#
#     ref_fasta                                           - Reference FASTA file.
#     ref_fasta_index                                     - Reference FASTA file index.
#     ref_fasta_dict                                      - Reference FASTA file sequence dictionary.
#     variant_vcf_to_funcotate                            - Variant Context File (VCF) containing the variants to annotate.
#     variant_vcf_to_funcotate_index                      - Index for the Variant Context File (VCF) containing the variants to annotate.
#     reference_version                                   - Version of the reference being used.  Either `hg19` or `hg38`.
#     output_file_base_name                               - Path to desired output file without extension.
#     output_format                                       - Format for the output file.  Either`VCF` or `MAF`.
#
#   Optional:
#     gatk4_jar_override                                  - Override Jar file containing GATK 4.0.  Use this when overriding the docker JAR or when using a backend without docker.
#
#     data_sources_tar_gz                                 - Path to tar.gz containing the data sources for Funcotator to create annotations.
#     transcript_selection_mode                           - Method of detailed transcript selection.  This will select the transcript for detailed annotation (either `CANONICAL` or `BEST_EFFECT`).
#     transcript_selection_list                           - Set of transcript IDs to use for annotation to override selected transcript.
#     annotation_defaults                                 - Annotations to include in all annotated variants if the annotation is not specified in the data sources (in the format <ANNOTATION>:<VALUE>).  This will add the specified annotation to every annotated variant if it is not already present.
#     annotation_overrides                                - Override values for annotations (in the format <ANNOTATION>:<VALUE>).  Replaces existing annotations of the given name with given values.
#     filter_funcotations                                 - True if filtered variants should be funcotated and included in output.
#     interval_list                                       - Intervals to be used for traversal.  If specified will only traverse the given intervals.
#     compress_funcotator_output                          - Set to true to compress the output file (only applicable for VCF files).
#
#     allow_hg19_gencode_b37_contig_matching              - Set to true to perform "fuzzy matching" should be used between contigs using references hg19 and b37.  This will allow B37 contigs to match hg19 references and visa-versa.
#     allow_hg19_gencode_b37_contig_matching_force_unsafe - Use in conjunction with `allow_hg19_gencode_b37_contig_matching`.  If you also this flag to true, no check that your input reference is b37 is actually performed.  This will allow ANY reference to potentially match ANY VCF.  Use at your own peril.
#
#     funcotator_extra_args                               - Extra command-line arguments to pass through to Funcotator.
#
# =============================================================================
# Outputs:
#   funcotated_output
#
# This WDL needs to decide whether to use the ``gatk_jar`` or ``gatk_jar_override`` for the jar location.  As of cromwell-0.24,
# this logic *must* go into each task.  Therefore, there is a lot of duplicated code.  This allows users to specify a jar file
# independent of what is in the docker file.  See the README.md for more info.
#
workflow Funcotator {
    String gatk_docker

    File ref_fasta
    File ref_fasta_index
    File ref_fasta_dict
    File variant_vcf_to_funcotate
    File variant_vcf_to_funcotate_index
    String reference_version
    String output_file_base_name
    String output_format

    File? gatk4_jar_override

    File? data_sources_tar_gz
    String? transcript_selection_mode
    Array[String]? transcript_selection_list
    Array[String]? annotation_defaults
    Array[String]? annotation_overrides
    Boolean? filter_funcotations
    File? interval_list
    Boolean? compress_funcotator_output

    Boolean? allow_hg19_gencode_b37_contig_matching
    Boolean? allow_hg19_gencode_b37_contig_matching_force_unsafe
    Int? lookahead_cache_bp

    String? funcotator_extra_args

    call Funcotate {
        input:
            ref_fasta                                       = ref_fasta,
            ref_fasta_index                                 = ref_fasta_index,
            ref_dict                                        = ref_fasta_dict,
            input_vcf                                       = variant_vcf_to_funcotate,
            input_vcf_idx                                   = variant_vcf_to_funcotate_index,
            reference_version                               = reference_version,
            output_file_base_name                           = output_file_base_name,
            output_format                                   = output_format,

            data_sources_tar_gz                             = data_sources_tar_gz,
            transcript_selection_mode                       = transcript_selection_mode,
            transcript_selection_list                       = transcript_selection_list,
            annotation_defaults                             = annotation_defaults,
            annotation_overrides                            = annotation_overrides,
            filter_funcotations                             = filter_funcotations,
            interval_list                                   = interval_list,
            compress                                        = compress_funcotator_output,

            allow_hg19_gencode_b37_contig_matching          = allow_hg19_gencode_b37_contig_matching,
            allow_hg19_gencode_b37_contig_matching_force_unsafe = allow_hg19_gencode_b37_contig_matching_force_unsafe,
            lookahead_cache_bp                              = lookahead_cache_bp,

            extra_args                                      = funcotator_extra_args,

            gatk_override                                   = gatk4_jar_override,
            gatk_docker                                     = gatk_docker
    }

    output {
        File funcotated_output     = Funcotate.funcotated_output
        File funcotated_output_idx = Funcotate.funcotated_index
    }
}


task Funcotate {

    # =======================================
    # Required Inputs:
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File input_vcf
    File input_vcf_idx
    String reference_version
    String output_file_base_name
    String output_format

    # =======================================
    # Optional Inputs:
    File? data_sources_tar_gz
    String? transcript_selection_mode
    Array[String]? transcript_selection_list
    Array[String]? annotation_defaults
    Array[String]? annotation_overrides
    Boolean? filter_funcotations
    File? interval_list
    Boolean? compress

    Boolean? allow_hg19_gencode_b37_contig_matching
    Boolean? allow_hg19_gencode_b37_contig_matching_force_unsafe
    Int? lookahead_cache_bp

    String? extra_args

    # =======================================
    # Process input args:

    String output_file = output_file_base_name + if output_format == "VCF" then if compress then ".vcf.gz" else ".vcf" else ".maf"
    String output_index = if output_format == "VCF" then output_file +  if compress then ".tbi" else ".idx" else ""

    String transcript_selection_list_arg = if defined(transcript_selection_list) then " --transcript-list " else ""
    String annotation_def_arg = if defined(annotation_defaults) then " --annotation-default " else ""
    String annotation_over_arg = if defined(annotation_overrides) then " --annotation-override " else ""
    String filter_funcotations_args = if (filter_funcotations) then " --remove-filtered-variants " else ""
    String interval_list_arg = if defined(interval_list) then " -L " else ""

    String allow_hg19_gencode_b37_contig_matching_arg = if defined(allow_hg19_gencode_b37_contig_matching) then "--allow-hg19-gencode-b37-contig-matching" else ""
    String allow_hg19_gencode_b37_contig_matching_force_unsafe_arg = if defined(allow_hg19_gencode_b37_contig_matching_force_unsafe) then "--allow-hg19-gencode-b37-contig-matching-override" else ""
    String lookahead_cache_bp_arg = if defined(lookahead_cache_bp) then "--lookahead-cache-bp" else ""

    String extra_args_arg = select_first([extra_args, ""])

    # =======================================
    # Runtime setup:

    String gatk_docker
    File? gatk_override
    Int? mem
    Int? preemptible_attempts
    Int? disk_space_gb
    Int? cpu

    Boolean use_ssd = false

    # This should be updated when a new version of the data sources is released
    # TODO: Make this dynamically chosen in the command.
    String default_datasources_version = "funcotator_dataSources.v1.3.20180531"

    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 3000
    # WARNING: In the workflow, you should calculate the disk space as an input to this task (disk_space_gb).  Please see [TODO: Link from Jose] for examples.
    Int default_disk_space_gb = 100

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem *1000 else default_ram_mb
    Int command_mem = machine_mem - 1000

    command <<<
        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        DATA_SOURCES_TAR_GZ=${data_sources_tar_gz}
        if [[ ! -e $DATA_SOURCES_TAR_GZ ]] ; then
            # We have to download the data sources:
            echo "Data sources gzip does not exist: $DATA_SOURCES_TAR_GZ"
            echo "Downloading default data sources..."
            wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/funcotator/${default_datasources_version}.tar.gz
            tar -zxf ${default_datasources_version}.tar.gz
            DATA_SOURCES_FOLDER=${default_datasources_version}
        else
            # Extract the tar.gz:
            mkdir datasources_dir
            tar zxvf ${data_sources_tar_gz} -C datasources_dir --strip-components 1
            DATA_SOURCES_FOLDER="$PWD/datasources_dir"
        fi

        gatk --java-options "-Xmx${command_mem}m" Funcotator \
            --data-sources-path $DATA_SOURCES_FOLDER \
            --ref-version ${reference_version} \
            --output-file-format ${output_format} \
            -R ${ref_fasta} \
            -V ${input_vcf} \
            -O ${output_file} \
            ${interval_list_arg} ${default="" interval_list} \
            ${"--transcript-selection-mode " + transcript_selection_mode} \
            ${transcript_selection_list_arg}${default="" sep=" --transcript-list " transcript_selection_list} \
            ${annotation_def_arg}${default="" sep=" --annotation-default " annotation_defaults} \
            ${annotation_over_arg}${default="" sep=" --annotation-override " annotation_overrides} \
            ${filter_funcotations_args} \
            ${allow_hg19_gencode_b37_contig_matching_arg} \
            ${allow_hg19_gencode_b37_contig_matching_force_unsafe_arg} \
            ${lookahead_cache_bp_arg}${lookahead_cache_bp} \
            ${extra_args_arg}
    >>>

    runtime {
        docker: gatk_docker
        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + if use_ssd then " SSD" else " HDD"
        preemptible: select_first([preemptible_attempts, 3])
        cpu: select_first([cpu, 1])
    }

    output {
        File funcotated_output = "${output_file}"
        File funcotated_index  = "${output_index}"
    }
 }