# Run Funcotator on a set of called variants.
#
# Description of inputs:
#
#   Required:
#     String gatk_docker                       - GATK Docker image in which to run
#     File ref_fasta                           - Reference FASTA file.
#     File ref_fasta_index                     - Reference FASTA file index.
#     File ref_fasta_dict                      - Reference FASTA file sequence dictionary.
#     File variant_vcf_to_funcotate            - Variant Context File (VCF) containing the variants to annotate.
#     File variant_vcf_to_funcotate_index      - Index file corresponding to the input Variant Context File (VCF) containing the variants to annotate.
#     String reference_version                 - Version of the reference being used.  Either `hg19` or `hg38`.
#     String output_file_name                  - Path to desired output file.
#     String output_format                     - Output file format (either VCF or MAF).
#     Boolean compress				           - Whether to compress the resulting output file.
#     Boolean use_gnomad                       - If true, will enable the gnomAD data sources in the data source tar.gz, if they exist.
#
#   Optional:
#     File? interval_list                      - Intervals to be used for traversal.  If specified will only traverse the given intervals.
#     File? data_sources_tar_gz                - Path to tar.gz containing the data sources for Funcotator to create annotations.
#     String? transcript_selection_mode        - Method of detailed transcript selection.  This will select the transcript for detailed annotation (either `CANONICAL` or `BEST_EFFECT`).
#     Array[String]? transcript_selection_list - Set of transcript IDs to use for annotation to override selected transcript.
#     Array[String]? annotation_defaults       - Annotations to include in all annotated variants if the annotation is not specified in the data sources (in the format <ANNOTATION>:<VALUE>).  This will add the specified annotation to every annotated variant if it is not already present.
#     Array[String]? annotation_overrides      - Override values for annotations (in the format <ANNOTATION>:<VALUE>).  Replaces existing annotations of the given name with given values.
#     File? gatk4_jar_override                 - Override Jar file containing GATK 4.0.  Use this when overriding the docker JAR or when using a backend without docker.
#     String? funcotator_extra_args            - Extra command-line arguments to pass through to Funcotator.  (e.g. " --exclude-field foo_field --exclude-field bar_field ")
#     String? gcs_project_for_requester_pays   - GCS Project specifier if the input data are in a Requester Pays Bucket.
#
# This WDL needs to decide whether to use the ``gatk_jar`` or ``gatk_jar_override`` for the jar location.  As of cromwell-0.24,
# this logic *must* go into each task.  Therefore, there is a lot of duplicated code.  This allows users to specify a jar file
# independent of what is in the docker file.  See the README.md for more info.
#
workflow Funcotator {
    String gatk_docker
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File variant_vcf_to_funcotate
    File variant_vcf_to_funcotate_index
    String reference_version
    String output_file_base_name
    String output_format
    Boolean compress
    Boolean use_gnomad

    File? interval_list
    File? data_sources_tar_gz
    String? transcript_selection_mode
    Array[String]? transcript_selection_list
    Array[String]? annotation_defaults
    Array[String]? annotation_overrides
    String? funcotator_extra_args
    String? gcs_project_for_requester_pays

    File? gatk4_jar_override

    call Funcotate {
        input:
            gatk_docker                    = gatk_docker,
            ref_fasta                      = ref_fasta,
            ref_fasta_index                = ref_fasta_index,
            ref_dict                       = ref_dict,
            input_vcf                      = variant_vcf_to_funcotate,
            input_vcf_idx                  = variant_vcf_to_funcotate_index,
            reference_version              = reference_version,
            output_file_base_name          = output_file_base_name,
            output_format                  = output_format,
            compress                       = compress,
            use_gnomad                     = use_gnomad,

            interval_list                  = interval_list,
            data_sources_tar_gz            = data_sources_tar_gz,
            transcript_selection_mode      = transcript_selection_mode,
            transcript_selection_list      = transcript_selection_list,
            annotation_defaults            = annotation_defaults,
            annotation_overrides           = annotation_overrides,
            extra_args                     = funcotator_extra_args,
            gcs_project_for_requester_pays = gcs_project_for_requester_pays,

            gatk_override                  = gatk4_jar_override
    }

    output {
        File funcotated_file_out = Funcotate.funcotated_output_file
        File funcotated_file_out_idx = Funcotate.funcotated_output_file_index
    }
}

################################################################################

task Funcotate {

    # ==============
    # Inputs
    File ref_fasta
    File ref_fasta_index
    File ref_dict

    File input_vcf
    File input_vcf_idx

    String reference_version

    String output_file_base_name
    String output_format

    Boolean compress
    Boolean use_gnomad

    # This should be updated when a new version of the data sources is released
    # TODO: Make this dynamically chosen in the command.
    File? data_sources_tar_gz = "gs://broad-public-datasets/funcotator/funcotator_dataSources.v1.6.20190124s.tar.gz"

    String? control_id
    String? case_id
    String? sequencing_center
    String? sequence_source
    String? transcript_selection_mode
    File? transcript_selection_list
    Array[String]? annotation_defaults
    Array[String]? annotation_overrides
    Array[String]? funcotator_excluded_fields
    Boolean? filter_funcotations
    File? interval_list

    String? extra_args

    String? gcs_project_for_requester_pays

    # ==============
    # Process input args:

    String output_maf = output_file_base_name + ".maf"
    String output_maf_index = output_maf + ".idx"

    String output_vcf = output_file_base_name + if compress then ".vcf.gz" else ".vcf"
    String output_vcf_idx = output_vcf +  if compress then ".tbi" else ".idx"

    String output_file = if output_format == "MAF" then output_maf else output_vcf
    String output_file_index = if output_format == "MAF" then output_maf_index else output_vcf_idx

    String transcript_selection_arg = if defined(transcript_selection_list) then " --transcript-list " else ""
    String annotation_def_arg = if defined(annotation_defaults) then " --annotation-default " else ""
    String annotation_over_arg = if defined(annotation_overrides) then " --annotation-override " else ""
    String filter_funcotations_args = if defined(filter_funcotations) && (filter_funcotations) then " --remove-filtered-variants " else ""
    String excluded_fields_args = if defined(funcotator_excluded_fields) then " --exclude-field " else ""

    String interval_list_arg = if defined(interval_list) then " -L " else ""

    String extra_args_arg = select_first([extra_args, ""])

    String requester_pays_arg = if defined(gcs_project_for_requester_pays)  then " --gcs-project-for-requester-pays " else ""

    # ==============
    # Runtime options:
    String gatk_docker

    File? gatk_override
    Int? mem
    Int? preemptible_attempts
    Int? max_retries
    Int? disk_space_gb
    Int? cpu

    Boolean use_ssd = false

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int default_ram_mb = 1024 * 3
    Int machine_mem = if defined(mem) then mem *1024 else default_ram_mb
    Int command_mem = machine_mem - 1024

    # Calculate disk size:
    Float ref_size_gb = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
    Float vcf_size_gb = size(input_vcf, "GiB") + size(input_vcf_idx, "GiB")
    Float ds_size_gb = size(data_sources_tar_gz, "GiB")

    Int default_disk_space_gb = ceil( ref_size_gb + (ds_size_gb * 2) + (vcf_size_gb * 10) ) + 20

    # Silly hack to allow us to use the dollar sign in the command section:
    String dollar = "$"

    command <<<
        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        # =======================================
        # Hack to validate our WDL inputs:
        #
        # NOTE: This happens here so that we don't waste time copying down the data sources if there's an error.

        if [[ "${output_format}" != "MAF" ]] && [[ "${output_format}" != "VCF" ]] ; then
            echo "ERROR: Output format must be MAF or VCF."
        fi

        # =======================================
        # Handle our data sources:

        # Extract the tar.gz:
        echo "Extracting data sources tar/gzip file..."
        mkdir datasources_dir
        tar zxvf ${data_sources_tar_gz} -C datasources_dir --strip-components 1
        DATA_SOURCES_FOLDER="$PWD/datasources_dir"

        # Handle gnomAD:
        if ${use_gnomad} ; then
            echo "Enabling gnomAD..."
            for potential_gnomad_gz in gnomAD_exome.tar.gz gnomAD_genome.tar.gz ; do
                if [[ -f ${dollar}{DATA_SOURCES_FOLDER}/${dollar}{potential_gnomad_gz} ]] ; then
                    cd ${dollar}{DATA_SOURCES_FOLDER}
                    tar -zvxf ${dollar}{potential_gnomad_gz}
                    cd -
                else
                    echo "ERROR: Cannot find gnomAD folder: ${dollar}{potential_gnomad_gz}" 1>&2
                    false
                fi
            done
        fi

        # =======================================
        # Run Funcotator:
        gatk --java-options "-Xmx${command_mem}m" Funcotator \
            --data-sources-path $DATA_SOURCES_FOLDER \
            --ref-version ${reference_version} \
            --output-file-format ${output_format} \
            -R ${ref_fasta} \
            -V ${input_vcf} \
            -O ${output_file} \
            ${interval_list_arg} ${default="" interval_list} \
            --annotation-default normal_barcode:${default="Unknown" control_id} \
            --annotation-default tumor_barcode:${default="Unknown" case_id} \
            --annotation-default Center:${default="Unknown" sequencing_center} \
            --annotation-default source:${default="Unknown" sequence_source} \
            ${"--transcript-selection-mode " + transcript_selection_mode} \
            ${transcript_selection_arg}${default="" sep=" --transcript-list " transcript_selection_list} \
            ${annotation_def_arg}${default="" sep=" --annotation-default " annotation_defaults} \
            ${annotation_over_arg}${default="" sep=" --annotation-override " annotation_overrides} \
            ${excluded_fields_args}${default="" sep=" --exclude-field " funcotator_excluded_fields} \
            ${filter_funcotations_args} \
            ${extra_args_arg} \
            ${requester_pays_arg}${default="" gcs_project_for_requester_pays}

        # =======================================
        # Make sure we have a placeholder index for MAF files so this workflow doesn't fail:
        if [[ "${output_format}" == "MAF" ]] ; then
            touch ${output_maf_index}
        fi
    >>>

    runtime {
        docker: gatk_docker
        bootDiskSizeGb: 20
        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + if use_ssd then " SSD" else " HDD"
        preemptible: select_first([preemptible_attempts, 3])
        maxRetries: select_first([max_retries, 0])
        cpu: select_first([cpu, 1])
    }

    output {
        File funcotated_output_file = "${output_file}"
        File funcotated_output_file_index = "${output_file_index}"
    }
}
