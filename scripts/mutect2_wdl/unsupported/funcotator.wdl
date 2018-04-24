# Run Funcotator on a set of called variants from Mutect 2.
#
# Description of inputs:
#
#   Required:
#     gatk_docker                -  GATK Docker image in which to run
#     ref_fasta                  -  Reference FASTA file.
#     ref_fasta_index            -  Reference FASTA file index.
#     ref_fasta_dict             -  Reference FASTA file sequence dictionary.
#     variant_vcf_to_funcotate   -  Variant Context File (VCF) containing the variants to annotate.
#     reference_version          -  Version of the reference being used.  Either `hg19` or `hg38`.
#     output_file_name           -  Path to desired output file.
#     output_format              -  Format of the output file (`VCF` or `MAF`).
#
#   Optional:
#     data_sources_tar_gz        -  Path to tar.gz containing the data sources for Funcotator to create annotations.
#     transcript_selection_mode  -  Method of detailed transcript selection.  This will select the transcript for detailed annotation (either `CANONICAL` or `BEST_EFFECT`).
#     transcript_selection_list  -  Set of transcript IDs to use for annotation to override selected transcript.
#     annotation_defaults        -  Annotations to include in all annotated variants if the annotation is not specified in the data sources (in the format <ANNOTATION>:<VALUE>).  This will add the specified annotation to every annotated variant if it is not already present.
#     annotation_overrides       -  Override values for annotations (in the format <ANNOTATION>:<VALUE>).  Replaces existing annotations of the given name with given values.
#     gatk4_jar_override         -  Override Jar file containing GATK 4.0.  Use this when overriding the docker JAR or when using a backend without docker.
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
    File variant_vcf_to_annotate
    String reference_version
    String output_file_name
    String output_format

    String? data_sources_tar_gz
    String? transcript_selection_mode
    Array[String]? transcript_selection_list
    Array[String]? annotation_defaults
    Array[String]? annotation_overrides
    File? gatk4_jar_override

    call MakeItFunky {
        input:
            ref_fasta                 = ref_fasta,
            ref_fasta_index           = ref_fasta_index,
            ref_dict                  = ref_dict,
            variant_vcf_to_annotate   = variant_vcf_to_annotate,
            reference_version         = reference_version,
            output_file_name          = output_file_name,
            output_format             = output_format,
            data_sources_tar_gz       = data_sources_tar_gz,
            transcript_selection_mode = transcript_selection_mode,
            transcript_selection_list = transcript_selection_list,
            annotation_defaults       = annotation_defaults,
            annotation_overrides      = annotation_overrides,
            gatk4_jar_override        = gatk4_jar_override,
            gatk_docker               = gatk_docker
    }

    output {
        File annotated_out = MakeItFunky.annotated_out
        File annotated_out_idx = MakeItFunky.annotated_out_idx
    }
}


task MakeItFunky {
    # Inputs for this task
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File variant_vcf_to_annotate
    String reference_version
    String output_file_name
    String output_format

    File? data_sources_tar_gz
    String? transcript_selection_mode
    Array[String]? transcript_selection_list
    Array[String]? annotation_defaults
    Array[String]? annotation_overrides

    # ==============
    # Process input args:
    String transcript_selection_arg = if defined(transcript_selection_list) then " --transcript-list " else ""
    String annotation_def_arg = if defined(annotation_defaults) then " --annotation-default " else ""
    String annotation_over_arg = if defined(annotation_overrides) then " --annotation-override " else ""

    # ==============
    # Runtime parameters

    String gatk_docker

    File? gatk4_jar_override
    Int? mem_gb
    Int? preemptible_attempts
    Int? disk_space_gb
    Int? cpu

    Boolean use_ssd = false

    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 3000
    # WARNING: In the workflow, you should calculate the disk space as an input to this task (disk_space_gb).  Please see [TODO: Link from Jose] for examples.
    Int default_disk_space_gb = 100

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem_gb) then mem_gb *1000 else default_ram_mb
    Int command_mem = machine_mem - 1000

command <<<
    set -e
    export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk4_jar_override}

    DATA_SOURCES_TAR_GZ=${data_sources_tar_gz}
    if [[ ! -e $DATA_SOURCES_TAR_GZ ]] ; then
        # We have to download the data sources:
        echo "Data sources gzip does not exist: $DATA_SOURCES_TAR_GZ"
        echo "Downloading default data sources..."
        wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/funcotator/funcotator_dataSources.v1.0.20180105.tar.gz
        tar -zxf funcotator_dataSources.v1.0.20180105.tar.gz
        DATA_SOURCES_FOLDER=funcotator_dataSources.v1.0.20180105
    else
        # Extract the tar.gz:
        mkdir datasources_dir
        tar zxvf ${data_sources_tar_gz} -C datasources_dir --strip-components 1
        DATA_SOURCES_FOLDER="$PWD/datasources_dir"
    fi

    /gatk/gatk --java-options "-Xmx${command_mem}m" \
        Funcotator --data-sources-path $DATA_SOURCES_FOLDER \
            --ref-version ${reference_version} \
            -R ${ref_fasta} \
            -V ${variant_vcf_to_annotate} \
            -O ${output_file_name} \
            --output-file-format ${output_format}
            ${"--transcript-selection-mode " + transcript_selection_mode} \
            ${transcript_selection_arg}${default="" sep=" --transcript-list " transcript_selection_list} \
            ${annotation_def_arg}${default="" sep=" --annotation-default " annotation_defaults} \
            ${annotation_over_arg}${default="" sep=" --annotation-override " annotation_overrides}
>>>

  runtime {
    docker: gatk_docker
    memory: machine_mem + " MB"
    # Note that the space before SSD and HDD should be included.
    disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + if use_ssd then " SSD" else " HDD"
    preemptible: select_first([preemptible_attempts, 3])
    cpu: select_first([cpu, 1])
  }

  output {
    File annotated_out = "${output_file_name}"
    File annotated_out_idx = "${output_file_name}.idx"
  }
}