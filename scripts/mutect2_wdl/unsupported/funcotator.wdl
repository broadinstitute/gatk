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
#     funcotator_extra_args          -  Extra command-line arguments to pass through to Funcotator.
#
# This WDL needs to decide whether to use the ``gatk_jar`` or ``gatk_jar_override`` for the jar location.  As of cromwell-0.24,
# this logic *must* go into each task.  Therefore, there is a lot of duplicated code.  This allows users to specify a jar file
# independent of what is in the docker file.  See the README.md for more info.
#
# NOTE: This only does VCF output right now!
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
#    String output_format

    File? interval_list
    File? data_sources_tar_gz
    String? transcript_selection_mode
    Array[String]? transcript_selection_list
    Array[String]? annotation_defaults
    Array[String]? annotation_overrides
    File? gatk4_jar_override

    String? funcotator_extra_args

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
     # inputs
     File ref_fasta
     File ref_fasta_index
     File ref_dict
     File input_vcf
     File input_vcf_idx
     String reference_version
     String output_file_base_name
     String output_format
     Boolean compress
     String output_vcf = output_file_base_name + if compress then ".vcf.gz" else ".vcf"
     String output_vcf_index = output_vcf +  if compress then ".tbi" else ".idx"

     File? data_sources_tar_gz
     String? transcript_selection_mode
     Array[String]? transcript_selection_list
     Array[String]? annotation_defaults
     Array[String]? annotation_overrides
     Boolean filter_funcotations
     File? interval_list

     String? extra_args

     # ==============
     # Process input args:
     String transcript_selection_arg = if defined(transcript_selection_list) then " --transcript-list " else ""
     String annotation_def_arg = if defined(annotation_defaults) then " --annotation-default " else ""
     String annotation_over_arg = if defined(annotation_overrides) then " --annotation-override " else ""
     String filter_funcotations_args = if (filter_funcotations) then " --remove-filtered-variants " else ""
     String interval_list_arg = if defined(interval_list) then " -L " else ""
     String extra_args_arg = select_first([extra_args, ""])
     # ==============

     # runtime

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
             -O ${output_vcf} \
             ${interval_list_arg} ${default="" interval_list} \
             ${"--transcript-selection-mode " + transcript_selection_mode} \
             ${transcript_selection_arg}${default="" sep=" --transcript-list " transcript_selection_list} \
             ${annotation_def_arg}${default="" sep=" --annotation-default " annotation_defaults} \
             ${annotation_over_arg}${default="" sep=" --annotation-override " annotation_overrides} \
             ${filter_funcotations_args} \
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
         File funcotated_vcf = "${output_vcf}"
         File funcotated_vcf_index = "${output_vcf_index}"
     }
 }