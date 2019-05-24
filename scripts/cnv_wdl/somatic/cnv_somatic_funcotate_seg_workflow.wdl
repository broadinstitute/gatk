workflow CNVFuncotateSegmentsWorkflow {

    File input_seg_file
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    String funcotator_ref_version
    File? gatk4_jar_override
    File? funcotator_data_sources_tar_gz
    String? transcript_selection_mode
    File? transcript_selection_list
    Array[String]? annotation_defaults
    Array[String]? annotation_overrides
    Array[String]? funcotator_excluded_fields
    File? interval_list
    String? extra_args

    # Set to true when running local or on-prem
    Boolean? is_removing_untared_datasources

    # Runtime parameters
    String gatk_docker
    Int? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? cpu
    Int? preemptible_attempts

    call FuncotateSegments {
        input:
            input_seg_file = input_seg_file,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            ref_fasta_dict = ref_fasta_dict,
            funcotator_ref_version = funcotator_ref_version,
            gatk4_jar_override = gatk4_jar_override,
            funcotator_data_sources_tar_gz = funcotator_data_sources_tar_gz,
            transcript_selection_mode = transcript_selection_mode,
            transcript_selection_list = transcript_selection_list,
            annotation_defaults = annotation_defaults,
            annotation_overrides = annotation_overrides,
            funcotator_excluded_fields = funcotator_excluded_fields,
            interval_list = interval_list,
            extra_args = extra_args,
            is_removing_untared_datasources = is_removing_untared_datasources,
            gatk_docker = gatk_docker,
            mem_gb = mem_gb,
            disk_space_gb = disk_space_gb,
            use_ssd = use_ssd,
            cpu = cpu,
            preemptible_attempts = preemptible_attempts
    }

    output {
        File funcotated_seg_simple_tsv = FuncotateSegments.funcotated_seg_simple_tsv
        File funcotated_gene_list_tsv = FuncotateSegments.funcotated_gene_list_tsv
    }
}

task FuncotateSegments {

    File input_seg_file
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    String funcotator_ref_version
    File? gatk4_jar_override
    File? funcotator_data_sources_tar_gz = "gs://broad-public-datasets/funcotator/funcotator_dataSources.v1.6.20190124s.tar.gz"
    String? transcript_selection_mode = "CANONICAL"
    File? transcript_selection_list
    Array[String]? annotation_defaults
    Array[String]? annotation_overrides
    Array[String]? funcotator_excluded_fields
    File? interval_list
    String? extra_args

    # Set to true when running local or on-prem
    Boolean? is_removing_untared_datasources

    # Runtime parameters
    String gatk_docker
    Int? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? cpu
    Int? preemptible_attempts

    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 3000
    # WARNING: In the workflow, you should calculate the disk space as an input to this task (disk_space_gb).
    Int default_disk_space_gb = 100

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem_mb = if defined(mem_gb) then mem_gb *1000 else default_ram_mb
    Int command_mem_mb = machine_mem_mb - 1000

    ## Process input args
    String transcript_selection_arg = if defined(transcript_selection_list) then " --transcript-list " else ""
    String annotation_def_arg = if defined(annotation_defaults) then " --annotation-default " else ""
    String annotation_over_arg = if defined(annotation_overrides) then " --annotation-override " else ""
    String excluded_fields_args = if defined(funcotator_excluded_fields) then " --exclude-field " else ""
    String interval_list_arg = if defined(interval_list) then " -L " else ""
    String extra_args_arg = select_first([extra_args, ""])
    Boolean is_removing_untared_datasources_final = select_first([is_removing_untared_datasources, true])
    String removing_untared_datasources = if is_removing_untared_datasources_final then " rm -Rf $DATA_SOURCES_FOLDER " else " echo Not bothering to remove datasources"
    String basename_input_seg_file = basename(input_seg_file)

    command <<<
        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk4_jar_override}

         # Extract our data sources:
         echo "Extracting data sources zip file..."
         mkdir datasources_dir
         tar zxvf ${funcotator_data_sources_tar_gz} -C datasources_dir --strip-components 1
         DATA_SOURCES_FOLDER="$PWD/datasources_dir"

         # Run FuncotateSegments:
         gatk --java-options "-Xmx${command_mem_mb}m" FuncotateSegments \
             --data-sources-path $DATA_SOURCES_FOLDER \
             --ref-version ${funcotator_ref_version} \
             --output-file-format SEG \
             -R ${ref_fasta} \
              --segments ${input_seg_file} \
             -O ${basename_input_seg_file}.funcotated.tsv \
             ${interval_list_arg} ${default="" interval_list} \
             ${"--transcript-selection-mode " + transcript_selection_mode} \
             ${transcript_selection_arg}${default="" sep=" --transcript-list " transcript_selection_list} \
             ${annotation_def_arg}${default="" sep=" --annotation-default " annotation_defaults} \
             ${annotation_over_arg}${default="" sep=" --annotation-override " annotation_overrides} \
             ${excluded_fields_args}${default="" sep=" --exclude-field " funcotator_excluded_fields} \
             ${extra_args_arg}


         echo "Removing $DATA_SOURCES_FOLDER"
         ${removing_untared_datasources}
    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File funcotated_seg_simple_tsv = "${basename_input_seg_file}.funcotated.tsv"
        File funcotated_gene_list_tsv = "${basename_input_seg_file}.funcotated.tsv.gene_list.txt"
    }
}