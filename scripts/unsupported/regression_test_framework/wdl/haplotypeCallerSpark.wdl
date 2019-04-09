# Calls variants on the given germline BAM file using HaplotypeCallerSpark.
#
# Description of inputs:
#
#   Required:
#     String gatk_docker                                -  GATK Docker image in which to run
#
#     File input_bam                                    -  Input reads over which to call small variants with Haplotype Caller.
#     File input_bam_index                              -  Index for the input BAM file.
#     File ref_fasta                                    -  Reference FASTA file.
#     File ref_fasta_dict                               -  Reference FASTA file dictionary.
#     File ref_fasta_index                              -  Reference FASTA file index.
#
#     String spark_runner                               -  Which type of spark cluster to run on.
#     String cluster                                    -  Server head of the cluster on which to run.
#     String project                                    -  Project on the server under which to run.
#
#   Optional:
#
#     File? interval_list                               -  Interval list over which to call variants on the given BAM file.
#     Boolean? gvcf_mode                                -  Whether to run in GVCF mode (default: false).
#     Float? contamination                              -  Contamination threshold for HaplotypeCaller (default: 0.0).
#
#     Int? num_executors                                -  Number of machines on which to execute this job.
#     Int? executor_cores                               -  Number of cores per machine.
#     Int? executor_memory_gb                           -  Amount of memory to give to each machine (in gb).
#     String? conf_options                              -  Miscellaneous spark configuration options.
#
#     File gatk4_jar_override                           -  Override Jar file containing GATK 4.  Use this when overriding the docker JAR or when using a backend without docker.
#     Int  mem_gb                                       -  Amount of memory to give to the machine running each task in this workflow (in gb).
#     Int  preemptible_attempts                         -  Number of times to allow each task in this workflow to be preempted.
#     Int  disk_space_gb                                -  Amount of storage disk space (in Gb) to give to each machine running each task in this workflow.
#     Int  cpu                                          -  Number of CPU cores to give to each machine running each task in this workflow.
#     Int  boot_disk_size_gb                            -  Amount of boot disk space (in Gb) to give to each machine running each task in this workflow.
#
# This WDL needs to decide whether to use the ``gatk_jar`` or ``gatk_jar_override`` for the jar location.  As of cromwell-0.24,
# this logic *must* go into each task.  Therefore, there is a lot of duplicated code.  This allows users to specify a jar file
# independent of what is in the docker file.  See the README.md for more info.
workflow HaplotypeCallerSpark {

    # ------------------------------------------------
    # Input args:
    String gatk_docker

    File input_bam
    File input_bam_index
    File ref_fasta
    File ref_fasta_dict
    File ref_fasta_index

    String out_vcf_name

    File? interval_list
    Boolean? gvcf_mode
    Float? contamination
    Int? interval_padding

    String spark_runner
    String cluster
    String project
    Int? num_executors
    Int? executor_cores
    Int? executor_memory_gb
    String? conf_options

    File? gatk4_jar_override
    Int?  mem_gb
    Int? preemptible_attempts
    Int? disk_space_gb
    Int? cpu
    Int? boot_disk_size_gb

    # ------------------------------------------------
    # Call our tasks:
    call HaplotypeCallerSparkTask {
        input:
            input_bam                 = input_bam,
            input_bam_index           = input_bam_index,
            ref_fasta                 = ref_fasta,
            ref_fasta_dict            = ref_fasta_dict,
            ref_fasta_index           = ref_fasta_index,

            interval_list             = interval_list,
            gvcf_mode                 = gvcf_mode,
            contamination             = contamination,
            interval_padding          = interval_padding,

            out_file_name             = out_vcf_name,

            spark_runner              = spark_runner,
            cluster                   = cluster,
            project                   = project,
            num_executors             = num_executors,
            executor_cores            = executor_cores,
            executor_memory_gb        = executor_memory_gb,
            conf_options              = conf_options,

            gatk_docker               = gatk_docker,
            gatk_override             = gatk4_jar_override,
            mem                       = mem_gb,
            preemptible_attempts      = preemptible_attempts,
            disk_space_gb             = disk_space_gb,
            cpu                       = cpu,
            boot_disk_size_gb         = boot_disk_size_gb
    }

    # ------------------------------------------------
    # Outputs:
    output {
        File vcf_out     = HaplotypeCallerSparkTask.output_vcf
        File vcf_out_idx = HaplotypeCallerSparkTask.output_vcf_index
        File timingInfo  = HaplotypeCallerSparkTask.timing_info
    }
}

# ==========================================================================================
# ==========================================================================================
# ==========================================================================================

task HaplotypeCallerSparkTask {

    # ------------------------------------------------
    # Input args:

    # Required:
    File input_bam
    File input_bam_index
    File ref_fasta
    File ref_fasta_dict
    File ref_fasta_index

    File? interval_list
    Boolean? gvcf_mode
    Float? contamination
    Int? interval_padding

    # Output Names:
    String out_file_name

    # Spark Options:
    String spark_runner
    String cluster
    String project
    Int? num_executors
    Int? executor_cores
    Int? executor_memory_gb
    String? conf_options

    # Runtime Options:
    String gatk_docker
    File? gatk_override
    Int? mem
    Int? preemptible_attempts
    Int? disk_space_gb
    Int? cpu
    Int? boot_disk_size_gb

    # ------------------------------------------------
    # Process input args:
    String interval_list_arg = if defined(interval_list) then " -L " else ""
    String contamination_arg = if defined(contamination) then " --contamination " else ""
    String interval_padding_arg = if defined(interval_padding) then " --interval-padding " else ""

    String index_format = if sub(out_file_name, ".*\\.", "") == "vcf" then "idx" else "tbi"

    String timing_output_file = out_file_name + ".timingInformation.txt"

    # ------------------------------------------------
    # Get machine settings:
    Boolean use_ssd = false

    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 3 * 1024
    # WARNING: In the workflow, you should calculate the disk space as an input to this task (disk_space_gb).  Please see [TODO: Link from Jose] for examples.
    Int default_disk_space_gb = 100

    Int default_boot_disk_size_gb = 15

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1024 else default_ram_mb
    Int command_mem = machine_mem - 1024

    # ------------------------------------------------
    # Process optional spark args:
    Int num_executors_default = 5
    Int executor_cores_default = 4
    String conf_options_default = ""

    Int num_executors_val = if defined(num_executors) then num_executors * 1 else num_executors_default
    Int executor_cores_val = if defined(executor_cores) then executor_cores * 1 else executor_cores_default
    Int executor_memory_gb_val = if defined(executor_memory_gb) then executor_memory_gb * 1 else machine_mem

    String conf_options_arg = if defined(conf_options) then " --conf " else ""

    # ------------------------------------------------
    # Run our command:
    command <<<
        set -e

        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        startTime=`date +%s.%N`
        echo "StartTime: $startTime" > ${timing_output_file}

        gatk --java-options "-Xmx${command_mem}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
            HaplotypeCallerSpark \
                -I ${input_bam} \
                -O ${out_file_name} \
                -R ${ref_fasta} \
                --read-validation-stringency STRICT \
                ${interval_list_arg}${default="" sep=" -L " interval_list} \
                ${contamination_arg}${default="" sep=" --contamination " contamination} \
                ${interval_padding_arg}${default="" sep=" --interval-padding " interval_padding} \
                ${true="-ERC GVCF" false="" gvcf_mode} \
                -- \
                --spark-runner ${spark_runner} \
                --cluster ${cluster} \
                --project ${project} \
                --num-executors ${num_executors_val} \
                --executor-cores ${executor_cores_val} \
                --executor-memory ${executor_memory_gb_val}g \
                ${conf_options_arg}${default="" sep=" --conf " conf_options} \

        endTime=`date +%s.%N`
        echo "EndTime: $endTime" >> ${timing_output_file}
        elapsedTime=`python -c "print( $endTime - $startTime )"`
        echo "Elapsed Time: $elapsedTime" >> ${timing_output_file}
    >>>

    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: gatk_docker
        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + if use_ssd then " SSD" else " HDD"
        bootDiskSizeGb: select_first([boot_disk_size_gb, default_boot_disk_size_gb])
        preemptible: select_first([preemptible_attempts, 0])
        cpu: select_first([cpu, 1])
    }

    # ------------------------------------------------
    # Outputs:
    output {
        File output_vcf       = out_file_name
        File output_vcf_index = "${out_file_name}.${index_format}"
        File timing_info      = timing_output_file
    }
}
