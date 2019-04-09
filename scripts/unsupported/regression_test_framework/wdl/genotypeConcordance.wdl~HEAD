# Define the GenotypeConcordance task, to run GenotypeConcordance on a VCF file.
#
# Description of inputs:
#
#   Required:
#     gatk_docker                    - GATK Docker image in which to run
#
#     call_vcf                       - The produced data set.  Variant Context File (VCF) containing the variants.
#     call_index                     - The produced data set index.  Variant Context File (VCF) index for the given actual_vcf.
#     call_sample                    - Sample in the VCF to compare.
#
#     truth_vcf                      - The truth data set.  Variant Context File (VCF) containing the variants.
#     truth_index                    - The truth data set index.  Variant Context File (VCF) index for the given expected_vcf.
#     truth_sample                   - Sample in the VCF to use as truth data.
#
#     output_base_name               - Base file name for generated output files.
#
#   Optional:
#     File? interval_list            -  Interval list over which to call variants on the given BAM file.
#
#     gatk4_jar_override             - Override Jar file containing GATK 4.0.  Use this when overriding the docker JAR or when using a backend without docker.
#     mem_gb                         - Amount of memory to give the runtime environment.
#     disk_space_gb                  - Amount of disk space to give the runtime environment.
#     cpu                            - Number of CPUs to give the runtime environment.
#     boot_disk_size_gb              - Amount of boot disk space to give the runtime environment.
#     preemptible_attempts           - Number of times the comparison can be preempted.
#
# This WDL needs to decide whether to use the ``gatk_jar`` or ``gatk_jar_override`` for the jar location.  As of cromwell-0.24,
# this logic *must* go into each task.  Therefore, there is a lot of duplicated code.  This allows users to specify a jar file
# independent of what is in the docker file.  See the README.md for more info.

workflow GenotypeConcordance {

    # ------------------------------------------------
    # Input args:
    String gatk_docker

    File call_vcf
    File call_index
    String call_sample

    File truth_vcf
    File truth_index
    String truth_sample

    String output_base_name

    File? interval_list

    File? gatk4_jar_override
    Int?  mem_gb
    Int? preemptible_attempts
    Int? disk_space_gb
    Int? cpu
    Int? boot_disk_size_gb

    # ------------------------------------------------
    # Call our tasks:
    call GenotypeConcordanceTask {
        input:
            call_vcf                  = call_vcf,
            call_index                = call_index,
            call_sample               = call_sample,

            truth_vcf                 = truth_vcf,
            truth_index               = truth_index,
            truth_sample              = truth_sample,

            interval_list             = interval_list,

            output_base_name          = output_base_name,

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
        File summary_metrics                = GenotypeConcordanceTask.summary_metrics
        File detail_metrics                 = GenotypeConcordanceTask.detail_metrics
        File contingency_metrics            = GenotypeConcordanceTask.contingency_metrics
        File timing_info                    = GenotypeConcordanceTask.timing_info
    }

}

# ==========================================================================================
# ==========================================================================================
# ==========================================================================================

task GenotypeConcordanceTask {

    ####################################################################################
    # Inputs:
    File call_vcf
    File call_index
    String call_sample

    File truth_vcf
    File truth_index
    String truth_sample

    File? interval_list

    String output_base_name

    ####################################################################################
    # Runtime Inputs:
    String gatk_docker

    File? gatk_override
    Int? mem
    Int? preemptible_attempts
    Int? disk_space_gb
    Int? cpu
    Int? boot_disk_size_gb

    # ------------------------------------------------
    # Process input args:
    String interval_list_arg = if defined(interval_list) then " --INTERVALS " else ""

    String timing_output_file = "GenotypeConcordanceTask" + ".timingInformation.txt"

    ####################################################################################
    # Define default values and set up values for running:
    Boolean use_ssd = false

    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 3 * 1024
    # WARNING: In the workflow, you should calculate the disk space as an input to this task (disk_space_gb).  Please see [TODO: Link from Jose] for examples.
    Int default_disk_space_gb = 100

    Int default_boot_disk_size_gb = 15

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1024 else default_ram_mb
    Int command_mem = machine_mem - 1024

    ####################################################################################
    # Do the work:
    command {
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        startTime=`date +%s.%N`
        echo "StartTime: $startTime" > ${timing_output_file}

        gatk --java-options "-Xmx${command_mem}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
            GenotypeConcordance \
                --CALL_VCF ${call_vcf} \
                --CALL_SAMPLE ${call_sample} \
                --TRUTH_VCF ${truth_vcf} \
                --TRUTH_SAMPLE ${truth_sample} \
                -O ${output_base_name} \
                ${interval_list_arg}${default="" sep=" --INTERVALS " interval_list} \
                --INTERSECT_INTERVALS true \

        endTime=`date +%s.%N`
        echo "EndTime: $endTime" >> ${timing_output_file}
        elapsedTime=`python -c "print( $endTime - $startTime )"`
        echo "Elapsed Time: $elapsedTime" >> ${timing_output_file}
    }

    ####################################################################################
    # Runtime Environment:
    runtime {
        cpu: select_first([cpu, 1])
        memory: machine_mem + " MB"
        bootDiskSizeGb: select_first([disk_space_gb, default_disk_space_gb])
        disks: "local-disk " + select_first([boot_disk_size_gb, default_boot_disk_size_gb]) + if use_ssd then " SSD" else " HDD"
        docker: "${gatk_docker}"
        preemptible: select_first([preemptible_attempts, 0])
    }

    ####################################################################################
    # Outputs:
    output {
        File summary_metrics     = "${output_base_name}.genotype_concordance_summary_metrics"
        File detail_metrics      = "${output_base_name}.genotype_concordance_detail_metrics"
        File contingency_metrics = "${output_base_name}.genotype_concordance_contingency_metrics"
        File timing_info         = timing_output_file
    }
}
