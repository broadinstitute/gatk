# Compute the difference between the elapsed times of the truth and the call timing values.
#
# Description of inputs:
#
#   Required:
#     gatk_docker                    - GATK Docker image in which to run
#
#     call_timing_file               - The produced data set.  Variant Context File (VCF) containing the variants.
#     truth_timing_file              - The truth data set.  Variant Context File (VCF) containing the variants.
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

workflow CompareTiming {

    # ------------------------------------------------
    # Input args:
    String gatk_docker = "broadinstitute/gatk:gatkbase-2.0.2"

    File call_timing_file
    File truth_timing_file

    String? base_timing_output_name

    Int?  mem_gb
    Int? preemptible_attempts
    Int? disk_space_gb
    Int? cpu
    Int? boot_disk_size_gb

    # ------------------------------------------------
    # Call our tasks:
    call CompareTimingTask {
        input:
            truth_timing_file         = truth_timing_file,
            call_timing_file          = call_timing_file,

            base_timing_output_name   = base_timing_output_name,

            gatk_docker               = gatk_docker,
            mem                       = mem_gb,
            preemptible_attempts      = preemptible_attempts,
            disk_space_gb             = disk_space_gb,
            cpu                       = cpu,
            boot_disk_size_gb         = boot_disk_size_gb
    }

    # ------------------------------------------------
    # Outputs:
    output {
        File summary_metrics = CompareTimingTask.timing_diff
    }

}

# ==========================================================================================
# ==========================================================================================
# ==========================================================================================

task CompareTimingTask {

    ####################################################################################
    # Inputs:
    File truth_timing_file
    File call_timing_file

    String? base_timing_output_name
    String? bam_name

    ####################################################################################
    # Runtime Inputs:
    String gatk_docker

    Int? mem
    Int? preemptible_attempts
    Int? disk_space_gb
    Int? cpu
    Int? boot_disk_size_gb

    ####################################################################################
    # Default values:
    String timing_diff_file_name = if defined(base_timing_output_name) then basename(base_timing_output_name) + ".timingDiff.txt" else "timingDiff.txt"
    String timing_csv_file_name = if defined(base_timing_output_name) then basename(base_timing_output_name) + ".timingDiff.csv" else "timingDiff.csv"

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
        truthElapsed=$( grep "Elapsed" ${truth_timing_file} | sed 's#.*[ \t]##')
        callElapsed=$( grep "Elapsed" ${call_timing_file} | sed 's#.*[ \t]##')

        echo "Control,${bam_name},$truthElapsed" >> ${timing_csv_file_name}
        echo "Test,${bam_name},$truthElapsed" >> ${timing_csv_file_name}

        echo "truthElapsed = $truthElapsed"
        echo "callElapsed = $callElapsed"

        echo "truthElapsed = $truthElapsed" >> ${timing_diff_file_name}
        echo "callElapsed = $callElapsed" >> ${timing_diff_file_name}

        timeDiff=$( python -c "print($callElapsed - $truthElapsed)" )
        timeRatio=$( python -c "print($callElapsed/$truthElapsed)" )

        echo "# timeDiff = callElapsed - truthElapsed"
        echo "timeDiff = $timeDiff"

        echo "# timeRatio = callElapsed/truthElapsed"
        echo "timeRatio = $timeRatio"

        echo "# timeDiff = callElapsed - truthElapsed" >> ${timing_diff_file_name}
        echo "timeDiff: $timeDiff" >> ${timing_diff_file_name}

        echo "# timeRatio = callElapsed/truthElapsed" >> ${timing_diff_file_name}
        echo "timeRatio: $timeRatio" >> ${timing_diff_file_name}
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
        File timing_csv          = timing_csv_file_name
        File timing_diff         = timing_diff_file_name
        String run_title         = call_timing_file
    }
}
