# Calculate read counts of different types.
#
# Description of inputs:
#
#   Required:
#     String docker_image                         - Docker image in which to run the analysis.
#     Array[String] bam_files                     - Array of paths to bam files to read.  Assumes these are GS paths.
#
#   Optional:
#     String extra_args                           - Extra arguments to pass through to the tool.
#     Int  mem                                    - Amount of memory to give to the machine running each task in this workflow.
#     Int  preemptible_attempts                   - Number of times to allow each task in this workflow to be preempted.
#     Int  disk_space_gb                          - Amount of storage disk space (in Gb) to give to each machine running each task in this workflow.
#     Int  cpu                                    - Number of CPU cores to give to each machine running each task in this workflow.
#     Int  boot_disk_size_gb                      - Amount of boot disk space (in Gb) to give to each machine running each task in this workflow.
#
workflow LRMAReadCounter {
    String docker_image
    Array[String] bam_files

    Int? mem
    Int? preemptible_attempts
    Int? disk_space_gb
    Int? cpu
    Int? boot_disk_size_gb

    String? extra_args

    scatter (bam_file in bam_files) {
        call LRMAReadCounterTask {
            input:
                docker               = docker_image,
                bam_file_name        = bam_file,
                extra_args           = extra_args,
                mem                  = mem,
                preemptible_attempts = preemptible_attempts,
                disk_space_gb        = disk_space_gb,
                cpu                  = cpu,
                boot_disk_size_gb    = boot_disk_size_gb
        }
    }

    output {
        Array[File] reads_count_files = LRMAReadCounterTask.reads_count_file
        Array[File] timing_files = LRMAReadCounterTask.timing_file
    }
}

################################################################################

task LRMAReadCounterTask {

    # ------------------------------------------------
    # Input args:

    String docker
    String bam_file_name

    String? extra_args

    Int? mem
    Int? preemptible_attempts
    Int? disk_space_gb
    Int? cpu
    Int? boot_disk_size_gb

    # ------------------------------------------------
    # Process input args:

    String output_csv_file = sub(sub(bam_file_name, '.bam$', ''), '^.*/', '') + '_read_stats.csv'
    String timing_output_file = "timingInformation.txt"

    # ------------------------------------------------
    # Get machine settings:
    Boolean use_ssd = false

    Int default_ram_mb = 4 * 1024
    Float reads_size_gb = 1024
    Int default_disk_space_gb = ceil((reads_size_gb * 2) + 20)

    Int default_boot_disk_size_gb = 15

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1024 else default_ram_mb
    Int command_mem = machine_mem - 1024

    # ------------------------------------------------
    # Run our command:
     command <<<
       set -e

       startTime=`date +%s.%N`
       echo "StartTime: $startTime" > ${timing_output_file}

        gatk --java-options "-Xmx${command_mem}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
            LRMAReadCounter \
            -I ${bam_file_name} \
            --use-file-as-row-header \
            --output-csv-file ${output_csv_file} \
            ${extra_args} \


        endTime=`date +%s.%N`
        echo "EndTime: $endTime" >> ${timing_output_file}
        elapsedTime=`python -c "print( $endTime - $startTime )"`
        echo "Elapsed Time: $elapsedTime" >> ${timing_output_file}
     >>>

    # ------------------------------------------------
    # Runtime settings:
     runtime {
         docker: docker
         memory: machine_mem + " MB"
         disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + if use_ssd then " SSD" else " HDD"
         bootDiskSizeGb: select_first([boot_disk_size_gb, default_boot_disk_size_gb])
         preemptible: 0
         cpu: select_first([cpu, 1])
     }

    # ------------------------------------------------
    # Outputs:
     output {
        File reads_count_file = output_csv_file
        File timing_file = timing_output_file
     }
 }

