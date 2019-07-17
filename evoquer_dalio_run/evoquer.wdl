# Run Evoquer on a list of intervals.
#
# Description of inputs:
#
#  Required:
#    String gatk_docker                 - GATK Docker image in which to run
#
#    File reference                     - Reference fasta
#
#    Array[String] intervals            - Variant Context File (VCF) containing the variants to annotate.
#
#    String output_file                 - Base name of desired output file WITHOUT extension.
#
#  Optional:
#
#     File gatk4_jar_override           -  Override Jar file containing GATK 4.  Use this when overriding the docker JAR or when using a backend without docker.
#     Int  mem_gb                       -  Amount of memory to give to the machine running each task in this workflow (in gb).
#     Int  preemptible_attempts         -  Number of times to allow each task in this workflow to be preempted.
#     Int  disk_space_gb                -  Amount of storage disk space (in Gb) to give to each machine running each task in this workflow.
#     Int  cpu                          -  Number of CPU cores to give to each machine running each task in this workflow.
#     Int  boot_disk_size_gb            -  Amount of boot disk space (in Gb) to give to each machine running each task in this workflow.
#
# This WDL needs to decide whether to use the ``gatk_jar`` or ``gatk_jar_override`` for the jar location.  As of cromwell-0.24,
# this logic *must* go into each task.  Therefore, there is a lot of duplicated code.  This allows users to specify a jar file
# independent of what is in the docker file.  See the README.md for more info.
#
workflow Evoquer {
    String gatk_docker

    File reference
    File reference_index
    File reference_dict

    String project_id

    File dataset_map

    File interval_file

    String output_file_base_name

    String batch_mode

    String disable_gnarly_genotyper

    File? gatk4_jar_override
    Int?  mem_gb
    Int?  preemptible_attempts
    Int?  disk_space_gb
    Int?  cpu
    Int?  boot_disk_size_gb

    Array[String] intervals = read_lines(interval_file)

    scatter(interval in intervals) {
        call EvoquerTask {
            input:
                reference             = reference,
                reference_index       = reference_index,
                reference_dict        = reference_dict,
                interval              = interval,
                output_file           = "${output_file_base_name}_${interval}.vcf.gz",
                project_id            = project_id,
                dataset_map           = dataset_map,
                batch_mode            = batch_mode,
                disable_gnarly_genotyper = disable_gnarly_genotyper,
                gatk_docker           = gatk_docker,
                gatk_override         = gatk4_jar_override,
                mem                   = mem_gb,
                preemptible_attempts  = preemptible_attempts,
                disk_space_gb         = disk_space_gb,
                cpu                   = cpu,
                boot_disk_size_gb     = boot_disk_size_gb
        }
    }

    output {
        Array[File] genotyped_variant_calls_files = EvoquerTask.evoked_variants
        Array[File] timing_output_file = EvoquerTask.timing_output
    }
}

################################################################################

task EvoquerTask {

    # ------------------------------------------------
    # Input args:
    File reference
    File reference_index
    File reference_dict

    String interval
    String output_file

    String project_id

    File dataset_map

    String batch_mode

    String disable_gnarly_genotyper

    # Runtime Options:
    String gatk_docker
    File? gatk_override
    Int? mem
    Int? preemptible_attempts
    Int? disk_space_gb
    Int? cpu
    Int? boot_disk_size_gb

    # ------------------------------------------------
    # Misc settings:
    String dollar = "$"

    # ------------------------------------------------
    # Process input args:

    # Timing info:
    String timing_output_file = "Evoquer.timingInformation.txt"

    # ------------------------------------------------
    # Get machine settings:
    Boolean use_ssd = false

    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 1024 * 3
    # WARNING: In the workflow, you should calculate the disk space as an input to this task (disk_space_gb).  Please see [TODO: Link from Jose] for examples.
    Int default_disk_space_gb = 100

    Int default_boot_disk_size_gb = 15

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1024 else default_ram_mb
    Int command_mem = machine_mem - 1024

    # ------------------------------------------------
    # Run our command:
    command <<<
        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        startTime=`date +%s.%N`
        echo "StartTime: $startTime" > ${timing_output_file}

        gatk --java-options "-Xmx${command_mem}m" \
            Evoquer \
                -R "${reference}" \
                -O "${output_file}" \
                --project-id "${project_id}" \
                --dataset-map "${dataset_map}" \
                -L "${interval}" \
                --run-query-in-batch-mode "${batch_mode}" \
                --disable-gnarly-genotyper "${disable_gnarly_genotyper}"

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
        preemptible: if defined(preemptible_attempts) then preemptible_attempts else 0
        cpu: select_first([cpu, 1])
    }

    # ------------------------------------------------
    # Outputs:
    output {
        File evoked_variants = output_file
        File timing_output = timing_output_file
    }
 }
