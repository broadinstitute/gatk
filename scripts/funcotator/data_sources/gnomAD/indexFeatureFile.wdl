# Run IndexFeatureFile on an array of VCF files.
#
# Description of inputs:
#
#   Required:
#     String gatk_docker           -  GATK Docker image in which to run
#     Array[File] variant_vcfs     -  Array of Variant Context Files (VCFs) containing the variants to index.
#
#   Optional:
#     File gatk4_jar_override      -  Override Jar file containing GATK 4.  Use this when overriding the docker JAR or when using a backend without docker.
#     Int  mem                     -  Amount of memory to give to the machine running each task in this workflow.
#     Int  preemptible_attempts    -  Number of times to allow each task in this workflow to be preempted.
#     Int  disk_space_gb           -  Amount of storage disk space (in Gb) to give to each machine running each task in this workflow.
#     Int  cpu                     -  Number of CPU cores to give to each machine running each task in this workflow.
#     Int  boot_disk_size_gb       -  Amount of boot disk space (in Gb) to give to each machine running each task in this workflow.
#
# This WDL needs to decide whether to use the ``gatk_jar`` or ``gatk_jar_override`` for the jar location.  As of cromwell-0.24,
# this logic *must* go into each task.  Therefore, there is a lot of duplicated code.  This allows users to specify a jar file
# independent of what is in the docker file.  See the README.md for more info.
#
workflow IndexFeatureFile {
    String gatk_docker
    Array[File] variant_vcfs

    File? gatk4_jar_override
    Int?  mem
    Int? preemptible_attempts
    Int? disk_space_gb
    Int? cpu
    Int? boot_disk_size_gb

    scatter ( vcf in variant_vcfs ) {
        call DoIndex {
            input:
                input_vcf                 = vcf,
                gatk_docker               = gatk_docker,
                gatk_override             = gatk4_jar_override,
                mem                       = mem,
                preemptible_attempts      = preemptible_attempts,
                disk_space_gb             = disk_space_gb,
                cpu                       = cpu,
                boot_disk_size_gb         = boot_disk_size_gb
        }
    }

    output {
        Array[File] vcf_out_idxs = DoIndex.vcf_index
    }
}


task DoIndex {

    # ------------------------------------------------
    # Input args:

    # Required:
     File input_vcf

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
     String index_format = if sub(input_vcf, ".*\\.", "") == "vcf" then "idx" else "tbi"
     String timing_output_file = basename(input_vcf) + ".timingInformation.txt"

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
    # Run our command:
     command <<<
         set -e

         startTime=`date +%s.%N`
         echo "StartTime: $startTime" > ${timing_output_file}

         export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

         gatk --java-options "-Xmx${command_mem}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
            IndexFeatureFile \
             -F ${input_vcf} \

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
         preemptible: 0 
         cpu: select_first([cpu, 1])
     }

    # ------------------------------------------------
    # Outputs:
     output {
         File vcf_index = "${input_vcf}.${index_format}"
         File timing_info = timing_output_file
     }
 }
