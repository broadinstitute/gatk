# Removes a given list of INFO field annotations from the given VCF files.
#
# Description of inputs:
#
#   Required:
#     String gatk_docker                                -  GATK Docker image in which to run
#
#     Array[File] variant_vcfs                          -  Array of Variant Context Files (VCFs) from which to remove INFO field annotations.
#                                                          Assumes that index files are in the same folder as VCF files and of the correct extension.
#     Array[String] info_annotations_to_remove          -  Array of strings with each being the name of an annotation to remove from the INFO field.
#
#   Optional:
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

workflow RemoveInfoFieldAnnotationsFromVcf {
    String gatk_docker

    Array[File] variant_vcfs
    Array[String] info_annotations_to_remove

    File? gatk4_jar_override
    Int?  mem_gb
    Int?  preemptible_attempts
    Int?  disk_space_gb
    Int?  cpu
    Int?  boot_disk_size_gb

    # Run liftover on each input VCF:
    scatter ( vcf_file in variant_vcfs ) {

      # Get the name of this run's index file:
      String index_format = if sub(vcf_file, ".*\\.", "") == "vcf" then "idx" else "tbi"
      File vcf_index = vcf_file + "." + index_format

      # Get the name of this run's VCF file:
      String vcf_extension = sub(vcf_file, "^.*.vcf", ".vcf" )
      String vcf_base_name = basename( vcf_file, vcf_extension )
      String output_vcf_file = vcf_base_name + ".INFO_ANNOTATIONS_FIXED" + vcf_extension

      call SelectVariantsTask {
          input:
              input_vcf_file                              = vcf_file,
              input_vcf_file_index                        = vcf_index,
              output_vcf_file                             = output_vcf_file,
              info_annotations_to_remove                  = info_annotations_to_remove,

              gatk_docker                                 = gatk_docker,
              gatk_override                               = gatk4_jar_override,
              mem                                         = mem_gb,
              preemptible_attempts                        = preemptible_attempts,
              disk_space_gb                               = disk_space_gb,
              cpu                                         = cpu,
              boot_disk_size_gb                           = boot_disk_size_gb
      }
    }

    output {
        Array[File] vcf_file_with_clean_info_fields       = SelectVariantsTask.vcf_with_cleaned_info_field
        Array[File] vcf_file_with_clean_info_fields_index = SelectVariantsTask.vcf_with_cleaned_info_field_index
    }
}

task SelectVariantsTask {

     # ------------------------------------------------
     # Input args:
     File input_vcf_file
     File input_vcf_file_index
     String output_vcf_file
     Array[String] info_annotations_to_remove

     # Runtime Options:
     String gatk_docker
     File? gatk_override
     Int? mem
     Int? preemptible_attempts
     Int? disk_space_gb
     Int? cpu
     Int? boot_disk_size_gb

     String index_format = if sub(input_vcf_file, ".*\\.", "") == "vcf" then "idx" else "tbi"
     String timing_output_file = basename(input_vcf_file) + ".timingInformation.txt"

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

        echo "Disk Space:"
        df -h

        startTime=`date +%s.%N`
        echo "StartTime: $startTime" > ${timing_output_file}

         gatk --java-options "-Xmx${command_mem}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
            SelectVariants \
                -V ${input_vcf_file} \
                -O ${output_vcf_file} \
                --drop-info-annotation ${sep=" --drop-info-annotation " info_annotations_to_remove}

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
         File vcf_with_cleaned_info_field       = "${output_vcf_file}"
         File vcf_with_cleaned_info_field_index = "${output_vcf_file}.${index_format}"
         File timing_info                       = timing_output_file
     }
 }

