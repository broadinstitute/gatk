# Lifts over a VCF file from one reference version to another.
#
# Description of inputs:
#
#   Required:
#     String gatk_docker                                -  GATK Docker image in which to run
#
#     Array[File] variant_vcfs                          -  Array of Variant Context Files (VCFs) containing the variants.
#     File chain_file                                   -  Chain file specifying the conversion from the variant reference to the target reference.
#     File target_reference_sequence_fasta_file         -  Reference FASTA file for the target reference (not the reference for the given VCF file).
#     File target_reference_sequence_fasta_file_index   -  Index for reference FASTA file for the target reference (not the reference for the given VCF file).
#     File target_reference_sequence_fasta_file_dict    -  Sequence dictionary for reference FASTA file for the target reference (not the reference for the given VCF file).
#     String lifted_over_vcf_name                       -  Output name for the lifted over VCF file.
#     String lifted_over_rejects_vcf_name               -  Output name for the lifted over rejects file.
#
#   Optional:
#     Boolean warn_on_missing_contig                    -  Whether to create a warning message when a contig is missing.
#     Boolean write_original_position                   -  Whether to write the original position as an annotation in the resulting lifted over VCF.
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

workflow LiftoverVcf {
    String gatk_docker

    Array[File] variant_vcfs
    File chain_file
    File target_reference_sequence_fasta_file
    File target_reference_sequence_fasta_file_index
    File target_reference_sequence_fasta_file_dict
    String lifted_over_vcf_name
    String lifted_over_rejects_vcf_name

    Boolean? warn_on_missing_contig
    Boolean? write_original_position

    File? gatk4_jar_override
    Int?  mem_gb
    Int?  preemptible_attempts
    Int?  disk_space_gb
    Int?  cpu
    Int?  boot_disk_size_gb

    # Run liftover on each input VCF:
    scatter ( vcf in variant_vcfs ) {
      call LiftoverVcfTask {
          input:
              input_vcf_file                              = vcf,
              chain_file                                  = chain_file,
              target_reference_sequence_fasta_file        = target_reference_sequence_fasta_file,
              target_reference_sequence_fasta_file_index  = target_reference_sequence_fasta_file_index,
              target_reference_sequence_fasta_file_dict   = target_reference_sequence_fasta_file_dict,
              lifted_over_vcf_name                        = lifted_over_vcf_name,
              lifted_over_rejects_vcf_name                = lifted_over_rejects_vcf_name,

              warn_on_missing_contig                      = warn_on_missing_contig,
              write_original_position                     = write_original_position,

              gatk_docker                                 = gatk_docker,
              gatk_override                               = gatk4_jar_override,
              mem                                         = mem_gb,
              preemptible_attempts                        = preemptible_attempts,
              disk_space_gb                               = disk_space_gb,
              cpu                                         = cpu,
              boot_disk_size_gb                           = boot_disk_size_gb
      }

      call IndexFeatureFileTask as IndexContigFiles {
        input:
          input_vcf            = LiftoverVcfTask.vcf_file,

          gatk_docker          = gatk_docker,
          gatk_override        = gatk4_jar_override,
          mem                  = mem_gb,
          preemptible_attempts = preemptible_attempts,
          disk_space_gb        = disk_space_gb,
          cpu                  = cpu,
          boot_disk_size_gb    = boot_disk_size_gb
      }
    }

    # Consolidate output VCFs:
    scatter ( vcf_out_pair in [ [LiftoverVcfTask.lifted_over_vcf, lifted_over_vcf_name], [LiftoverVcfTask.lifted_over_rejects_vcf , lifted_over_rejects_vcf_name] ] ) {
      call MergeVcfsTask {
        input:
          input_vcfs           = vcf_out_pair[0],
          output_vcf_file      = vcf_out_pair[1],

          gatk_docker          = gatk_docker,
          gatk_override        = gatk4_jar_override,
          mem                  = mem_gb,
          preemptible_attempts = preemptible_attempts,
          disk_space_gb        = disk_space_gb,
          cpu                  = cpu,
          boot_disk_size_gb    = boot_disk_size_gb
      }
      call IndexFeatureFileTask {
        input:
          input_vcf            = MergeVcfsTask.vcf_file,

          gatk_docker          = gatk_docker,
          gatk_override        = gatk4_jar_override,
          mem                  = mem_gb,
          preemptible_attempts = preemptible_attempts,
          disk_space_gb        = disk_space_gb,
          cpu                  = cpu,
          boot_disk_size_gb    = boot_disk_size_gb
     }
    }

    output {
        File lifted_over_vcf_file               = MergeVcfsTask.vcf_file[0]
        File lifted_over_vcf_file_index         = IndexFeatureFileTask.vcf_index[0]
        File lifted_over_vcf_rejects_file       = MergeVcfsTask.vcf_file[1]
        File lifted_over_vcf_rejects_file_index = IndexFeatureFileTask.vcf_index[1]

        Array[File] lifted_over_vcf_contig_files        = LiftoverVcfTask.vcf_file
        Array[File] lifted_over_vcf_contig_file_indexes = IndexContigFiles.vcf_index
    }
}

task indexFastaFileTask {
  
  File input_fasta_file

  # runtime
  String gatk_docker
  File? gatk_override
  Int? mem
  Int? preemptible_attempts
  Int? disk_space_gb
  Int? cpu
  Int? boot_disk_size_gb

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
    samtools faidx ${input_fasta_file}

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
      File vcf_index         = "${input_fasta_file}.idx"
  }

}


task LiftoverVcfTask {
     # inputs
     File input_vcf_file
     File chain_file
     File target_reference_sequence_fasta_file
     File target_reference_sequence_fasta_file_index
     File target_reference_sequence_fasta_file_dict

     String lifted_over_vcf_name
     String lifted_over_rejects_vcf_name

     Boolean? warn_on_missing_contig
     Boolean? write_original_position

     # runtime
     String gatk_docker
     File? gatk_override
     Int? mem
     Int? preemptible_attempts
     Int? disk_space_gb
     Int? cpu
     Int? boot_disk_size_gb

     # ------------------------------------------------
     # Process input args:
     String warn_on_missing_contig_arg = if defined(warn_on_missing_contig) then " --WARN_ON_MISSING_CONTIG " else ""
     String write_original_position_arg = if defined(write_original_position) then " --WRITE_ORIGINAL_POSITION " else ""

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
        echo "StartTime: $startTime" > timingInformation.txt

         gatk --java-options "-Xmx${command_mem}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
            LiftoverVcf \
                -I ${input_vcf_file} \
                -O ${lifted_over_vcf_name} \
                -CHAIN ${chain_file} \
                -REJECT ${lifted_over_rejects_vcf_name} \
                -R ${target_reference_sequence_fasta_file} \
                --MAX_RECORDS_IN_RAM 1000 \
                ${warn_on_missing_contig_arg}${default="" sep=" --WARN_ON_MISSING_CONTIG " warn_on_missing_contig} \
                ${write_original_position_arg}${default="" sep=" --WRITE_ORIGINAL_POSITION " write_original_position}

        endTime=`date +%s.%N`
        echo "EndTime: $endTime" >> timingInformation.txt
        elapsedTime=`echo "scale=5;$endTime - $startTime" | bc`
        echo "Elapsed Time: $elapsedTime" >> timingInformation.txt
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
         File lifted_over_vcf         = "${lifted_over_vcf_name}"
         File lifted_over_rejects_vcf = "${lifted_over_rejects_vcf_name}"
     }
 }

task MergeVcfsTask {
     # inputs
     Array[File] input_vcfs

     # outputs
     String output_vcf_file

     # runtime
     String gatk_docker
     File? gatk_override
     Int? mem
     Int? preemptible_attempts
     Int? disk_space_gb
     Int? cpu
     Int? boot_disk_size_gb

     Boolean use_ssd = false

     # You may have to change the following two parameter values depending on the task requirements
     Int default_ram_mb = 3000
     # WARNING: In the workflow, you should calculate the disk space as an input to this task (disk_space_gb).  Please see [TODO: Link from Jose] for examples.
     Int default_disk_space_gb = 100

     Int default_boot_disk_size_gb = 15

     # Mem is in units of GB but our command and memory runtime values are in MB
     Int machine_mem = if defined(mem) then mem *1000 else default_ram_mb
     Int command_mem = machine_mem - 1000

     command <<<
         set -e
         export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        startTime=`date +%s.%N`
        echo "StartTime: $startTime" > timingInformation.txt

         gatk --java-options "-Xmx${command_mem}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
            MergeVcfs \
             -I ${sep=' -I ' input_vcfs} \
             -O ${output_vcf_file}

        endTime=`date +%s.%N`
        echo "EndTime: $endTime" >> timingInformation.txt
        elapsedTime=`echo "scale=5;$endTime - $startTime" | bc`
        echo "Elapsed Time: $elapsedTime" >> timingInformation.txt
     >>>

     runtime {
         docker: gatk_docker
         memory: machine_mem + " MB"
         disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + if use_ssd then " SSD" else " HDD"
         bootDiskSizeGb: select_first([boot_disk_size_gb, default_boot_disk_size_gb])
         preemptible: 0 
         cpu: select_first([cpu, 1])
     }

     output {
         File vcf_file = "${output_vcf_file}"
     }
 }


task IndexFeatureFileTask {
     # inputs
     File input_vcf

     # runtime
     String gatk_docker
     File? gatk_override
     Int? mem
     Int? preemptible_attempts
     Int? disk_space_gb
     Int? cpu
     Int? boot_disk_size_gb

     String index_format = if sub(input_vcf, ".*\\.", "") == "vcf" then "idx" else "tbi"

     Boolean use_ssd = false

     # You may have to change the following two parameter values depending on the task requirements
     Int default_ram_mb = 3000
     # WARNING: In the workflow, you should calculate the disk space as an input to this task (disk_space_gb).  Please see [TODO: Link from Jose] for examples.
     Int default_disk_space_gb = 100

     Int default_boot_disk_size_gb = 15

     # Mem is in units of GB but our command and memory runtime values are in MB
     Int machine_mem = if defined(mem) then mem *1000 else default_ram_mb
     Int command_mem = machine_mem - 1000

     command <<<
        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        startTime=`date +%s.%N`
        echo "StartTime: $startTime" > timingInformation.txt

        gatk --java-options "-Xmx${command_mem}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
            IndexFeatureFile \
                -F ${input_vcf} \

        endTime=`date +%s.%N`
        echo "EndTime: $endTime" >> timingInformation.txt
        elapsedTime=`echo "scale=5;$endTime - $startTime" | bc`
        echo "Elapsed Time: $elapsedTime" >> timingInformation.txt
     >>>

     runtime {
         docker: gatk_docker
         memory: machine_mem + " MB"
         disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + if use_ssd then " SSD" else " HDD"
         bootDiskSizeGb: select_first([boot_disk_size_gb, default_boot_disk_size_gb])
         preemptible: 0 
         cpu: select_first([cpu, 1])
     }

     output {
         File vcf_index = "${input_vcf}.${index_format}"
     }
} 
