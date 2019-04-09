# Run IndexFeatureFile on a VCF file. 
#
# Description of inputs:
#
#   Required:
#     gatk_docker                    -  GATK Docker image in which to run
#     variant_vcf                    -  Variant Context File (VCF) containing the variants.
#
#   Optional:
#     gatk4_jar_override             -  Override Jar file containing GATK 4.0.  Use this when overriding the docker JAR or when using a backend without docker.
#
# This WDL needs to decide whether to use the ``gatk_jar`` or ``gatk_jar_override`` for the jar location.  As of cromwell-0.24,
# this logic *must* go into each task.  Therefore, there is a lot of duplicated code.  This allows users to specify a jar file
# independent of what is in the docker file.  See the README.md for more info.
workflow IndexFeatureFile {
    String gatk_docker
    File variant_vcf
    String output_vcf_index

    File? gatk4_jar_override
    Int?  mem
    Int? preemptible_attempts
    Int? disk_space_gb
    Int? cpu
    Int? boot_disk_size_gb

    call DoIndex {
        input:
            input_vcf                 = variant_vcf,
            output_vcf_index          = output_vcf_index,
            gatk_docker               = gatk_docker,
            gatk_override             = gatk4_jar_override,
            mem                       = mem,
            preemptible_attempts      = preemptible_attempts,
            disk_space_gb             = disk_space_gb,
            cpu                       = cpu,
            boot_disk_size_gb          = boot_disk_size_gb
    }

    output {
        File vcf_out_idx = DoIndex.vcf_index
    }
}


task DoIndex {
     # inputs
     File input_vcf

     # outputs
     String output_vcf_index 

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

         gatk --java-options "-Xmx${command_mem}m" \
                IndexFeatureFile \
                 -F ${input_vcf} \
                 -O ${output_vcf_index}
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
         File vcf_index = "${output_vcf_index}"
     }
 }
