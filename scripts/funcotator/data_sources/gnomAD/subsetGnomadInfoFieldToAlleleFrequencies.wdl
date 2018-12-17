# Create a shortened version of gnomAD that only contains sub-population allele frequency data
# from the INFO field.
#
# Description of inputs:
#
#   Required:
#     String gatk_docker        -  GATK Docker image in which to run
#
#     File gnomad_vcf           -  Single-file version of gnomAD to truncate.
#     String gnomad_version     -  Version of gnomAD used.
#
#     File gatk4_jar_override   -  Override Jar file containing GATK 4.  Use this when overriding the docker JAR or when using a backend without docker.
#     Int  mem_gb               -  Amount of memory to give to the machine running each task in this workflow (in gb).
#     Int  preemptible_attempts -  Number of times to allow each task in this workflow to be preempted.
#     Int  disk_space_gb        -  Amount of storage disk space (in Gb) to give to each machine running each task in this workflow.
#     Int  cpu                  -  Number of CPU cores to give to each machine running each task in this workflow.
#     Int  boot_disk_size_gb    -  Amount of boot disk space (in Gb) to give to each machine running each task in this workflow.
#
# This WDL needs to decide whether to use the ``gatk_jar`` or ``gatk_jar_override`` for the jar location.  As of cromwell-0.24,
# this logic *must* go into each task.  Therefore, there is a lot of duplicated code.  This allows users to specify a jar file
# independent of what is in the docker file.  See the README.md for more info.

workflow SubsetGenomicInfoFieldToAlleleFrequencies {
    # inputs
    File gnomad_vcf
    String gnomad_version

    # runtime
    String gatk_docker
    File? gatk4_jar_override
    Int?  mem_gb
    Int? preemptible_attempts
    Int? disk_space_gb
    Int? cpu
    Int? boot_disk_size_gb

    call MakeAlleleFrequencyOnlyGnomadVcf {
        input:
            input_vcf                 = gnomad_vcf,
            output_base_name          = "gnomad-v${gnomad_version}_AF_Info_Only",
            gatk_docker               = gatk_docker,
            gatk_override             = gatk4_jar_override,
            mem_gb                    = mem_gb,
            preemptible_attempts      = preemptible_attempts,
            disk_space_gb             = disk_space_gb,
            cpu                       = cpu,
            boot_disk_size_gb         = boot_disk_size_gb
    }

    output {
        File gnomad_info_subset_vcf = MakeAlleleFrequencyOnlyGnomadVcf.output_vcf
        File gnomad_info_subset_vcf_idx = MakeAlleleFrequencyOnlyGnomadVcf.output_vcf_idx
    }
}

task MakeAlleleFrequencyOnlyGnomadVcf {

    # ------------------------------------------------
    # Input args:
    File input_vcf
    String output_base_name

    File? gatk_override

    # Runtime Options:
    String gatk_docker
    Int? mem_gb
    Int? preemptible_attempts
    Int? disk_space_gb
    Int? cpu
    Int? boot_disk_size_gb

    # ------------------------------------------------
    # Get machine settings:
    Boolean use_ssd = false

    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 3 * 1024
    # WARNING: In the workflow, you should calculate the disk space as an input to this task (disk_space_gb).  Please see [TODO: Link from Jose] for examples.
    Int default_disk_space_gb = 100

    Int default_boot_disk_size_gb = 15

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem_gb) then mem_gb * 1024 else default_ram_mb
    Int command_mem = machine_mem - 1024

    Int runtime_cpu = select_first([cpu, 4])

    # ------------------------------------------------
    # Run our command:
    command {
        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        # Save off the header for later:
        grep '^#' ${input_vcf} > header &

        # Get all lines in the file except the header:
        # Preserve all fields before INFO, Grab only the AF annotation from the INFO Field
        # Replace QUAL (6th) column with '.' (empty):
        grep -v "^#" ${input_vcf} | sed -e 's#\(.*\t.*\t.*\t.*\t.*\)\t.*\t\(.*\)\t.*;\(AF_afr=[0-9e\.+-]*\);.*;\(AF_amr=[0-9e\.+-]*\);.*;\(AF_eas=[0-9e\.+-]*\);.*;\(AF_nfe=[0-9e\.+-]*\);.*;\(AF_fin=[0-9e\.+-]*\);.*;\(AF_asj=[0-9e\.+-]*\);.*;\(AF_oth=[0-9e\.+-]*\);.*#\1\t.\t\2\t\3;\4;\5;\6;\7;\8;\9#g' > simplified_body &

        # Wait for background processes to finish:
        wait

        # Consolidate files:
        cat header simplified_body > ${output_base_name}.vcf

        # Zip the VCF:
        bgzip --threads ${runtime_cpu} ${output_base_name}.vcf

        # Index output file:
        gatk --java-options "-Xmx${command_mem}g" IndexFeatureFile -F ${output_base_name}.vcf.gz

        # Cleanup:
        rm -f header body simplified_info simplified_body simplified.vcf simplified.vcf.idx
    }

    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: gatk_docker
        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + if use_ssd then " SSD" else " HDD"
        bootDiskSizeGb: select_first([boot_disk_size_gb, default_boot_disk_size_gb])
        preemptible: 0
        cpu: runtime_cpu
    }

    # ------------------------------------------------
    # Outputs:
    output {
        File output_vcf = "${output_base_name}.vcf.gz"
        File output_vcf_idx = "${output_base_name}.vcf.gz.tbi"
    }
}
