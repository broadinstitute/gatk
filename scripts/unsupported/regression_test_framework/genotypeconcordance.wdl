# Define the GenotypeConcordance task, to run GenotypeConcordance on a VCF file.
#
# Description of inputs:
#
#   Required:
#     gatk_docker                    - GATK Docker image in which to run
#     call_vcf                       - The produced data set.  Variant Context File (VCF) containing the variants.
#     call_index                     - The produced data set index.  Variant Context File (VCF) index for the given actual_vcf.
#     call_sample                    - Sample in the VCF to compare.
#     truth_vcf                      - The truth data set.  Variant Context File (VCF) containing the variants.
#     truth_index                    - The truth data set index.  Variant Context File (VCF) index for the given expected_vcf.
#     truth_sample                   - Sample in the VCF to use as truth data.
#     intervals                      - File containing intervals over which the comparison should be run.
#     output_base_name               - Base file name for generated output files.
#
#   Optional:
#     gatk4_jar_override             - Override Jar file containing GATK 4.0.  Use this when overriding the docker JAR or when using a backend without docker.
#     mem                            - Amount of memory to give the runtime environment.
#     disk_space_gb                  - Amount of disk space to give the runtime environment.
#     cpu                            - Number of CPUs to give the runtime environment.
#     boot_disk_size_gb              - Amount of boot disk space to give the runtime environment.
#     use_ssd                        - Whether or not to mandate the use of Solid State Drives in the runtime environment.
#     preemptible_attempts           - Number of times the comparison can be preempted.
#
# This WDL needs to decide whether to use the ``gatk_jar`` or ``gatk_jar_override`` for the jar location.  As of cromwell-0.24,
# this logic *must* go into each task.  Therefore, there is a lot of duplicated code.  This allows users to specify a jar file
# independent of what is in the docker file.  See the README.md for more info.
task GenotypeConcordance {

    ####################################################################################
    # Inputs:
    File call_vcf
    File call_index
    String call_sample

    File truth_vcf
    File truth_index
    String truth_sample

    File intervals

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
    Boolean? use_ssd

    ####################################################################################
    # Define default values and set up values for running:
    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 3000
    # WARNING: In the workflow, you should calculate the disk space as an input to this task (disk_space_gb).  Please see [TODO: Link from Jose] for examples.
    Int default_disk_space_gb = 100

    Int default_boot_disk_size_gb = 15

    Int default_min_preemptable_attempts = 2

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem *1000 else default_ram_mb
    Int command_mem = machine_mem - 1000

    ####################################################################################
    # Do the work:
    command {
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        gatk --java-options "-Xmx${command_mem}m" \
            GenotypeConcordance \
                CALL_VCF=${call_vcf} \
                CALL_SAMPLE=${call_sample} \
                TRUTH_VCF=${truth_vcf} \
                TRUTH_SAMPLE=${truth_sample} \
                ${" INTERVALS=" + intervals} \
                INTERSECT_INTERVALS=true \
                OUTPUT_VCF=true
                O=${output_base_name}
    }

    ####################################################################################
    # Runtime Environment:
    runtime {
        cpu: select_first([cpu, 1])
        memory: machine_mem + " MB"
        bootDiskSizeGb: select_first([disk_space_gb, default_disk_space_gb])
        disks: "local-disk " + select_first([boot_disk_size_gb, default_boot_disk_size_gb]) + if use_ssd then " SSD" else " HDD"
        docker: "${gatk_docker}"
        preemptible: select_first([preemptible_attempts, default_min_preemptable_attempts])
    }

    ####################################################################################
    # Outputs:
    output {
        File summary_metrics = "${output_base_name}.genotype_concordance_summary_metrics"
        File detail_metrics = "${output_base_name}.genotype_concordance_detail_metrics"
        File contingency_metrics = "${output_base_name}.genotype_concordance_contingency_metrics"
        File output_vcf = "${output_base_name}.genotype_concordance.vcf.gz"
        File output_vcf_index = "${output_base_name}.genotype_concordance.vcf.gz.tbi"
    }
}
