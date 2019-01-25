# Calls variants on the given germline BAM file using HaplotypeCaller.
#
# Description of inputs:
#
#   Required:
#     String gatk_docker                                -  GATK Docker image in which to run
#
#     Array[File] input_bams                            -  Input reads over which to call small variants with Haplotype Caller.
#
#     File ref_fasta                                    -  Reference FASTA file.
#     File ref_fasta_dict                               -  Reference FASTA file dictionary.
#     File ref_fasta_index                              -  Reference FASTA file index.
#
#   Optional:
#
#     File? interval_list                               -  Interval list over which to call variants on the given BAM file.
#     Boolean? gvcf_mode                                -  Whether to run in GVCF mode (default: false).
#     Float? contamination                              -  Contamination threshold for HaplotypeCaller (default: 0.0).
#     Int? interval_padding                             -  Amount of padding (in bp) to add to each interval you are including (default: 0).
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
workflow HaplotypeCaller {

    # ------------------------------------------------
    # Input args:
    String gatk_docker

    Array[File] input_bams

    File ref_fasta
    File ref_fasta_dict
    File ref_fasta_index

    File? interval_list
    Boolean? gvcf_mode
    Float? contamination
    Int? interval_padding

    File? gatk4_jar_override
    Int?  mem_gb
    Int? preemptible_attempts
    Int? disk_space_gb
    Int? cpu
    Int? boot_disk_size_gb

    # ------------------------------------------------
    # Call our tasks:

    scatter (input_bam in input_bams) {

        File input_base_name = basename(input_bam, ".bam")
        File input_bam_index = sub(input_bam, ".bam$", ".bai")
        String output_name = if (defined(gvcf_mode) && gvcf_mode) then input_base_name + ".HC.g.vcf" else input_base_name + ".HC.vcf"

        Boolean isUsingIntervals = sub(input_base_name, "Pond.*", "") == "Nex"
        File? garbageFile
        File? interval_list_final = if !isUsingIntervals then interval_list else garbageFile

        call HaplotypeCallerTask {
            input:
                input_bam                 = input_bam,
                input_bam_index           = input_bam_index,

                out_file_name             = output_name,

                ref_fasta                 = ref_fasta,
                ref_fasta_dict            = ref_fasta_dict,
                ref_fasta_index           = ref_fasta_index,

                interval_list             = interval_list_final,
                gvcf_mode                 = gvcf_mode,
                contamination             = contamination,
                interval_padding          = interval_padding,

                gatk_docker               = gatk_docker,
                gatk_override             = gatk4_jar_override,
                mem                       = mem_gb,
                preemptible_attempts      = preemptible_attempts,
                disk_space_gb             = disk_space_gb,
                cpu                       = cpu,
                boot_disk_size_gb         = boot_disk_size_gb
        }
    }

    # ------------------------------------------------
    # Outputs:
    output {
        Array[File] vcf_out     = HaplotypeCallerTask.output_vcf
        Array[File] vcf_out_idx = HaplotypeCallerTask.output_vcf_index
        Array[File] timingInfo  = HaplotypeCallerTask.timing_info
    }
}

# ==========================================================================================
# ==========================================================================================
# ==========================================================================================

task HaplotypeCallerTask {

    # ------------------------------------------------
    # Input args:

    # Required:
    File input_bam
    File input_bam_index

    # Output Names:
    String out_file_name
#    String out_file_dir
#    String out_file_name = basename(out_file_dir)

    File ref_fasta
    File ref_fasta_dict
    File ref_fasta_index

    File? interval_list
    Boolean? gvcf_mode
    Float? contamination
    Int? interval_padding

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
    String interval_list_arg = if defined(interval_list) then " -L " else ""
    String contamination_arg = if defined(contamination) then " --contamination " else ""
    String interval_padding_arg = if defined(interval_padding) then " --interval-padding " else ""

    String index_format = if sub(out_file_name, ".*\\.", "") == "vcf" then "idx" else "tbi"

    String timing_output_file = basename(out_file_name) + ".timingInformation.txt"
#    String timing_output_file_gs_dir = out_file_dir + ".timingInformation.txt"

    # ------------------------------------------------
    # Get machine settings:
    Boolean use_ssd = false

    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 3 * 1024
    # WARNING: In the workflow, you should calculate the disk space as an input to this task (disk_space_gb).  Please see [TODO: Link from Jose] for examples.
    Int default_disk_space_gb = 400

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

        gatk --java-options "-Xmx${command_mem}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
            HaplotypeCaller \
                -I ${input_bam} \
                -O ${out_file_name} \
                -R ${ref_fasta} \
                ${interval_list_arg}${default="" sep=" -L " interval_list} \
                ${contamination_arg}${default="" sep=" --contamination " contamination} \
                ${interval_padding_arg}${default="" sep=" --interval-padding " interval_padding} \
                ${true="-ERC GVCF" false="" gvcf_mode}

        endTime=`date +%s.%N`
        echo "EndTime: $endTime" >> ${timing_output_file}
        elapsedTime=`python -c "print( $endTime - $startTime )"`
        echo "Elapsed Time: $elapsedTime" >> ${timing_output_file}

#        # Push to Google Cloud repo for test runs:
#        gsutil cp ${timing_output_file} $timing_output_file_gs_dir}
#        gsutil cp  ${out_file_name} $out_file_dir}
    >>>

    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: gatk_docker
        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + if use_ssd then " SSD" else " HDD"
        bootDiskSizeGb: select_first([boot_disk_size_gb, default_boot_disk_size_gb])
        preemptible: select_first([preemptible_attempts, 0])
        cpu: select_first([cpu, 1])
    }

    # ------------------------------------------------
    # Outputs:
    output {
        File output_vcf       = out_file_name
        File output_vcf_index = "${out_file_name}.${index_format}"
        File timing_info      = timing_output_file
    }
}
