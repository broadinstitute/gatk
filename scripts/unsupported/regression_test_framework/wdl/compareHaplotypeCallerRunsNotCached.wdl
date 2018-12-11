# Compare runs done on the haplotype caller.
#
# Description of inputs:
#
#   Required:
#     String gatk_docker                                -  GATK Docker image in which to run, which contains the new HaplotypeCaller to be tested.
#
#
#   Optional:
#
#     HaplotypeCaller:
#       File? interval_list                             -  Interval list over which to call variants on the given BAM file.
#       Boolean? gvcf_mode                              -  Whether to run in GVCF mode (default: false).
#       Float? contamination                            -  Contamination threshold for HaplotypeCaller (default: 0.0).
#       Int? interval_padding                           -  Amount of padding (in bp) to add to each interval you are including (default: 0).
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

import "haplotypeCaller.wdl" as tool_wdl
import "genotypeConcordance.wdl" as analysis_1_wdl
import "variantCallerConcordance.wdl" as analysis_2_wdl
import "compareTiming.wdl" as analysis_3_wdl

workflow ToolComparisonWdl {

    # ------------------------------------------------
    # Input args:
    String gatk_docker
    String baseline_docker = "broadinstitute/gatk:4.0.11.0"
    # Default input files for HC and comparison:
    String truth_bucket_location = "gs://haplotypecallerspark-evaluation/groundTruth/"
    String input_bucket_location = "gs://haplotypecallerspark-evaluation/inputData/"
    Array[String] input_bams = [ "G94982.NA12878.bam", "G96830.NA12878.bam", "G96831.NA12878.bam", "G96832.NA12878.bam", "NexPond-359781.bam", "NexPond-412726.bam", "NexPond-445394.bam", "NexPond-472246.bam", "NexPond-506817.bam", "NexPond-538834.bam", "NexPond-572804.bam", "NexPond-603388.bam", "NexPond-633960.bam", "NexPond-656480.bam", "NexPond-679060.bam" ]

    File ref_fasta       = "gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta"
    File ref_fasta_dict  = "gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai"
    File ref_fasta_index = "gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict"

    # Haplotype Caller args:
    File? interval_list = "gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list"
    Boolean gvcf_mode
    Float contamination
    Int interval_padding

    # Output bucket name:
    String output_bucket_base_location = "gs://haplotypecallerspark-evaluation/testSets/"

    File? gatk4_jar_override
    Int?  mem_gb
    Int? preemptible_attempts
    Int? disk_space_gb
    Int? cpu
    Int? boot_disk_size_gb

    # ------------------------------------------------
    # Create a folder in our output area for this run:
    String output_folder_base = output_bucket_base_location + sub(sub(gatk_docker, "/", "-"), ":", "_") + "/"
    String baseline_output_folder_base = output_bucket_base_location + sub(sub(baseline_docker, "/", "-"), ":", "_") + "/"

    # ------------------------------------------------
    # Call our tool:
    scatter (i in range(length(input_bams))) {

        # Set up variables for this loop:
        File inputBaseName = basename(input_bams[i])
        File indexFile = inputBaseName + ".bai"

        String outputName = if gvcf_mode then sub(basename(input_bams[i]), ".*\\.", "") + ".g.vcf" else sub(basename(input_bams[i]), ".*\\.", "") + ".vcf"

        File truthBaseName = outputName
        File truthVcf = truth_bucket_location + truthBaseName
        File truthIndex = truth_bucket_location + truthBaseName + ".idx"

        call tool_wdl.HaplotypeCallerTask as baseline_run {
            input:
                input_bam                 = input_bucket_location + input_bams[i],
                input_bam_index           = input_bucket_location + indexFile,
                ref_fasta                 = ref_fasta,
                ref_fasta_dict            = ref_fasta_dict,
                ref_fasta_index           = ref_fasta_index,

                interval_list             = interval_list,
                gvcf_mode                 = gvcf_mode,
                contamination             = contamination,
                interval_padding          = interval_padding,

                out_file_name             = baseline_output_folder_base + outputName,

                gatk_docker               = baseline_docker,
                gatk_override             = gatk4_jar_override,
                mem                       = mem_gb,
                preemptible_attempts      = preemptible_attempts,
                disk_space_gb             = disk_space_gb,
                cpu                       = cpu,
                boot_disk_size_gb         = boot_disk_size_gb
        }


        call tool_wdl.HaplotypeCallerTask {
            input:
                input_bam                 = input_bucket_location + input_bams[i],
                input_bam_index           = input_bucket_location + indexFile,
                ref_fasta                 = ref_fasta,
                ref_fasta_dict            = ref_fasta_dict,
                ref_fasta_index           = ref_fasta_index,

                interval_list             = interval_list,
                gvcf_mode                 = gvcf_mode,
                contamination             = contamination,
                interval_padding          = interval_padding,

                out_file_name             = output_folder_base + outputName,

                gatk_docker               = gatk_docker,
                gatk_override             = gatk4_jar_override,
                mem                       = mem_gb,
                preemptible_attempts      = preemptible_attempts,
                disk_space_gb             = disk_space_gb,
                cpu                       = cpu,
                boot_disk_size_gb         = boot_disk_size_gb
        }

#        call analysis_1_wdl.GenotypeConcordanceTask {
#            input:
#                call_vcf                  = HaplotypeCallerTask.output_vcf,
#                call_index                = HaplotypeCallerTask.output_vcf_index,
#                call_sample               = "NA12878",
#
#                truth_vcf                 = truthVcf,
#                truth_index               = truthIndex,
#                truth_sample              = "NA12878",
#
#                interval_list             = interval_list,
#
#                output_base_name          = output_folder_base + outputName,
#
#                gatk_docker               = gatk_docker,
#                gatk_override             = gatk4_jar_override,
#                mem                       = mem_gb,
#                preemptible_attempts      = preemptible_attempts,
#                disk_space_gb             = disk_space_gb,
#                cpu                       = cpu,
#                boot_disk_size_gb         = boot_disk_size_gb
#        }
#
#        call analysis_2_wdl.Concordance {
#            input:
#                eval_vcf = HaplotypeCallerTask.output_vcf,
#                eval_vcf_idx = HaplotypeCallerTask.output_vcf_index,
#                truth_vcf = truthVcf,
#                truth_vcf_idx = truthIndex,
#                intervals = interval_list,
#a
#                gatk_docker = gatk_docker,
#                gatk_override = gatk4_jar_override,
#                mem_gb = mem_gb,
#                disk_space_gb = disk_space_gb,
#                cpu = cpu,
#                boot_disk_size_gb = boot_disk_size_gb,
#                preemptible_attempts = preemptible_attempts
#        }

        call analysis_3_wdl.CompareTimingTask {
            input:
                truth_timing_file = baseline_run.timing_info,
                call_timing_file = HaplotypeCallerTask.timing_info,
                base_timing_output_name = output_folder_base + "timingDiff.txt"
        }
    }

    # ------------------------------------------------
    # Outputs:
    output {
#        File vcf_out     = HaplotypeCallerTask.output_vcf
#        File vcf_out_idx = HaplotypeCallerTask.output_vcf_index
        Array[File] timingMetrics  = CompareTimingTask.timing_diff
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
    File ref_fasta
    File ref_fasta_dict
    File ref_fasta_index

    File? interval_list
    Boolean? gvcf_mode
    Float? contamination
    Int? interval_padding

    # Output Names:
    String out_file_name

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

        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        startTime=`date +%s.%N`
        echo "StartTime: $startTime" > timingInformation.txt

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
        echo "EndTime: $endTime" >> timingInformation.txt
        #elapsedTime=`echo "scale=5;$endTime - $startTime" | bc`
        #echo "Elapsed Time: $elapsedTime" >> timingInformation.txt
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
        File output_vcf       = "${out_file_name}"
        File output_vcf_index = "${out_file_name}.${index_format}"
        File timing_info      = "timingInformation.txt"
    }
}
