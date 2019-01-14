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

    # Default input files for HC and comparison:
    String truth_bucket_location = "gs://haplotypecallerspark-evaluation/groundTruth/"
    String input_bucket_location = "gs://haplotypecallerspark-evaluation/inputData/"
    Array[String] input_bams = [ "G94982.NA12878.bam", "G96830.NA12878.bam", "G96831.NA12878.bam", "G96832.NA12878.bam", "NexPond-359781.bam", "NexPond-412726.bam", "NexPond-445394.bam", "NexPond-472246.bam", "NexPond-506817.bam", "NexPond-538834.bam", "NexPond-572804.bam", "NexPond-603388.bam", "NexPond-633960.bam", "NexPond-656480.bam", "NexPond-679060.bam", "NexPond-412726.chr15.bam"]

    File ref_fasta       = "gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta"
    File ref_fasta_dict  = "gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai"
    File ref_fasta_index = "gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict"

    # Haplotype Caller args:
    File interval_list = "gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list"
    Boolean gvcf_mode
    Float contamination
    Int interval_padding

    # Output bucket name:
    String output_bucket_base_location = "https://console.cloud.google.com/storage/browser/haplotypecallerspark-evaluation/testSets/"

    File? gatk4_jar_override
    Int? mem_gb
    Int? preemptible_attempts
    Int? disk_space_gb
    Int? cpu
    Int? boot_disk_size_gb

    # ------------------------------------------------
    # Call our tool:
    scatter (i in range(length(input_bams))) {

        # Set up variables for this loop:
        File inputBaseName = basename(input_bams[i], ".bam")
        File indexFile = inputBaseName + ".bai"

        String outputName = if gvcf_mode then inputBaseName + ".HC.g.vcf" else inputBaseName + ".HC.vcf"

        String truthBaseName = outputName
        File truthVcf = truth_bucket_location + truthBaseName
        File truthIndex = truth_bucket_location + truthBaseName + ".idx"

        # This is kind of a total hack
        Boolean isUsingIntervals = sub(inputBaseName, "Pond.*", "") == "Nex"
        File? garbageFile
        File? interval_list_final = if !isUsingIntervals then interval_list else garbageFile

        # Only use intervals if we're running on an Exome:
        call tool_wdl.HaplotypeCallerTask {
            input:
                input_bam                 = input_bucket_location + input_bams[i],
                input_bam_index           = input_bucket_location + indexFile,
                ref_fasta                 = ref_fasta,
                ref_fasta_dict            = ref_fasta_dict,
                ref_fasta_index           = ref_fasta_index,

                interval_list             = interval_list_final,
                gvcf_mode                 = gvcf_mode,
                contamination             = contamination,
                interval_padding          = interval_padding,

                out_file_name             = outputName,

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
#                interval_list             = interval_list_final,
#
#                output_base_name          = outputName,
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
#                intervals = interval_list_final,
#
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
                gatk_docker = gatk_docker,
                truth_timing_file = truthVcf + ".timingInformation.txt",
                call_timing_file = HaplotypeCallerTask.timing_info
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
