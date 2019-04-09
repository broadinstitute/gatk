# Compare runs done on the haplotype caller.
#
# Description of inputs:
#
#   Required:
#     String gatk_docker                                -  GATK Docker image in which to run, which contains the new HaplotypeCaller to be tested.
#
#   Optional:
#
#      String? analysis_docker                          -  GATK Docker image containing analysis packages to be run.
#                                                          If not provided, will default to gatk_docker
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

    # Set the analysis docker:
    String? analysis_docker
    String analysis_docker_final = if (defined(analysis_docker)) then "" + analysis_docker else gatk_docker

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
        File input_base_name = basename(input_bams[i], ".bam")
        File input_bam_index = input_base_name + ".bai"

        String outputName = if gvcf_mode then input_base_name + ".HC.g.vcf" else input_base_name + ".HC.vcf"

        String truthBaseName = input_base_name + ".vcf.gz"
        File truthVcf = truth_bucket_location + truthBaseName
        File truthIndex = truth_bucket_location + truthBaseName + ".tbi"

        # This is kind of a total hack
        Boolean isUsingIntervals = sub(input_base_name, "Pond.*", "") == "Nex"
        File? emptyFileVariable
        File? interval_list_final = if !isUsingIntervals then interval_list else emptyFileVariable

        # ================================================================================================
        # Run the tool itself:
        # -------------------------------------------------

        # Only use intervals if we're running on an Exome:
        call tool_wdl.HaplotypeCallerTask {
            input:
                input_bam                = input_bucket_location + input_bams[i],
                input_bam_index          = input_bucket_location + input_bam_index,

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

        # ================================================================================================
        # Run the metrics on tool output and the evaluation of the metrics between this run and
        # the "truth"/baseline data:
        # -------------------------------------------------

        call analysis_1_wdl.GenotypeConcordanceTask {
            input:
                call_vcf                  = HaplotypeCallerTask.output_vcf,
                call_index                = HaplotypeCallerTask.output_vcf_index,
                call_sample               = "NA12878",

                truth_vcf                 = truthVcf,
                truth_index               = truthIndex,
                truth_sample              = "NA12878",

                interval_list             = interval_list_final,

                output_base_name          = outputName,

                gatk_docker               = analysis_docker_final,
                #gatk_override             = gatk4_jar_override,
                mem                       = mem_gb,
                preemptible_attempts      = preemptible_attempts,
                disk_space_gb             = disk_space_gb,
                cpu                       = cpu,
                boot_disk_size_gb         = boot_disk_size_gb
        }

        call analysis_2_wdl.Concordance {
            input:
                eval_vcf = HaplotypeCallerTask.output_vcf,
                eval_vcf_idx = HaplotypeCallerTask.output_vcf_index,
                truth_vcf = truthVcf,
                truth_vcf_idx = truthIndex,
                intervals = interval_list_final,

                gatk_docker = analysis_docker_final,
                #gatk_override = gatk4_jar_override,
                mem_gb = mem_gb,
                disk_space_gb = disk_space_gb,
                cpu = cpu,
                boot_disk_size_gb = boot_disk_size_gb,
                preemptible_attempts = preemptible_attempts
        }

#        call analysis_3_wdl.CompareTimingTask {
#            input:
#                gatk_docker = analysis_docker_final,
#                truth_timing_file = truthVcf + ".timingInformation.txt",
#                call_timing_file = HaplotypeCallerTask.timing_info
#        }
    }

    # ------------------------------------------------
    # Outputs:
    output {
        Array[File] vcf_out                                 = HaplotypeCallerTask.output_vcf
        Array[File] vcf_out_idx                             = HaplotypeCallerTask.output_vcf_index

        Array[File] genotypeConcordance_summary_metrics     = GenotypeConcordanceTask.summary_metrics
        Array[File] genotypeConcordance_detail_metrics      = GenotypeConcordanceTask.detail_metrics
        Array[File] genotypeConcordance_contingency_metrics = GenotypeConcordanceTask.contingency_metrics

        Array[File] variantCallerConcordance_fn              = Concordance.fn
        Array[File] variantCallerConcordance_fn_idx          = Concordance.fn_idx
        Array[File] variantCallerConcordance_fp              = Concordance.fp
        Array[File] variantCallerConcordance_fp_idx          = Concordance.fp_idx
        Array[File] variantCallerConcordance_tp              = Concordance.tp
        Array[File] variantCallerConcordance_tp_idx          = Concordance.tp_idx
        Array[File] variantCallerConcordance_ffn             = Concordance.ffn
        Array[File] variantCallerConcordance_ffn_idx         = Concordance.ffn_idx
        Array[File] variantCallerConcordance_summary         = Concordance.summary
        Array[File] variantCallerConcordance_filter_analysis = Concordance.filter_analysis

#        Array[File] timingMetrics                            = CompareTimingTask.timing_diff
    }
}
