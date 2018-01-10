# Workflow for running the GATK CNV pipeline on a matched pair. Supports both WGS and WES.
#
# Notes:
#
# - The interval-list file is required for both WGS and WES workflows and should be a Picard or GATK-style interval list.
#   These intervals will be padded on both sides by the amount specified by PreprocessIntervals.padding (default 250)
#   and split into bins of length specified by PreprocessIntervals.bin_length (default 1000; specify 0 to skip binning,
#   e.g. for WES).  For WGS, the intervals should simply cover the autosomal chromosomes (sex chromosomes may be
#   included, but care should be taken to 1) avoid creating panels of mixed sex, and 2) denoise case samples only
#   with panels containing individuals of the same sex as the case samples).
#
# - The sites file (common_sites) should be a Picard or GATK-style interval list.  This is a list of sites
#   of known variation at which allelic counts will be collected for use in modeling minor-allele fractions.
#
# - Example invocation:
#
#       java -jar cromwell.jar run cnv_somatic_pair_workflow.wdl -i myParameters.json
#
#   See cnv_somatic_pair_workflow_template.json for a template json file to modify with your own parameters (please save
#   your modified version with a different filename and do not commit to the gatk repository).
#
#############

import "cnv_common_tasks.wdl" as CNVTasks
import "cnv_somatic_oncotator_workflow.wdl" as CNVOncotator

workflow CNVSomaticPairWorkflow {
    File common_sites
    File intervals
    File tumor_bam
    File tumor_bam_idx
    File normal_bam
    File normal_bam_idx
    File read_count_pon
    File ref_fasta_dict
    File ref_fasta_fai
    File ref_fasta
    String gatk_docker
    File? gatk4_jar_override

    # For running oncotator
    Boolean is_run_oncotator = false

    # Ignored if not running oncotator
    String oncotator_docker = "broadinstitute/oncotator:1.9.5.0-eval-gatk-protected"

    # Use as a last resort to increase the disk given to every task in case of ill behaving data
    Int? emergency_extra_disk

    Int ref_size = ceil(size(ref_fasta, "GB") + size(ref_fasta_dict, "GB") + size(ref_fasta_fai, "GB"))
    Int read_count_pon_size = ceil(size(read_count_pon, "GB"))
    Int tumor_bam_size = ceil(size(tumor_bam, "GB") + size(tumor_bam_idx, "GB"))
    Int normal_bam_size = ceil(size(normal_bam, "GB") + size(normal_bam_idx, "GB"))

    Int gatk4_override_size = if defined(gatk4_jar_override) then ceil(size(gatk4_jar_override, "GB")) else 0
    # This is added to every task as padding, should increase if systematically you need more disk for every call
    Int disk_pad = 20 + ceil(size(intervals, "GB")) + ceil(size(common_sites, "GB")) + gatk4_override_size + select_first([emergency_extra_disk,0])

    Int process_disk = ref_size + disk_pad
    call CNVTasks.PreprocessIntervals {
        input:
            intervals = intervals,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            ref_fasta_dict = ref_fasta_dict,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            disk_space_gb = process_disk
    }

    Int collect_counts_tumor_disk = tumor_bam_size + ceil(size(PreprocessIntervals.preprocessed_intervals, "GB")) + disk_pad
    call CNVTasks.CollectCounts as CollectCountsTumor {
        input:
            intervals = PreprocessIntervals.preprocessed_intervals,
            bam = tumor_bam,
            bam_idx = tumor_bam_idx,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            disk_space_gb = collect_counts_tumor_disk
    }

    Int collect_counts_normal_disk = normal_bam_size + ceil(size(PreprocessIntervals.preprocessed_intervals, "GB")) + disk_pad
    call CNVTasks.CollectCounts as CollectCountsNormal {
        input:
            intervals = PreprocessIntervals.preprocessed_intervals,
            bam = normal_bam,
            bam_idx = normal_bam_idx,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            disk_space_gb = collect_counts_normal_disk
    }

    Int collect_allelic_counts_tumor_disk = tumor_bam_size + ref_size + disk_pad
    call CNVTasks.CollectAllelicCounts as CollectAllelicCountsTumor {
        input:
            common_sites = common_sites,
            bam = tumor_bam,
            bam_idx = tumor_bam_idx,
            ref_fasta = ref_fasta,
            ref_fasta_dict = ref_fasta_dict,
            ref_fasta_fai = ref_fasta_fai,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            disk_space_gb = collect_allelic_counts_tumor_disk
    }

    Int collect_allelic_counts_normal_disk = normal_bam_size + ref_size + disk_pad
    call CNVTasks.CollectAllelicCounts as CollectAllelicCountsNormal {
        input:
            common_sites = common_sites,
            bam = normal_bam,
            bam_idx = normal_bam_idx,
            ref_fasta = ref_fasta,
            ref_fasta_dict = ref_fasta_dict,
            ref_fasta_fai = ref_fasta_fai,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            disk_space_gb = collect_allelic_counts_normal_disk
    }

    Int denoise_read_counts_tumor_disk = read_count_pon_size + ceil(size(CollectCountsTumor.counts, "GB")) + disk_pad
    call DenoiseReadCounts as DenoiseReadCountsTumor {
        input:
            entity_id = CollectCountsTumor.entity_id,
            read_counts = CollectCountsTumor.counts,
            read_count_pon = read_count_pon,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            disk_space_gb = denoise_read_counts_tumor_disk
    }

    Int denoise_read_counts_normal_disk = read_count_pon_size + ceil(size(CollectCountsNormal.counts, "GB")) + disk_pad
    call DenoiseReadCounts as DenoiseReadCountsNormal {
        input:
            entity_id = CollectCountsNormal.entity_id,
            read_counts = CollectCountsNormal.counts,
            read_count_pon = read_count_pon,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            disk_space_gb = denoise_read_counts_normal_disk
    }

    Int model_segments_disk = ceil(size(DenoiseReadCountsTumor.denoised_copy_ratios, "GB")) + ceil(size(CollectAllelicCountsTumor.allelic_counts, "GB")) + ceil(size(CollectAllelicCountsNormal.allelic_counts, "GB")) + disk_pad
    call ModelSegments as ModelSegmentsTumor {
        input:
            entity_id = CollectCountsTumor.entity_id,
            denoised_copy_ratios = DenoiseReadCountsTumor.denoised_copy_ratios,
            allelic_counts = CollectAllelicCountsTumor.allelic_counts,
            normal_allelic_counts = CollectAllelicCountsNormal.allelic_counts,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            disk_space_gb = model_segments_disk
    }

    call ModelSegments as ModelSegmentsNormal {
        input:
            entity_id = CollectCountsNormal.entity_id,
            denoised_copy_ratios = DenoiseReadCountsNormal.denoised_copy_ratios,
            allelic_counts = CollectAllelicCountsNormal.allelic_counts,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            disk_space_gb = model_segments_disk
    }

    Int copy_ratio_segments_tumor_disk = ceil(size(DenoiseReadCountsTumor.denoised_copy_ratios, "GB")) + ceil(size(ModelSegmentsTumor.copy_ratio_only_segments, "GB")) + disk_pad
    call CallCopyRatioSegments as CallCopyRatioSegmentsTumor {
        input:
            entity_id = CollectCountsTumor.entity_id,
            copy_ratio_segments = ModelSegmentsTumor.copy_ratio_only_segments,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            disk_space_gb = copy_ratio_segments_tumor_disk
    }

    Int copy_ratio_segments_normal_disk = ceil(size(DenoiseReadCountsNormal.denoised_copy_ratios, "GB")) + ceil(size(ModelSegmentsNormal.copy_ratio_only_segments, "GB")) + disk_pad
    call CallCopyRatioSegments as CallCopyRatioSegmentsNormal {
        input:
            entity_id = CollectCountsNormal.entity_id,
            copy_ratio_segments = ModelSegmentsNormal.copy_ratio_only_segments,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            disk_space_gb = copy_ratio_segments_normal_disk
    }

    # The F=files from other tasks are small enough to just combine into one disk variable and pass to the tumor plotting tasks
    Int plot_tumor_disk = ref_size + ceil(size(DenoiseReadCountsTumor.standardized_copy_ratios, "GB")) + ceil(size(DenoiseReadCountsTumor.denoised_copy_ratios, "GB")) + ceil(size(ModelSegmentsTumor.het_allelic_counts, "GB")) + ceil(size(ModelSegmentsTumor.modeled_segments, "GB")) + disk_pad
    call PlotDenoisedCopyRatios as PlotDenoisedCopyRatiosTumor {
        input:
            entity_id = CollectCountsTumor.entity_id,
            standardized_copy_ratios = DenoiseReadCountsTumor.standardized_copy_ratios,
            denoised_copy_ratios = DenoiseReadCountsTumor.denoised_copy_ratios,
            ref_fasta_dict = ref_fasta_dict,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            disk_space_gb = plot_tumor_disk
    }

    # The files from other tasks are small enough to just combine into one disk variable and pass to the normal plotting tasks
    Int plot_normal_disk = ref_size + ceil(size(DenoiseReadCountsNormal.standardized_copy_ratios, "GB")) + ceil(size(DenoiseReadCountsNormal.denoised_copy_ratios, "GB")) + ceil(size(ModelSegmentsNormal.het_allelic_counts, "GB")) + ceil(size(ModelSegmentsNormal.modeled_segments, "GB")) + disk_pad
    call PlotDenoisedCopyRatios as PlotDenoisedCopyRatiosNormal {
        input:
            entity_id = CollectCountsNormal.entity_id,
            standardized_copy_ratios = DenoiseReadCountsNormal.standardized_copy_ratios,
            denoised_copy_ratios = DenoiseReadCountsNormal.denoised_copy_ratios,
            ref_fasta_dict = ref_fasta_dict,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            disk_space_gb = plot_normal_disk
    }

    call PlotModeledSegments as PlotModeledSegmentsTumor {
        input:
            entity_id = CollectCountsTumor.entity_id,
            denoised_copy_ratios = DenoiseReadCountsTumor.denoised_copy_ratios,
            het_allelic_counts = ModelSegmentsTumor.het_allelic_counts,
            modeled_segments = ModelSegmentsTumor.modeled_segments,
            ref_fasta_dict = ref_fasta_dict,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            disk_space_gb = plot_tumor_disk
    }

    call PlotModeledSegments as PlotModeledSegmentsNormal {
        input:
            entity_id = CollectCountsNormal.entity_id,
            denoised_copy_ratios = DenoiseReadCountsNormal.denoised_copy_ratios,
            het_allelic_counts = ModelSegmentsNormal.het_allelic_counts,
            modeled_segments = ModelSegmentsNormal.modeled_segments,
            ref_fasta_dict = ref_fasta_dict,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            disk_space_gb = plot_normal_disk
    }

    if (is_run_oncotator) {
        call CNVOncotator.CNVOncotatorWorkflow as CNVOncotatorWorkflow {
            input:
                 called_file = CallCopyRatioSegmentsTumor.called_copy_ratio_segments,
                 oncotator_docker = oncotator_docker
        }
    }

    output {
        File preprocessed_intervals = PreprocessIntervals.preprocessed_intervals

        File read_counts_tumor = CollectCountsTumor.counts
        File read_counts_entity_id_tumor = CollectCountsTumor.entity_id
        File allelic_counts_tumor = CollectAllelicCountsTumor.allelic_counts
        File allelic_counts_entity_id_tumor = CollectAllelicCountsTumor.entity_id
        File denoised_copy_ratios_tumor = DenoiseReadCountsTumor.denoised_copy_ratios
        File standardized_copy_ratios_tumor = DenoiseReadCountsTumor.standardized_copy_ratios
        File het_allelic_counts_tumor = ModelSegmentsTumor.het_allelic_counts
        File normal_het_allelic_counts_tumor = ModelSegmentsTumor.normal_het_allelic_counts
        File copy_ratio_only_segments_tumor = ModelSegmentsTumor.copy_ratio_only_segments
        File modeled_segments_begin_tumor = ModelSegmentsTumor.modeled_segments_begin
        File copy_ratio_parameters_begin_tumor = ModelSegmentsTumor.copy_ratio_parameters_begin
        File allele_fraction_parameters_begin_tumor = ModelSegmentsTumor.allele_fraction_parameters_begin
        File modeled_segments_tumor = ModelSegmentsTumor.modeled_segments
        File copy_ratio_parameters_tumor = ModelSegmentsTumor.copy_ratio_parameters
        File allele_fraction_parameters_tumor = ModelSegmentsTumor.allele_fraction_parameters
        File called_copy_ratio_segments_tumor = CallCopyRatioSegmentsTumor.called_copy_ratio_segments
        File denoised_copy_ratios_plot_tumor = PlotDenoisedCopyRatiosTumor.denoised_copy_ratios_plot
        File denoised_copy_ratios_lim_4_plot_tumor = PlotDenoisedCopyRatiosTumor.denoised_copy_ratios_lim_4_plot
        File standardized_MAD_tumor = PlotDenoisedCopyRatiosTumor.standardized_MAD
        File denoised_MAD_tumor = PlotDenoisedCopyRatiosTumor.denoised_MAD
        File delta_MAD_tumor = PlotDenoisedCopyRatiosTumor.delta_MAD
        File scaled_delta_MAD_tumor = PlotDenoisedCopyRatiosTumor.scaled_delta_MAD
        File modeled_segments_plot_tumor = PlotModeledSegmentsTumor.modeled_segments_plot

        File read_counts_normal = CollectCountsNormal.counts
        File read_counts_entity_id_normal = CollectCountsNormal.entity_id
        File allelic_counts_normal = CollectAllelicCountsNormal.allelic_counts
        File allelic_counts_entity_id_normal = CollectAllelicCountsNormal.entity_id
        File denoised_copy_ratios_normal = DenoiseReadCountsNormal.denoised_copy_ratios
        File standardized_copy_ratios_normal = DenoiseReadCountsNormal.standardized_copy_ratios
        File het_allelic_counts_normal = ModelSegmentsNormal.het_allelic_counts
        File normal_het_allelic_counts_normal = ModelSegmentsNormal.normal_het_allelic_counts
        File copy_ratio_only_segments_normal = ModelSegmentsNormal.copy_ratio_only_segments
        File modeled_segments_begin_normal = ModelSegmentsNormal.modeled_segments_begin
        File copy_ratio_parameters_begin_normal = ModelSegmentsNormal.copy_ratio_parameters_begin
        File allele_fraction_parameters_begin_normal = ModelSegmentsNormal.allele_fraction_parameters_begin
        File modeled_segments_normal = ModelSegmentsNormal.modeled_segments
        File copy_ratio_parameters_normal = ModelSegmentsNormal.copy_ratio_parameters
        File allele_fraction_parameters_normal = ModelSegmentsNormal.allele_fraction_parameters
        File called_copy_ratio_segments_normal = CallCopyRatioSegmentsNormal.called_copy_ratio_segments
        File denoised_copy_ratios_plot_normal = PlotDenoisedCopyRatiosNormal.denoised_copy_ratios_plot
        File denoised_copy_ratios_lim_4_plot_normal = PlotDenoisedCopyRatiosNormal.denoised_copy_ratios_lim_4_plot
        File standardized_MAD_normal = PlotDenoisedCopyRatiosNormal.standardized_MAD
        File denoised_MAD_normal = PlotDenoisedCopyRatiosNormal.denoised_MAD
        File delta_MAD_normal = PlotDenoisedCopyRatiosNormal.delta_MAD
        File scaled_delta_MAD_normal = PlotDenoisedCopyRatiosNormal.scaled_delta_MAD
        File modeled_segments_plot_normal = PlotModeledSegmentsNormal.modeled_segments_plot

        File oncotated_called_file_tumor = select_first([CNVOncotatorWorkflow.oncotated_called_file, "null"])
        File oncotated_called_gene_list_file_tumor = select_first([CNVOncotatorWorkflow.oncotated_called_gene_list_file, "null"])
    }
}

task DenoiseReadCounts {
    String entity_id
    File read_counts
    File read_count_pon
    Int? number_of_eigensamples #use all eigensamples in panel by default
    File? gatk4_jar_override

    # Runtime parameters
    Int? mem
    String gatk_docker
    Int? preemptible_attempts
    Int disk_space_gb

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1000 else 13000
    Int command_mem = machine_mem - 1000

    command <<<
        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk4_jar_override}

        gatk --java-options "-Xmx${command_mem}m" DenoiseReadCounts \
            --input ${read_counts} \
            --count-panel-of-normals ${read_count_pon} \
            ${"--number-of-eigensamples " + number_of_eigensamples} \
            --standardized-copy-ratios ${entity_id}.standardizedCR.tsv \
            --denoised-copy-ratios ${entity_id}.denoisedCR.tsv
    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: machine_mem + " MB"
        disks: "local-disk " + disk_space_gb + " HDD"
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File standardized_copy_ratios = "${entity_id}.standardizedCR.tsv"
        File denoised_copy_ratios = "${entity_id}.denoisedCR.tsv"
    }
}

task ModelSegments {
    String entity_id
    File denoised_copy_ratios
    File allelic_counts
    File? normal_allelic_counts
    Int? max_num_segments_per_chromosome
    Int? min_total_allele_count
    Float? genotyping_homozygous_log_ratio_threshold
    Float? genotyping_base_error_rate
    Float? kernel_variance_copy_ratio
    Float? kernel_variance_allele_fraction
    Float? kernel_scaling_allele_fraction
    Int? kernel_approximation_dimension
    Array[Int]+? window_sizes = [8, 16, 32, 64, 128, 256]
    Float? num_changepoints_penalty_factor
    Float? minor_allele_fraction_prior_alpha
    Int? num_samples_copy_ratio
    Int? num_burn_in_copy_ratio
    Int? num_samples_allele_fraction
    Int? num_burn_in_allele_fraction
    Float? smoothing_threshold_copy_ratio
    Float? smoothing_threshold_allele_fraction
    Int? max_num_smoothing_iterations
    Int? num_smoothing_iterations_per_fit
    String? output_dir
    File? gatk4_jar_override

    # Runtime parameters
    Int? mem
    String gatk_docker
    Int? preemptible_attempts
    Int disk_space_gb

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1000 else 13000
    # ModelSegments seems to need at least 3GB of overhead to run
    Int command_mem = machine_mem - 3000

    # If optional output_dir not specified, use "out"
    String output_dir_ = select_first([output_dir, "out"])

    command <<<
        set -e
        mkdir ${output_dir_}
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk4_jar_override}

        gatk --java-options "-Xmx${command_mem}m" ModelSegments \
            --denoised-copy-ratios ${denoised_copy_ratios} \
            --allelic-counts ${allelic_counts} \
            ${"--normal-allelic-counts " + normal_allelic_counts} \
            --minimum-total-allele-count ${default="30" min_total_allele_count} \
            --genotyping-homozygous-log-ratio-threshold ${default="-10.0" genotyping_homozygous_log_ratio_threshold} \
            --genotyping-base-error-rate ${default="0.05" genotyping_base_error_rate} \
            --maximum-number-of-segments-per-chromosome ${default="1000" max_num_segments_per_chromosome} \
            --kernel-variance-copy-ratio ${default="0.0" kernel_variance_copy_ratio} \
            --kernel-variance-allele-fraction ${default="0.025" kernel_variance_allele_fraction} \
            --kernel-scaling-allele-fraction ${default="1.0" kernel_scaling_allele_fraction} \
            --kernel-approximation-dimension ${default="100" kernel_approximation_dimension} \
            --window-size ${sep= " --window-size " window_sizes} \
            --number-of-changepoints-penalty-factor ${default="1.0" num_changepoints_penalty_factor} \
            --minor-allele-fraction-prior-alpha ${default="25.0" minor_allele_fraction_prior_alpha} \
            --number-of-samples-copy-ratio ${default=100 num_samples_copy_ratio} \
            --number-of-burn-in-samples-copy-ratio ${default=50 num_burn_in_copy_ratio} \
            --number-of-samples-allele-fraction ${default=100 num_samples_allele_fraction} \
            --number-of-burn-in-samples-allele-fraction ${default=50 num_burn_in_allele_fraction} \
            --smoothing-credible-interval-threshold-copy-ratio ${default="2.0" smoothing_threshold_copy_ratio} \
            --smoothing-credible-interval-threshold-allele-fraction ${default="2.0" smoothing_threshold_allele_fraction} \
            --maximum-number-of-smoothing-iterations ${default=10 max_num_smoothing_iterations} \
            --number-of-smoothing-iterations-per-fit ${default=0 num_smoothing_iterations_per_fit} \
            --output ${output_dir_} \
            --output-prefix ${entity_id}

        # We need to create the file even if the above command doesn't so we have something to delocalize
        # If no file is created by the above task then it will copy out an empty file
        touch ${output_dir_}/${entity_id}.hets.normal.tsv
    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: machine_mem + " MB"
        disks: "local-disk " + disk_space_gb + " HDD"
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File het_allelic_counts = "${output_dir_}/${entity_id}.hets.tsv"
        File normal_het_allelic_counts = "${output_dir_}/${entity_id}.hets.normal.tsv"
        File copy_ratio_only_segments = "${output_dir_}/${entity_id}.cr.seg"
        File modeled_segments_begin = "${output_dir_}/${entity_id}.modelBegin.seg"
        File copy_ratio_parameters_begin = "${output_dir_}/${entity_id}.modelBegin.cr.param"
        File allele_fraction_parameters_begin = "${output_dir_}/${entity_id}.modelBegin.af.param"
        File modeled_segments = "${output_dir_}/${entity_id}.modelFinal.seg"
        File copy_ratio_parameters = "${output_dir_}/${entity_id}.modelFinal.cr.param"
        File allele_fraction_parameters = "${output_dir_}/${entity_id}.modelFinal.af.param"
    }
}

task CallCopyRatioSegments {
    String entity_id
    File copy_ratio_segments
    Float? neutral_segment_copy_ratio_threshold
    Float? outlier_neutral_segment_copy_ratio_z_score_threshold
    Float? calling_copy_ratio_z_score_threshold
    File? gatk4_jar_override

    # Runtime parameters
    Int? mem
    String gatk_docker
    Int? preemptible_attempts
    Int disk_space_gb

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1000 else 7000
    Int command_mem = machine_mem - 1000

    command <<<
        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk4_jar_override}

        gatk --java-options "-Xmx${command_mem}m" CallCopyRatioSegments \
            --input ${copy_ratio_segments} \
            --neutral-segment-copy-ratio-threshold ${default="0.1" neutral_segment_copy_ratio_threshold} \
            --outlier-neutral-segment-copy-ratio-z-score-threshold ${default="2.0" outlier_neutral_segment_copy_ratio_z_score_threshold} \
            --calling-copy-ratio-z-score-threshold ${default="2.0" calling_copy_ratio_z_score_threshold} \
            --output ${entity_id}.called.seg
    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: machine_mem + " MB"
        disks: "local-disk " + disk_space_gb + " HDD"
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File called_copy_ratio_segments = "${entity_id}.called.seg"
    }
}

task PlotDenoisedCopyRatios {
    String entity_id
    File standardized_copy_ratios
    File denoised_copy_ratios
    File ref_fasta_dict
    Int? minimum_contig_length
    String? output_dir
    File? gatk4_jar_override

    # Runtime parameters
    Int? mem
    String gatk_docker
    Int? preemptible_attempts
    Int disk_space_gb

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1000 else 7000
    Int command_mem = machine_mem - 1000

    # If optional output_dir not specified, use "out"
    String output_dir_ = select_first([output_dir, "out"])

    command <<<
        set -e
        mkdir ${output_dir_}
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk4_jar_override}

        gatk --java-options "-Xmx${command_mem}m" PlotDenoisedCopyRatios \
            --standardized-copy-ratios ${standardized_copy_ratios} \
            --denoised-copy-ratios ${denoised_copy_ratios} \
            --sequence-dictionary ${ref_fasta_dict} \
            --minimum-contig-length ${default="1000000" minimum_contig_length} \
            --output ${output_dir_} \
            --output-prefix ${entity_id}
    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: machine_mem + " MB"
        disks: "local-disk " + disk_space_gb + " HDD"
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File denoised_copy_ratios_plot = "${output_dir_}/${entity_id}.denoised.png"
        File denoised_copy_ratios_lim_4_plot = "${output_dir_}/${entity_id}.denoisedLimit4.png"
        File standardized_MAD = "${output_dir_}/${entity_id}.standardizedMAD.txt"
        File denoised_MAD = "${output_dir_}/${entity_id}.denoisedMAD.txt"
        File delta_MAD = "${output_dir_}/${entity_id}.deltaMAD.txt"
        File scaled_delta_MAD = "${output_dir_}/${entity_id}.scaledDeltaMAD.txt"
    }
}

task PlotModeledSegments {
    String entity_id
    File denoised_copy_ratios
    File het_allelic_counts
    File modeled_segments
    File ref_fasta_dict
    Int? minimum_contig_length
    String? output_dir
    File? gatk4_jar_override

    # Runtime parameters
    Int? mem
    String gatk_docker
    Int? preemptible_attempts
    Int disk_space_gb

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1000 else 7000
    Int command_mem = machine_mem - 1000

    # If optional output_dir not specified, use "out"
    String output_dir_ = select_first([output_dir, "out"])

    command <<<
        set -e
        mkdir ${output_dir_}
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk4_jar_override}

        gatk --java-options "-Xmx${command_mem}m" PlotModeledSegments \
            --denoised-copy-ratios ${denoised_copy_ratios} \
            --allelic-counts ${het_allelic_counts} \
            --segments ${modeled_segments} \
            --sequence-dictionary ${ref_fasta_dict} \
            --minimum-contig-length ${default="1000000" minimum_contig_length} \
            --output ${output_dir_} \
            --output-prefix ${entity_id}
    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: machine_mem + " MB"
        disks: "local-disk " + disk_space_gb + " HDD"
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File modeled_segments_plot = "${output_dir_}/${entity_id}.modeled.png"
    }
}
