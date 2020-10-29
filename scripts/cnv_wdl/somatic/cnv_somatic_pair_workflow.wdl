# Workflow for running the GATK CNV pipeline on a matched pair. Supports both WGS and WES.
#
# Notes:
#
# - The intervals argument is required for both WGS and WES workflows and accepts formats compatible with the
#   GATK -L argument (see https://gatkforums.broadinstitute.org/gatk/discussion/11009/intervals-and-interval-lists).
#   These intervals will be padded on both sides by the amount specified by padding (default 250)
#   and split into bins of length specified by bin_length (default 1000; specify 0 to skip binning,
#   e.g., for WES).  For WGS, the intervals should simply cover the autosomal chromosomes (sex chromosomes may be
#   included, but care should be taken to 1) avoid creating panels of mixed sex, and 2) denoise case samples only
#   with panels containing only individuals of the same sex as the case samples).
#
# - Intervals can be blacklisted from coverage collection and all downstream steps by using the blacklist_intervals
#   argument, which accepts formats compatible with the GATK -XL argument
#   (see https://gatkforums.broadinstitute.org/gatk/discussion/11009/intervals-and-interval-lists).
#   This may be useful for excluding centromeric regions, etc. from analysis.  Alternatively, these regions may
#   be manually filtered from the final callset.
#
#  A reasonable blacklist for excluded intervals (-XL) can be found at:
#   hg19: gs://gatk-best-practices/somatic-b37/CNV_and_centromere_blacklist.hg19.list
#   hg38: gs://gatk-best-practices/somatic-hg38/CNV_and_centromere_blacklist.hg38liftover.list (untested)
#
# - The sites file (common_sites) should be a Picard or GATK-style interval list.  This is a list of sites
#   of known variation at which allelic counts will be collected for use in modeling minor-allele fractions.
#
# - If you opt to run FuncotateSegments (i.e. set `is_run_funcotator` to `true`), then please also ensure that you have
#       the correct value for `funcotator_ref_version`.  Treat `funcotator_ref_version` as required if
#       `is_run_funcotator` is `true`.  Valid values for `funcotator_ref_version` are `hg38` and `hg19`.
#       The latter includes GRCh37.
#
#
# - Example invocation:
#
#       java -jar cromwell.jar run cnv_somatic_pair_workflow.wdl -i my_parameters.json
#
#############

version 1.0

import "../cnv_common_tasks.wdl" as CNVTasks
import "cnv_somatic_oncotator_workflow.wdl" as CNVOncotator
import "cnv_somatic_funcotate_seg_workflow.wdl" as CNVFuncotateSegments

workflow CNVSomaticPairWorkflow {

    input {
      ##################################
      #### required basic arguments ####
      ##################################
      File common_sites
      File intervals
      File? blacklist_intervals
      File tumor_bam
      File tumor_bam_idx
      File? normal_bam
      File? normal_bam_idx
      File read_count_pon
      File ref_fasta_dict
      File ref_fasta_fai
      File ref_fasta
      String gatk_docker

      ##################################
      #### optional basic arguments ####
      ##################################
       # For running oncotator
      Boolean? is_run_oncotator
       # For running funcotator
      Boolean? is_run_funcotator

      File? gatk4_jar_override
      Int? preemptible_attempts
      # Use as a last resort to increase the disk given to every task in case of ill behaving data
      Int? emergency_extra_disk

      # Required if BAM/CRAM is in a requester pays bucket
      String? gcs_project_for_requester_pays

      ####################################################
      #### optional arguments for PreprocessIntervals ####
      ####################################################
      Int? padding
      Int? bin_length
      Int? mem_gb_for_preprocess_intervals

      ##############################################
      #### optional arguments for CollectCounts ####
      ##############################################
      String? collect_counts_format
      Int? mem_gb_for_collect_counts

      #####################################################
      #### optional arguments for CollectAllelicCounts ####
      #####################################################
      String? minimum_base_quality
      Int? mem_gb_for_collect_allelic_counts

      ##################################################
      #### optional arguments for DenoiseReadCounts ####
      ##################################################
      Int? number_of_eigensamples
      Int? mem_gb_for_denoise_read_counts

      ##############################################
      #### optional arguments for ModelSegments ####
      ##############################################
      Int? max_num_segments_per_chromosome
      Int? min_total_allele_count
      Int? min_total_allele_count_normal
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
      Int? mem_gb_for_model_segments

      ######################################################
      #### optional arguments for CallCopyRatioSegments ####
      ######################################################
      Float? neutral_segment_copy_ratio_lower_bound
      Float? neutral_segment_copy_ratio_upper_bound
      Float? outlier_neutral_segment_copy_ratio_z_score_threshold
      Float? calling_copy_ratio_z_score_threshold
      Int? mem_gb_for_call_copy_ratio_segments

      #########################################
      #### optional arguments for plotting ####
      #########################################
      Int? minimum_contig_length
      # If maximum_copy_ratio = Infinity, the maximum copy ratio will be automatically determined
      String? maximum_copy_ratio
      Float? point_size_copy_ratio
      Float? point_size_allele_fraction
      Int? mem_gb_for_plotting

      ##########################################
      #### optional arguments for Oncotator ####
      ##########################################
      String? additional_args_for_oncotator
      String? oncotator_docker
      Int? mem_gb_for_oncotator
      Int? boot_disk_space_gb_for_oncotator

      ##################################################
      #### optional arguments for FuncotateSegments ####
      ##################################################
      String? additional_args_for_funcotator
      String? funcotator_ref_version
      Int? mem_gb_for_funcotator
      File? funcotator_transcript_selection_list
      File? funcotator_data_sources_tar_gz
      String? funcotator_transcript_selection_mode
      Array[String]? funcotator_annotation_defaults
      Array[String]? funcotator_annotation_overrides
      Array[String]? funcotator_excluded_fields
      Boolean? funcotator_is_removing_untared_datasources
      Int? funcotator_disk_space_gb
      Boolean? funcotator_use_ssd
      Int? funcotator_cpu
    }

    Int ref_size = ceil(size(ref_fasta, "GB") + size(ref_fasta_dict, "GB") + size(ref_fasta_fai, "GB"))
    Int read_count_pon_size = ceil(size(read_count_pon, "GB"))
    Int tumor_bam_size = ceil(size(tumor_bam, "GB") + size(tumor_bam_idx, "GB"))
    Int normal_bam_size = if defined(normal_bam) then ceil(size(normal_bam, "GB") + size(normal_bam_idx, "GB")) else 0

    Int gatk4_override_size = if defined(gatk4_jar_override) then ceil(size(gatk4_jar_override, "GB")) else 0
    # This is added to every task as padding, should increase if systematically you need more disk for every call
    Int disk_pad = 20 + ceil(size(intervals, "GB")) + ceil(size(common_sites, "GB")) + gatk4_override_size + select_first([emergency_extra_disk, 0])

    File final_normal_bam = select_first([normal_bam, "null"])
    File final_normal_bam_idx = select_first([normal_bam_idx, "null"])

    Int preprocess_intervals_disk = ref_size + disk_pad
    call CNVTasks.PreprocessIntervals {
        input:
            intervals = intervals,
            blacklist_intervals = blacklist_intervals,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            ref_fasta_dict = ref_fasta_dict,
            padding = padding,
            bin_length = bin_length,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            mem_gb = mem_gb_for_preprocess_intervals,
            disk_space_gb = preprocess_intervals_disk,
            preemptible_attempts = preemptible_attempts
    }

    Int collect_counts_tumor_disk = tumor_bam_size + ceil(size(PreprocessIntervals.preprocessed_intervals, "GB")) + disk_pad
    call CNVTasks.CollectCounts as CollectCountsTumor {
        input:
            intervals = PreprocessIntervals.preprocessed_intervals,
            bam = tumor_bam,
            bam_idx = tumor_bam_idx,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            ref_fasta_dict = ref_fasta_dict,
            format = collect_counts_format,
            enable_indexing = false,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            mem_gb = mem_gb_for_collect_counts,
            disk_space_gb = collect_counts_tumor_disk,
            preemptible_attempts = preemptible_attempts,
            gcs_project_for_requester_pays = gcs_project_for_requester_pays
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
            minimum_base_quality =  minimum_base_quality,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            mem_gb = mem_gb_for_collect_allelic_counts,
            disk_space_gb = collect_allelic_counts_tumor_disk,
            preemptible_attempts = preemptible_attempts,
            gcs_project_for_requester_pays = gcs_project_for_requester_pays
    }

    Int denoise_read_counts_tumor_disk = read_count_pon_size + ceil(size(CollectCountsTumor.counts, "GB")) + disk_pad
    call DenoiseReadCounts as DenoiseReadCountsTumor {
        input:
            entity_id = CollectCountsTumor.entity_id,
            read_counts = CollectCountsTumor.counts,
            read_count_pon = read_count_pon,
            number_of_eigensamples = number_of_eigensamples,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            mem_gb = mem_gb_for_denoise_read_counts,
            disk_space_gb = denoise_read_counts_tumor_disk,
            preemptible_attempts = preemptible_attempts
    }

    Int model_segments_normal_portion = if defined(normal_bam) then ceil(size(CollectAllelicCountsNormal.allelic_counts, "GB")) else 0
    Int model_segments_tumor_disk = ceil(size(DenoiseReadCountsTumor.denoised_copy_ratios, "GB")) + ceil(size(CollectAllelicCountsTumor.allelic_counts, "GB")) + model_segments_normal_portion + disk_pad
    call ModelSegments as ModelSegmentsTumor {
        input:
            entity_id = CollectCountsTumor.entity_id,
            denoised_copy_ratios = DenoiseReadCountsTumor.denoised_copy_ratios,
            allelic_counts = CollectAllelicCountsTumor.allelic_counts,
            normal_allelic_counts = CollectAllelicCountsNormal.allelic_counts,
            max_num_segments_per_chromosome = max_num_segments_per_chromosome,
            min_total_allele_count = min_total_allele_count,
            min_total_allele_count_normal = min_total_allele_count_normal,
            genotyping_homozygous_log_ratio_threshold = genotyping_homozygous_log_ratio_threshold,
            genotyping_base_error_rate = genotyping_base_error_rate,
            kernel_variance_copy_ratio = kernel_variance_copy_ratio,
            kernel_variance_allele_fraction = kernel_variance_allele_fraction,
            kernel_scaling_allele_fraction = kernel_scaling_allele_fraction,
            kernel_approximation_dimension = kernel_approximation_dimension,
            window_sizes = window_sizes,
            num_changepoints_penalty_factor = num_changepoints_penalty_factor,
            minor_allele_fraction_prior_alpha = minor_allele_fraction_prior_alpha,
            num_samples_copy_ratio = num_samples_copy_ratio,
            num_burn_in_copy_ratio = num_burn_in_copy_ratio,
            num_samples_allele_fraction = num_samples_allele_fraction,
            num_burn_in_allele_fraction = num_burn_in_allele_fraction,
            smoothing_threshold_copy_ratio = smoothing_threshold_copy_ratio,
            smoothing_threshold_allele_fraction = smoothing_threshold_allele_fraction,
            max_num_smoothing_iterations = max_num_smoothing_iterations,
            num_smoothing_iterations_per_fit = num_smoothing_iterations_per_fit,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            mem_gb = mem_gb_for_model_segments,
            disk_space_gb = model_segments_tumor_disk,
            preemptible_attempts = preemptible_attempts
    }

    Int copy_ratio_segments_tumor_disk = ceil(size(DenoiseReadCountsTumor.denoised_copy_ratios, "GB")) + ceil(size(ModelSegmentsTumor.copy_ratio_only_segments, "GB")) + disk_pad
    call CallCopyRatioSegments as CallCopyRatioSegmentsTumor {
        input:
            entity_id = CollectCountsTumor.entity_id,
            copy_ratio_segments = ModelSegmentsTumor.copy_ratio_only_segments,
            neutral_segment_copy_ratio_lower_bound = neutral_segment_copy_ratio_lower_bound,
            neutral_segment_copy_ratio_upper_bound = neutral_segment_copy_ratio_upper_bound,
            outlier_neutral_segment_copy_ratio_z_score_threshold = outlier_neutral_segment_copy_ratio_z_score_threshold,
            calling_copy_ratio_z_score_threshold = calling_copy_ratio_z_score_threshold,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            mem_gb = mem_gb_for_call_copy_ratio_segments,
            disk_space_gb = copy_ratio_segments_tumor_disk,
            preemptible_attempts = preemptible_attempts
    }

    # The F=files from other tasks are small enough to just combine into one disk variable and pass to the tumor plotting tasks
    Int plot_tumor_disk = ref_size + ceil(size(DenoiseReadCountsTumor.standardized_copy_ratios, "GB")) + ceil(size(DenoiseReadCountsTumor.denoised_copy_ratios, "GB")) + ceil(size(ModelSegmentsTumor.het_allelic_counts, "GB")) + ceil(size(ModelSegmentsTumor.modeled_segments, "GB")) + disk_pad
    call PlotDenoisedCopyRatios as PlotDenoisedCopyRatiosTumor {
        input:
            entity_id = CollectCountsTumor.entity_id,
            standardized_copy_ratios = DenoiseReadCountsTumor.standardized_copy_ratios,
            denoised_copy_ratios = DenoiseReadCountsTumor.denoised_copy_ratios,
            ref_fasta_dict = ref_fasta_dict,
            minimum_contig_length = minimum_contig_length,
            maximum_copy_ratio = maximum_copy_ratio,
            point_size_copy_ratio = point_size_copy_ratio,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            mem_gb = mem_gb_for_plotting,
            disk_space_gb = plot_tumor_disk,
            preemptible_attempts = preemptible_attempts
    }

    call PlotModeledSegments as PlotModeledSegmentsTumor {
        input:
            entity_id = CollectCountsTumor.entity_id,
            denoised_copy_ratios = DenoiseReadCountsTumor.denoised_copy_ratios,
            het_allelic_counts = ModelSegmentsTumor.het_allelic_counts,
            modeled_segments = ModelSegmentsTumor.modeled_segments,
            ref_fasta_dict = ref_fasta_dict,
            minimum_contig_length = minimum_contig_length,
            maximum_copy_ratio = maximum_copy_ratio,
            point_size_copy_ratio = point_size_copy_ratio,
            point_size_allele_fraction = point_size_allele_fraction,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            mem_gb = mem_gb_for_plotting,
            disk_space_gb = plot_tumor_disk,
            preemptible_attempts = preemptible_attempts
    }

    Int collect_counts_normal_disk = normal_bam_size + ceil(size(PreprocessIntervals.preprocessed_intervals, "GB")) + disk_pad
    if (defined(normal_bam)) {
        call CNVTasks.CollectCounts as CollectCountsNormal {
            input:
                intervals = PreprocessIntervals.preprocessed_intervals,
                bam = final_normal_bam,
                bam_idx = final_normal_bam_idx,
                ref_fasta = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,
                ref_fasta_dict = ref_fasta_dict,
                format = collect_counts_format,
                enable_indexing = false,
                gatk4_jar_override = gatk4_jar_override,
                gatk_docker = gatk_docker,
                mem_gb = mem_gb_for_collect_counts,
                disk_space_gb = collect_counts_normal_disk,
                preemptible_attempts = preemptible_attempts,
                gcs_project_for_requester_pays = gcs_project_for_requester_pays
        }

        Int collect_allelic_counts_normal_disk = normal_bam_size + ref_size + disk_pad
        call CNVTasks.CollectAllelicCounts as CollectAllelicCountsNormal {
            input:
                common_sites = common_sites,
                bam = final_normal_bam,
                bam_idx = final_normal_bam_idx,
                ref_fasta = ref_fasta,
                ref_fasta_dict = ref_fasta_dict,
                ref_fasta_fai = ref_fasta_fai,
                minimum_base_quality =  minimum_base_quality,
                gatk4_jar_override = gatk4_jar_override,
                gatk_docker = gatk_docker,
                mem_gb = mem_gb_for_collect_allelic_counts,
                disk_space_gb = collect_allelic_counts_normal_disk,
                preemptible_attempts = preemptible_attempts,
                gcs_project_for_requester_pays = gcs_project_for_requester_pays
        }

        Int denoise_read_counts_normal_disk = read_count_pon_size + ceil(size(CollectCountsNormal.counts, "GB")) + disk_pad
        call DenoiseReadCounts as DenoiseReadCountsNormal {
            input:
                entity_id = CollectCountsNormal.entity_id,
                read_counts = CollectCountsNormal.counts,
                read_count_pon = read_count_pon,
                number_of_eigensamples = number_of_eigensamples,
                gatk4_jar_override = gatk4_jar_override,
                gatk_docker = gatk_docker,
                mem_gb = mem_gb_for_denoise_read_counts,
                disk_space_gb = denoise_read_counts_normal_disk,
                preemptible_attempts = preemptible_attempts
        }

        Int model_segments_normal_disk = ceil(size(DenoiseReadCountsNormal.denoised_copy_ratios, "GB")) + ceil(size(CollectAllelicCountsNormal.allelic_counts, "GB")) + disk_pad
        call ModelSegments as ModelSegmentsNormal {
            input:
                entity_id = CollectCountsNormal.entity_id,
                denoised_copy_ratios = DenoiseReadCountsNormal.denoised_copy_ratios,
                allelic_counts = CollectAllelicCountsNormal.allelic_counts,
                max_num_segments_per_chromosome = max_num_segments_per_chromosome,
                min_total_allele_count = min_total_allele_count_normal,
                genotyping_homozygous_log_ratio_threshold = genotyping_homozygous_log_ratio_threshold,
                genotyping_base_error_rate = genotyping_base_error_rate,
                kernel_variance_copy_ratio = kernel_variance_copy_ratio,
                kernel_variance_allele_fraction = kernel_variance_allele_fraction,
                kernel_scaling_allele_fraction = kernel_scaling_allele_fraction,
                kernel_approximation_dimension = kernel_approximation_dimension,
                window_sizes = window_sizes,
                num_changepoints_penalty_factor = num_changepoints_penalty_factor,
                minor_allele_fraction_prior_alpha = minor_allele_fraction_prior_alpha,
                num_samples_copy_ratio = num_samples_copy_ratio,
                num_burn_in_copy_ratio = num_burn_in_copy_ratio,
                num_samples_allele_fraction = num_samples_allele_fraction,
                num_burn_in_allele_fraction = num_burn_in_allele_fraction,
                smoothing_threshold_copy_ratio = smoothing_threshold_copy_ratio,
                smoothing_threshold_allele_fraction = smoothing_threshold_allele_fraction,
                max_num_smoothing_iterations = max_num_smoothing_iterations,
                num_smoothing_iterations_per_fit = num_smoothing_iterations_per_fit,
                gatk4_jar_override = gatk4_jar_override,
                gatk_docker = gatk_docker,
                mem_gb = mem_gb_for_model_segments,
                disk_space_gb = model_segments_normal_disk,
                preemptible_attempts = preemptible_attempts
        }

        Int copy_ratio_segments_normal_disk = ceil(size(DenoiseReadCountsNormal.denoised_copy_ratios, "GB")) + ceil(size(ModelSegmentsNormal.copy_ratio_only_segments, "GB")) + disk_pad
        call CallCopyRatioSegments as CallCopyRatioSegmentsNormal {
            input:
                entity_id = CollectCountsNormal.entity_id,
                copy_ratio_segments = ModelSegmentsNormal.copy_ratio_only_segments,
                neutral_segment_copy_ratio_lower_bound = neutral_segment_copy_ratio_lower_bound,
                neutral_segment_copy_ratio_upper_bound = neutral_segment_copy_ratio_upper_bound,
                outlier_neutral_segment_copy_ratio_z_score_threshold = outlier_neutral_segment_copy_ratio_z_score_threshold,
                calling_copy_ratio_z_score_threshold = calling_copy_ratio_z_score_threshold,
                gatk4_jar_override = gatk4_jar_override,
                gatk_docker = gatk_docker,
                mem_gb = mem_gb_for_call_copy_ratio_segments,
                disk_space_gb = copy_ratio_segments_normal_disk,
                preemptible_attempts = preemptible_attempts
        }

        # The files from other tasks are small enough to just combine into one disk variable and pass to the normal plotting tasks
        Int plot_normal_disk = ref_size + ceil(size(DenoiseReadCountsNormal.standardized_copy_ratios, "GB")) + ceil(size(DenoiseReadCountsNormal.denoised_copy_ratios, "GB")) + ceil(size(ModelSegmentsNormal.het_allelic_counts, "GB")) + ceil(size(ModelSegmentsNormal.modeled_segments, "GB")) + disk_pad
        call PlotDenoisedCopyRatios as PlotDenoisedCopyRatiosNormal {
            input:
                entity_id = CollectCountsNormal.entity_id,
                standardized_copy_ratios = DenoiseReadCountsNormal.standardized_copy_ratios,
                denoised_copy_ratios = DenoiseReadCountsNormal.denoised_copy_ratios,
                ref_fasta_dict = ref_fasta_dict,
                minimum_contig_length = minimum_contig_length,
                maximum_copy_ratio = maximum_copy_ratio,
                point_size_copy_ratio = point_size_copy_ratio,
                gatk4_jar_override = gatk4_jar_override,
                gatk_docker = gatk_docker,
                mem_gb = mem_gb_for_plotting,
                disk_space_gb = plot_normal_disk,
                preemptible_attempts = preemptible_attempts
        }

        call PlotModeledSegments as PlotModeledSegmentsNormal {
            input:
                entity_id = CollectCountsNormal.entity_id,
                denoised_copy_ratios = DenoiseReadCountsNormal.denoised_copy_ratios,
                het_allelic_counts = ModelSegmentsNormal.het_allelic_counts,
                modeled_segments = ModelSegmentsNormal.modeled_segments,
                ref_fasta_dict = ref_fasta_dict,
                minimum_contig_length = minimum_contig_length,
                maximum_copy_ratio = maximum_copy_ratio,
                point_size_copy_ratio = point_size_copy_ratio,
                point_size_allele_fraction = point_size_allele_fraction,
                gatk4_jar_override = gatk4_jar_override,
                gatk_docker = gatk_docker,
                mem_gb = mem_gb_for_plotting,
                disk_space_gb = plot_normal_disk,
                preemptible_attempts = preemptible_attempts
        }
    }

    if (select_first([is_run_oncotator, false])) {
        call CNVOncotator.CNVOncotatorWorkflow as CNVOncotatorWorkflow {
            input:
                 called_file = CallCopyRatioSegmentsTumor.called_copy_ratio_segments,
                 additional_args = additional_args_for_oncotator,
                 oncotator_docker = oncotator_docker,
                 mem_gb_for_oncotator = mem_gb_for_oncotator,
                 boot_disk_space_gb_for_oncotator = boot_disk_space_gb_for_oncotator,
                 preemptible_attempts = preemptible_attempts
        }
    }
    if (select_first([is_run_funcotator, false])) {
        call CNVFuncotateSegments.CNVFuncotateSegmentsWorkflow as CNVFuncotateSegmentsWorkflow {
            input:
                 input_seg_file = CallCopyRatioSegmentsTumor.called_copy_ratio_segments,
                 funcotator_ref_version = select_first([funcotator_ref_version, "hg19"]),
                 extra_args = additional_args_for_funcotator,
                 ref_fasta = ref_fasta,
                 ref_fasta_fai = ref_fasta_fai,
                 ref_fasta_dict = ref_fasta_dict,
                 transcript_selection_list = funcotator_transcript_selection_list,
                 funcotator_data_sources_tar_gz = funcotator_data_sources_tar_gz,
                 gatk4_jar_override = gatk4_jar_override,
                 gatk_docker = gatk_docker,
                 mem_gb = mem_gb_for_funcotator,
                 preemptible_attempts = preemptible_attempts,
                 transcript_selection_mode = funcotator_transcript_selection_mode,
                 annotation_defaults = funcotator_annotation_defaults,
                 annotation_overrides = funcotator_annotation_overrides,
                 funcotator_excluded_fields = funcotator_excluded_fields,
                 is_removing_untared_datasources = funcotator_is_removing_untared_datasources,
                 disk_space_gb = funcotator_disk_space_gb,
                 use_ssd = funcotator_use_ssd,
                 cpu = funcotator_cpu
        }
    }

    output {
        File preprocessed_intervals = PreprocessIntervals.preprocessed_intervals

        File read_counts_entity_id_tumor = CollectCountsTumor.entity_id
        File read_counts_tumor = CollectCountsTumor.counts
        File allelic_counts_entity_id_tumor = CollectAllelicCountsTumor.entity_id
        File allelic_counts_tumor = CollectAllelicCountsTumor.allelic_counts
        File denoised_copy_ratios_tumor = DenoiseReadCountsTumor.denoised_copy_ratios
        File standardized_copy_ratios_tumor = DenoiseReadCountsTumor.standardized_copy_ratios
        File het_allelic_counts_tumor = ModelSegmentsTumor.het_allelic_counts
        File normal_het_allelic_counts_tumor = ModelSegmentsTumor.normal_het_allelic_counts
        File copy_ratio_only_segments_tumor = ModelSegmentsTumor.copy_ratio_only_segments
        File copy_ratio_legacy_segments_tumor = ModelSegmentsTumor.copy_ratio_legacy_segments
        File allele_fraction_legacy_segments_tumor = ModelSegmentsTumor.allele_fraction_legacy_segments
        File modeled_segments_begin_tumor = ModelSegmentsTumor.modeled_segments_begin
        File copy_ratio_parameters_begin_tumor = ModelSegmentsTumor.copy_ratio_parameters_begin
        File allele_fraction_parameters_begin_tumor = ModelSegmentsTumor.allele_fraction_parameters_begin
        File modeled_segments_tumor = ModelSegmentsTumor.modeled_segments
        File copy_ratio_parameters_tumor = ModelSegmentsTumor.copy_ratio_parameters
        File allele_fraction_parameters_tumor = ModelSegmentsTumor.allele_fraction_parameters
        File called_copy_ratio_segments_tumor = CallCopyRatioSegmentsTumor.called_copy_ratio_segments
        File called_copy_ratio_legacy_segments_tumor = CallCopyRatioSegmentsTumor.called_copy_ratio_legacy_segments
        File denoised_copy_ratios_plot_tumor = PlotDenoisedCopyRatiosTumor.denoised_copy_ratios_plot
        File standardized_MAD_tumor = PlotDenoisedCopyRatiosTumor.standardized_MAD
        Float standardized_MAD_value_tumor = PlotDenoisedCopyRatiosTumor.standardized_MAD_value
        File denoised_MAD_tumor = PlotDenoisedCopyRatiosTumor.denoised_MAD
        Float denoised_MAD_value_tumor = PlotDenoisedCopyRatiosTumor.denoised_MAD_value
        File delta_MAD_tumor = PlotDenoisedCopyRatiosTumor.delta_MAD
        Float delta_MAD_value_tumor = PlotDenoisedCopyRatiosTumor.delta_MAD_value
        File scaled_delta_MAD_tumor = PlotDenoisedCopyRatiosTumor.scaled_delta_MAD
        Float scaled_delta_MAD_value_tumor = PlotDenoisedCopyRatiosTumor.scaled_delta_MAD_value
        File modeled_segments_plot_tumor = PlotModeledSegmentsTumor.modeled_segments_plot

        File? read_counts_entity_id_normal = CollectCountsNormal.entity_id
        File? read_counts_normal = CollectCountsNormal.counts
        File? allelic_counts_entity_id_normal = CollectAllelicCountsNormal.entity_id
        File? allelic_counts_normal = CollectAllelicCountsNormal.allelic_counts
        File? denoised_copy_ratios_normal = DenoiseReadCountsNormal.denoised_copy_ratios
        File? standardized_copy_ratios_normal = DenoiseReadCountsNormal.standardized_copy_ratios
        File? het_allelic_counts_normal = ModelSegmentsNormal.het_allelic_counts
        File? normal_het_allelic_counts_normal = ModelSegmentsNormal.normal_het_allelic_counts
        File? copy_ratio_only_segments_normal = ModelSegmentsNormal.copy_ratio_only_segments
        File? copy_ratio_legacy_segments_normal = ModelSegmentsNormal.copy_ratio_legacy_segments
        File? allele_fraction_legacy_segments_normal = ModelSegmentsNormal.allele_fraction_legacy_segments
        File? modeled_segments_begin_normal = ModelSegmentsNormal.modeled_segments_begin
        File? copy_ratio_parameters_begin_normal = ModelSegmentsNormal.copy_ratio_parameters_begin
        File? allele_fraction_parameters_begin_normal = ModelSegmentsNormal.allele_fraction_parameters_begin
        File? modeled_segments_normal = ModelSegmentsNormal.modeled_segments
        File? copy_ratio_parameters_normal = ModelSegmentsNormal.copy_ratio_parameters
        File? allele_fraction_parameters_normal = ModelSegmentsNormal.allele_fraction_parameters
        File? called_copy_ratio_segments_normal = CallCopyRatioSegmentsNormal.called_copy_ratio_segments
        File? called_copy_ratio_legacy_segments_normal = CallCopyRatioSegmentsNormal.called_copy_ratio_legacy_segments
        File? denoised_copy_ratios_plot_normal = PlotDenoisedCopyRatiosNormal.denoised_copy_ratios_plot
        File? standardized_MAD_normal = PlotDenoisedCopyRatiosNormal.standardized_MAD
        Float? standardized_MAD_value_normal = PlotDenoisedCopyRatiosNormal.standardized_MAD_value
        File? denoised_MAD_normal = PlotDenoisedCopyRatiosNormal.denoised_MAD
        Float? denoised_MAD_value_normal = PlotDenoisedCopyRatiosNormal.denoised_MAD_value
        File? delta_MAD_normal = PlotDenoisedCopyRatiosNormal.delta_MAD
        Float? delta_MAD_value_normal = PlotDenoisedCopyRatiosNormal.delta_MAD_value
        File? scaled_delta_MAD_normal = PlotDenoisedCopyRatiosNormal.scaled_delta_MAD
        Float? scaled_delta_MAD_value_normal = PlotDenoisedCopyRatiosNormal.scaled_delta_MAD_value
        File? modeled_segments_plot_normal = PlotModeledSegmentsNormal.modeled_segments_plot

        File oncotated_called_file_tumor = select_first([CNVOncotatorWorkflow.oncotated_called_file, "null"])
        File oncotated_called_gene_list_file_tumor = select_first([CNVOncotatorWorkflow.oncotated_called_gene_list_file, "null"])
        File funcotated_called_file_tumor = select_first([CNVFuncotateSegmentsWorkflow.funcotated_seg_simple_tsv, "null"])
        File funcotated_called_gene_list_file_tumor = select_first([CNVFuncotateSegmentsWorkflow.funcotated_gene_list_tsv, "null"])
    }
}

task DenoiseReadCounts {
    input {
      String entity_id
      File read_counts
      File read_count_pon
      Int? number_of_eigensamples #use all eigensamples in panel by default
      File? gatk4_jar_override

      # Runtime parameters
      String gatk_docker
      Int? mem_gb
      Int? disk_space_gb
      Boolean use_ssd = false
      Int? cpu
      Int? preemptible_attempts
    }

    Int machine_mem_mb = select_first([mem_gb, 13]) * 1000
    Int command_mem_mb = machine_mem_mb - 1000

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}

        gatk --java-options "-Xmx~{command_mem_mb}m" DenoiseReadCounts \
            --input ~{read_counts} \
            --count-panel-of-normals ~{read_count_pon} \
            ~{"--number-of-eigensamples " + number_of_eigensamples} \
            --standardized-copy-ratios ~{entity_id}.standardizedCR.tsv \
            --denoised-copy-ratios ~{entity_id}.denoisedCR.tsv
    >>>

    runtime {
        docker: "~{gatk_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + disk_space_gb + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File standardized_copy_ratios = "~{entity_id}.standardizedCR.tsv"
        File denoised_copy_ratios = "~{entity_id}.denoisedCR.tsv"
    }
}

task ModelSegments {
    input {
      String entity_id
      File denoised_copy_ratios
      File allelic_counts
      File? normal_allelic_counts
      Int? max_num_segments_per_chromosome
      Int? min_total_allele_count
      Int? min_total_allele_count_normal
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
      String gatk_docker
      Int? mem_gb
      Int? disk_space_gb
      Boolean use_ssd = false
      Int? cpu
      Int? preemptible_attempts
    }

    Int machine_mem_mb = select_first([mem_gb, 13]) * 1000
    # ModelSegments seems to need at least 3GB of overhead to run
    Int command_mem_mb = machine_mem_mb - 3000

    # If optional output_dir not specified, use "out"
    String output_dir_ = select_first([output_dir, "out"])

    # default values are min_total_allele_count_ = 0 in matched-normal mode
    #                                            = 30 in case-only mode
    Int default_min_total_allele_count = if defined(normal_allelic_counts) then 0 else 30
    Int min_total_allele_count_ = select_first([min_total_allele_count, default_min_total_allele_count])

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}

        gatk --java-options "-Xmx~{command_mem_mb}m" ModelSegments \
            --denoised-copy-ratios ~{denoised_copy_ratios} \
            --allelic-counts ~{allelic_counts} \
            ~{"--normal-allelic-counts " + normal_allelic_counts} \
            --minimum-total-allele-count-case ~{min_total_allele_count_} \
            --minimum-total-allele-count-normal ~{default="30" min_total_allele_count_normal} \
            --genotyping-homozygous-log-ratio-threshold ~{default="-10.0" genotyping_homozygous_log_ratio_threshold} \
            --genotyping-base-error-rate ~{default="0.05" genotyping_base_error_rate} \
            --maximum-number-of-segments-per-chromosome ~{default="1000" max_num_segments_per_chromosome} \
            --kernel-variance-copy-ratio ~{default="0.0" kernel_variance_copy_ratio} \
            --kernel-variance-allele-fraction ~{default="0.025" kernel_variance_allele_fraction} \
            --kernel-scaling-allele-fraction ~{default="1.0" kernel_scaling_allele_fraction} \
            --kernel-approximation-dimension ~{default="100" kernel_approximation_dimension} \
            --window-size ~{sep=" --window-size " window_sizes} \
            --number-of-changepoints-penalty-factor ~{default="1.0" num_changepoints_penalty_factor} \
            --minor-allele-fraction-prior-alpha ~{default="25.0" minor_allele_fraction_prior_alpha} \
            --number-of-samples-copy-ratio ~{default="100" num_samples_copy_ratio} \
            --number-of-burn-in-samples-copy-ratio ~{default="50" num_burn_in_copy_ratio} \
            --number-of-samples-allele-fraction ~{default="100" num_samples_allele_fraction} \
            --number-of-burn-in-samples-allele-fraction ~{default="50" num_burn_in_allele_fraction} \
            --smoothing-credible-interval-threshold-copy-ratio ~{default="2.0" smoothing_threshold_copy_ratio} \
            --smoothing-credible-interval-threshold-allele-fraction ~{default="2.0" smoothing_threshold_allele_fraction} \
            --maximum-number-of-smoothing-iterations ~{default="10" max_num_smoothing_iterations} \
            --number-of-smoothing-iterations-per-fit ~{default="0" num_smoothing_iterations_per_fit} \
            --output ~{output_dir_} \
            --output-prefix ~{entity_id}

        # We need to create the file even if the above command doesn't so we have something to delocalize
        # If no file is created by the above task then it will copy out an empty file
        touch ~{output_dir_}/~{entity_id}.hets.normal.tsv
    >>>

    runtime {
        docker: "~{gatk_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + disk_space_gb + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File het_allelic_counts = "~{output_dir_}/~{entity_id}.hets.tsv"
        File normal_het_allelic_counts = "~{output_dir_}/~{entity_id}.hets.normal.tsv"
        File copy_ratio_only_segments = "~{output_dir_}/~{entity_id}.cr.seg"
        File copy_ratio_legacy_segments = "~{output_dir_}/~{entity_id}.cr.igv.seg"
        File allele_fraction_legacy_segments = "~{output_dir_}/~{entity_id}.af.igv.seg"
        File modeled_segments_begin = "~{output_dir_}/~{entity_id}.modelBegin.seg"
        File copy_ratio_parameters_begin = "~{output_dir_}/~{entity_id}.modelBegin.cr.param"
        File allele_fraction_parameters_begin = "~{output_dir_}/~{entity_id}.modelBegin.af.param"
        File modeled_segments = "~{output_dir_}/~{entity_id}.modelFinal.seg"
        File copy_ratio_parameters = "~{output_dir_}/~{entity_id}.modelFinal.cr.param"
        File allele_fraction_parameters = "~{output_dir_}/~{entity_id}.modelFinal.af.param"
    }
}

task CallCopyRatioSegments {
    input {
      String entity_id
      File copy_ratio_segments
      Float? neutral_segment_copy_ratio_lower_bound
      Float? neutral_segment_copy_ratio_upper_bound
      Float? outlier_neutral_segment_copy_ratio_z_score_threshold
      Float? calling_copy_ratio_z_score_threshold
      File? gatk4_jar_override

      # Runtime parameters
      String gatk_docker
      Int? mem_gb
      Int? disk_space_gb
      Boolean use_ssd = false
      Int? cpu
      Int? preemptible_attempts
    }

    Int machine_mem_mb = select_first([mem_gb, 7]) * 1000
    Int command_mem_mb = machine_mem_mb - 1000

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}

        gatk --java-options "-Xmx~{command_mem_mb}m" CallCopyRatioSegments \
            --input ~{copy_ratio_segments} \
            --neutral-segment-copy-ratio-lower-bound ~{default="0.9" neutral_segment_copy_ratio_lower_bound} \
            --neutral-segment-copy-ratio-upper-bound ~{default="1.1" neutral_segment_copy_ratio_upper_bound} \
            --outlier-neutral-segment-copy-ratio-z-score-threshold ~{default="2.0" outlier_neutral_segment_copy_ratio_z_score_threshold} \
            --calling-copy-ratio-z-score-threshold ~{default="2.0" calling_copy_ratio_z_score_threshold} \
            --output ~{entity_id}.called.seg
    >>>

    runtime {
        docker: "~{gatk_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + disk_space_gb + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File called_copy_ratio_segments = "~{entity_id}.called.seg"
        File called_copy_ratio_legacy_segments = "~{entity_id}.called.igv.seg"
    }
}

task PlotDenoisedCopyRatios {
    input {
      String entity_id
      File standardized_copy_ratios
      File denoised_copy_ratios
      File ref_fasta_dict
      Int? minimum_contig_length
      String? maximum_copy_ratio
      Float? point_size_copy_ratio
      String? output_dir
      File? gatk4_jar_override

      # Runtime parameters
      String gatk_docker
      Int? mem_gb
      Int? disk_space_gb
      Boolean use_ssd = false
      Int? cpu
      Int? preemptible_attempts
    }

    Int machine_mem_mb = select_first([mem_gb, 7]) * 1000
    Int command_mem_mb = machine_mem_mb - 1000

    # If optional output_dir not specified, use "out"
    String output_dir_ = select_first([output_dir, "out"])

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}

        gatk --java-options "-Xmx~{command_mem_mb}m" PlotDenoisedCopyRatios \
            --standardized-copy-ratios ~{standardized_copy_ratios} \
            --denoised-copy-ratios ~{denoised_copy_ratios} \
            --sequence-dictionary ~{ref_fasta_dict} \
            --minimum-contig-length ~{default="1000000" minimum_contig_length} \
            --maximum-copy-ratio ~{default="4.0" maximum_copy_ratio} \
            --point-size-copy-ratio ~{default="0.2" point_size_copy_ratio} \
            --output ~{output_dir_} \
            --output-prefix ~{entity_id}
    >>>

    runtime {
        docker: "~{gatk_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + disk_space_gb + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File denoised_copy_ratios_plot = "~{output_dir_}/~{entity_id}.denoised.png"
        File standardized_MAD = "~{output_dir_}/~{entity_id}.standardizedMAD.txt"
        Float standardized_MAD_value = read_float(standardized_MAD)
        File denoised_MAD = "~{output_dir_}/~{entity_id}.denoisedMAD.txt"
        Float denoised_MAD_value = read_float(denoised_MAD)
        File delta_MAD = "~{output_dir_}/~{entity_id}.deltaMAD.txt"
        Float delta_MAD_value = read_float(delta_MAD)
        File scaled_delta_MAD = "~{output_dir_}/~{entity_id}.scaledDeltaMAD.txt"
        Float scaled_delta_MAD_value = read_float(scaled_delta_MAD)
    }
}

task PlotModeledSegments {
    input {
      String entity_id
      File denoised_copy_ratios
      File het_allelic_counts
      File modeled_segments
      File ref_fasta_dict
      Int? minimum_contig_length
      String? maximum_copy_ratio
      Float? point_size_copy_ratio
      Float? point_size_allele_fraction
      String? output_dir
      File? gatk4_jar_override

      # Runtime parameters
      String gatk_docker
      Int? mem_gb
      Int? disk_space_gb
      Boolean use_ssd = false
      Int? cpu
      Int? preemptible_attempts
    }

    Int machine_mem_mb = select_first([mem_gb, 7]) * 1000
    Int command_mem_mb = machine_mem_mb - 1000

    # If optional output_dir not specified, use "out"
    String output_dir_ = select_first([output_dir, "out"])

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}

        gatk --java-options "-Xmx~{command_mem_mb}m" PlotModeledSegments \
            --denoised-copy-ratios ~{denoised_copy_ratios} \
            --allelic-counts ~{het_allelic_counts} \
            --segments ~{modeled_segments} \
            --sequence-dictionary ~{ref_fasta_dict} \
            --minimum-contig-length ~{default="1000000" minimum_contig_length} \
            --maximum-copy-ratio ~{default="4.0" maximum_copy_ratio} \
            --point-size-copy-ratio ~{default="0.2" point_size_copy_ratio} \
            --point-size-allele-fraction ~{default="0.4" point_size_allele_fraction} \
            --output ~{output_dir_} \
            --output-prefix ~{entity_id}
    >>>

    runtime {
        docker: "~{gatk_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + disk_space_gb + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File modeled_segments_plot = "~{output_dir_}/~{entity_id}.modeled.png"
    }
}
