#
# A wrapper for the gCNV case workflow intended for lowering computing cost by making it feasible to use 
# preemptible cloud instances with low memory requirements. CPU, memory and disk requirements can be 
# lowered for GermlineCNVCaller and DetermineGermlineContigPloidy tasks.
# 
#
# - Example invocation:
#
#       java -jar cromwell.jar run cnv_germline_case_scattered_workflow.wdl -i my_parameters.json
#
####################

version 1.0

import "cnv_germline_case_workflow.wdl" as GermlineCNVCaseWorkflow
import "../cnv_common_tasks.wdl" as CNVTasks

workflow CNVGermlineCaseScatteredWorkflow {

    input {
      ##################################
      #### required basic arguments ####
      ##################################
      File intervals
      File? blacklist_intervals
      File filtered_intervals
      Array[String]+ normal_bams
      Array[String]+ normal_bais
      File contig_ploidy_model_tar
      Array[File]+ gcnv_model_tars
      Int num_intervals_per_scatter
      File ref_fasta_dict
      File ref_fasta_fai
      File ref_fasta
      String gatk_docker
      Int num_samples_per_scatter_block

      ##################################
      #### optional basic arguments ####
      ##################################
      File? gatk4_jar_override
      Int? preemptible_attempts

      # Required if BAM/CRAM is in a requester pays bucket
      String? gcs_project_for_requester_pays

      ####################################################
      #### optional arguments for PreprocessIntervals ####
      ####################################################
      Int? padding
      Int? bin_length

      ##############################################
      #### optional arguments for CollectCounts ####
      ##############################################
      Array[String]? disabled_read_filters_for_collect_counts
      String? collect_counts_format
      Boolean? collect_counts_enable_indexing
      Int? mem_gb_for_collect_counts

      ######################################################################
      #### optional arguments for DetermineGermlineContigPloidyCaseMode ####
      ######################################################################
      Float? ploidy_mapping_error_rate
      Float? ploidy_sample_psi_scale
      Int? mem_gb_for_determine_germline_contig_ploidy
      Int? cpu_for_determine_germline_contig_ploidy
      Int? disk_for_determine_germline_contig_ploidy

      ##########################################################
      #### optional arguments for GermlineCNVCallerCaseMode ####
      ##########################################################
      Float? gcnv_p_alt
      Float? gcnv_cnv_coherence_length
      Int? gcnv_max_copy_number
      Int? mem_gb_for_germline_cnv_caller
      Int? cpu_for_germline_cnv_caller
      Int? disk_for_germline_cnv_caller

      # optional arguments for germline CNV denoising model
      Float? gcnv_mapping_error_rate
      Float? gcnv_sample_psi_scale
      Float? gcnv_depth_correction_tau
      String? gcnv_copy_number_posterior_expectation_mode
      Int? gcnv_active_class_padding_hybrid_mode

      # optional arguments for Hybrid ADVI
      Float? gcnv_learning_rate
      Float? gcnv_adamax_beta_1
      Float? gcnv_adamax_beta_2
      Int? gcnv_log_emission_samples_per_round
      Float? gcnv_log_emission_sampling_median_rel_error
      Int? gcnv_log_emission_sampling_rounds
      Int? gcnv_max_advi_iter_first_epoch
      Int? gcnv_max_advi_iter_subsequent_epochs
      Int? gcnv_min_training_epochs
      Int? gcnv_max_training_epochs
      Float? gcnv_initial_temperature
      Int? gcnv_num_thermal_advi_iters
      Int? gcnv_convergence_snr_averaging_window
      Float? gcnv_convergence_snr_trigger_threshold
      Int? gcnv_convergence_snr_countdown_window
      Int? gcnv_max_calling_iters
      Float? gcnv_caller_update_convergence_threshold
      Float? gcnv_caller_internal_admixing_rate
      Float? gcnv_caller_external_admixing_rate
      Boolean? gcnv_disable_annealing

      ###################################################
      #### arguments for PostprocessGermlineCNVCalls ####
      ###################################################
      Int ref_copy_number_autosomal_contigs
      Array[String]? allosomal_contigs

      ##########################
      #### arguments for QC ####
      ##########################
      Int maximum_number_events_per_sample
      Int maximum_number_pass_events_per_sample
    }

    call CNVTasks.SplitInputArray as SplitInputBamsList {
        input:
            input_array = normal_bams,
            num_inputs_in_scatter_block = num_samples_per_scatter_block,
            gatk_docker = gatk_docker
    }

    call CNVTasks.SplitInputArray as SplitInputBaisList {
        input:
            input_array = normal_bais,
            num_inputs_in_scatter_block = num_samples_per_scatter_block,
            gatk_docker = gatk_docker
    }

    Array[Array[String]] split_bams = SplitInputBamsList.split_array
    Array[Array[String]] split_bais = SplitInputBaisList.split_array

    scatter (subarray_index in range(length(split_bams))) {
        call GermlineCNVCaseWorkflow.CNVGermlineCaseWorkflow {
            input:
                intervals = intervals,
                blacklist_intervals = blacklist_intervals,
                filtered_intervals = filtered_intervals,
                normal_bams = split_bams[subarray_index],
                normal_bais = split_bais[subarray_index],
                contig_ploidy_model_tar = contig_ploidy_model_tar,
                gcnv_model_tars = gcnv_model_tars,
                num_intervals_per_scatter = num_intervals_per_scatter,
                ref_fasta_dict = ref_fasta_dict,
                ref_fasta_fai = ref_fasta_fai,
                ref_fasta = ref_fasta,
                gatk_docker = gatk_docker,
                gatk4_jar_override = gatk4_jar_override,
                preemptible_attempts = preemptible_attempts,
                gcs_project_for_requester_pays = gcs_project_for_requester_pays,
                padding = padding,
                bin_length = bin_length,
                disabled_read_filters_for_collect_counts = disabled_read_filters_for_collect_counts,
                collect_counts_format = collect_counts_format,
                collect_counts_enable_indexing = collect_counts_enable_indexing,
                mem_gb_for_collect_counts = mem_gb_for_collect_counts,
                ploidy_mapping_error_rate = ploidy_mapping_error_rate,
                ploidy_sample_psi_scale = ploidy_sample_psi_scale,
                mem_gb_for_determine_germline_contig_ploidy = mem_gb_for_determine_germline_contig_ploidy,
                cpu_for_determine_germline_contig_ploidy = cpu_for_determine_germline_contig_ploidy,
                disk_for_determine_germline_contig_ploidy = disk_for_determine_germline_contig_ploidy,
                gcnv_p_alt = gcnv_p_alt,
                gcnv_cnv_coherence_length = gcnv_cnv_coherence_length,
                gcnv_max_copy_number = gcnv_max_copy_number,
                mem_gb_for_germline_cnv_caller = mem_gb_for_germline_cnv_caller,
                cpu_for_germline_cnv_caller = cpu_for_germline_cnv_caller,
                disk_for_germline_cnv_caller = disk_for_germline_cnv_caller,
                gcnv_mapping_error_rate = gcnv_mapping_error_rate,
                gcnv_sample_psi_scale = gcnv_sample_psi_scale,
                gcnv_depth_correction_tau = gcnv_depth_correction_tau,
                gcnv_copy_number_posterior_expectation_mode = gcnv_copy_number_posterior_expectation_mode,
                gcnv_active_class_padding_hybrid_mode = gcnv_active_class_padding_hybrid_mode,
                gcnv_learning_rate = gcnv_learning_rate,
                gcnv_adamax_beta_1 = gcnv_adamax_beta_1,
                gcnv_adamax_beta_2 = gcnv_adamax_beta_2,
                gcnv_log_emission_samples_per_round = gcnv_log_emission_samples_per_round,
                gcnv_log_emission_sampling_median_rel_error = gcnv_log_emission_sampling_median_rel_error,
                gcnv_log_emission_sampling_rounds = gcnv_log_emission_sampling_rounds,
                gcnv_max_advi_iter_first_epoch = gcnv_max_advi_iter_first_epoch,
                gcnv_max_advi_iter_subsequent_epochs = gcnv_max_advi_iter_subsequent_epochs,
                gcnv_min_training_epochs = gcnv_min_training_epochs,
                gcnv_max_training_epochs = gcnv_max_training_epochs,
                gcnv_initial_temperature = gcnv_initial_temperature,
                gcnv_num_thermal_advi_iters = gcnv_num_thermal_advi_iters,
                gcnv_convergence_snr_averaging_window = gcnv_convergence_snr_averaging_window,
                gcnv_convergence_snr_trigger_threshold = gcnv_convergence_snr_trigger_threshold,
                gcnv_convergence_snr_countdown_window = gcnv_convergence_snr_countdown_window,
                gcnv_max_calling_iters = gcnv_max_calling_iters,
                gcnv_caller_update_convergence_threshold = gcnv_caller_update_convergence_threshold,
                gcnv_caller_internal_admixing_rate = gcnv_caller_internal_admixing_rate,
                gcnv_caller_external_admixing_rate = gcnv_caller_external_admixing_rate,
                gcnv_disable_annealing = gcnv_disable_annealing,
                ref_copy_number_autosomal_contigs = ref_copy_number_autosomal_contigs,
                allosomal_contigs = allosomal_contigs,
                maximum_number_events_per_sample = maximum_number_events_per_sample,
                maximum_number_pass_events_per_sample = maximum_number_pass_events_per_sample
        }
    }

    output {
        File preprocessed_intervals = CNVGermlineCaseWorkflow.preprocessed_intervals[0]
        Array[File] read_counts_entity_id = flatten(CNVGermlineCaseWorkflow.read_counts_entity_id)
        Array[File] read_counts = flatten(CNVGermlineCaseWorkflow.read_counts)
        Array[File] sample_contig_ploidy_calls_tars = flatten(CNVGermlineCaseWorkflow.sample_contig_ploidy_calls_tars)
        Array[Array[File]] gcnv_calls_tars = flatten(CNVGermlineCaseWorkflow.gcnv_calls_tars)
        Array[File] gcnv_tracking_tars = flatten(CNVGermlineCaseWorkflow.gcnv_tracking_tars)
        Array[File] genotyped_intervals_vcfs = flatten(CNVGermlineCaseWorkflow.genotyped_intervals_vcfs)
        Array[File] genotyped_intervals_vcf_indexes = flatten(CNVGermlineCaseWorkflow.genotyped_intervals_vcf_indexes)
        Array[File] genotyped_segments_vcfs = flatten(CNVGermlineCaseWorkflow.genotyped_segments_vcfs)
        Array[File] genotyped_segments_vcf_indexes = flatten(CNVGermlineCaseWorkflow.genotyped_segments_vcf_indexes)
        Array[File] qc_status_files = flatten(CNVGermlineCaseWorkflow.qc_status_files)
        Array[String] qc_status_strings = flatten(CNVGermlineCaseWorkflow.qc_status_strings)
        Array[File] denoised_copy_ratios = flatten(CNVGermlineCaseWorkflow.denoised_copy_ratios)
    }
}

