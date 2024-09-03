version 1.0

import "cnv_germline_case_scattered_workflow.wdl" as CaseCNVScatterWorkflow


workflow CaseAndFilterVCFs {
    input {
        Array[File] vcfs  # Input array of VCF files
        String? filter_expression #= "QUAL < 50"  # Example filter criteria
        String ref_fasta
        String gatk_docker


        ##################################
        ##for case mode
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


    call CaseCNVScatterWorkflow.CNVGermlineCaseScatteredWorkflow as CNVGermlineCaseWorkflow{
        input:
            File intervals,
            File? blacklist_intervals,
            File filtered_intervals,
            Array[String]+ normal_bams,
            Array[String]+ normal_bais,
            File contig_ploidy_model_tar,
            Array[File]+ gcnv_model_tars,
            Int num_intervals_per_scatter,
            File ref_fasta_dict,
            File ref_fasta_fai,
            File ref_fasta,
            String gatk_docker,
            Int num_samples_per_scatter_block,

            ##################################
            #### optional basic arguments ####
            ##################################
            File? gatk4_jar_override,
            Int? preemptible_attempts,

            # Required if BAM/CRAM is in a requester pays bucket
            String? gcs_project_for_requester_pays,

            ####################################################
            #### optional arguments for PreprocessIntervals ####
            ####################################################
            Int? padding,
            Int? bin_length,

            ##############################################
            #### optional arguments for CollectCounts ####
            ##############################################
            Array[String]? disabled_read_filters_for_collect_counts,
            String? collect_counts_format,
            Boolean? collect_counts_enable_indexing,
            Int? mem_gb_for_collect_counts,

            ######################################################################
            #### optional arguments for DetermineGermlineContigPloidyCaseMode ####
            ######################################################################
            Float? ploidy_mapping_error_rate,
            Float? ploidy_sample_psi_scale,
            Int? mem_gb_for_determine_germline_contig_ploidy,
            Int? cpu_for_determine_germline_contig_ploidy,
            Int? disk_for_determine_germline_contig_ploidy,

            ##########################################################
            #### optional arguments for GermlineCNVCallerCaseMode ####
            ##########################################################
            Float? gcnv_p_alt,
            Float? gcnv_cnv_coherence_length,
            Int? gcnv_max_copy_number,
            Int? mem_gb_for_germline_cnv_caller,
            Int? cpu_for_germline_cnv_caller,
            Int? disk_for_germline_cnv_caller,

            # optional arguments for germline CNV denoising model
            Float? gcnv_mapping_error_rate,
            Float? gcnv_sample_psi_scale,
            Float? gcnv_depth_correction_tau,
            String? gcnv_copy_number_posterior_expectation_mode,
            Int? gcnv_active_class_padding_hybrid_mode,

            # optional arguments for Hybrid ADVI
            Float? gcnv_learning_rate,
            Float? gcnv_adamax_beta_1,
            Float? gcnv_adamax_beta_2,
            Int? gcnv_log_emission_samples_per_round,
            Float? gcnv_log_emission_sampling_median_rel_error,
            Int? gcnv_log_emission_sampling_rounds,
            Int? gcnv_max_advi_iter_first_epoch,
            Int? gcnv_max_advi_iter_subsequent_epochs,
            Int? gcnv_min_training_epochs,
            Int? gcnv_max_training_epochs,
            Float? gcnv_initial_temperature,
            Int? gcnv_num_thermal_advi_iters,
            Int? gcnv_convergence_snr_averaging_window,
            Float? gcnv_convergence_snr_trigger_threshold,
            Int? gcnv_convergence_snr_countdown_window,
            Int? gcnv_max_calling_iters,
            Float? gcnv_caller_update_convergence_threshold,
            Float? gcnv_caller_internal_admixing_rate,
            Float? gcnv_caller_external_admixing_rate,
            Boolean? gcnv_disable_annealing,

            ###################################################
            #### arguments for PostprocessGermlineCNVCalls ####
            ###################################################
            Int ref_copy_number_autosomal_contigs,
            Array[String]? allosomal_contigs,

            ##########################
            #### arguments for QC ####
            ##########################
            Int maximum_number_events_per_sample,
            Int maximum_number_pass_events_per_sample
    }



    scatter (vcf in CNVGermlineCaseWorkflow.genotyped_segments_vcfs) {
        call FilterVCF{
            input:
                vcf_file = vcf,
                filter_expression = filter_expression,
                ref_fasta = ref_fasta,
                gatk_docker = gatk_docker
        }
    }


    output {
        Array[File] filtered_vcfs = FilterVCF.filtered_vcf

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

task FilterVCF {
    input {
        File vcf_file
        String? filter_expression
        String ref_fasta
        String gatk_docker
        String? filter_name
    }
    String filter_name_str = select_first([filter_name, "LowQual"])
    String filter_expr = select_first([filter_expression, "(QUAL < 50 and CN>2) or (QUAL < 100 and CN<2)"])

    command <<<
        samplename=$(basename ${vcf_file} | sed -E 's/genotyped-segments-(.*)\.cram\.vcf\.gz/\1/')
        gatk --java-options "-Xmx4g" VariantFiltration \
            -R ${ref_fasta} \
            -V ${vcf_file} \
            --filter-expression "${filter_expression}" \
            --filter-name "LowQual" \
            -O ${samplename}.filtered.genotyped-segments.vcf.gz
    >>>

    output {
        File filtered_vcf = ${samplename}.filtered.genotyped-segments.vcf.gz #"filtered_${basename(vcf_file)}"
    }

    runtime {
        docker: gatk_docker #"biocontainers/bcftools:v1.10.2-1-deb_cv1"
        memory: "8G"
        cpu: 2
    }
}