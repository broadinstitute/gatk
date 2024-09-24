version 1.0

import "cnv_germline_case_scattered_workflow.wdl" as CaseCNVScatterWorkflow


workflow CaseAndFilterVCFs {
    input {
        String? filter_expression #= "QUAL < 50"  # Example filter criteria


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


    call CaseCNVScatterWorkflow.CNVGermlineCaseScatteredWorkflow as CNVGermlineCaseWorkflow {
        input:
            intervals = intervals,
            blacklist_intervals=blacklist_intervals,
            filtered_intervals=filtered_intervals,
            normal_bams=normal_bams,
            normal_bais=normal_bais,
            contig_ploidy_model_tar=contig_ploidy_model_tar,
            gcnv_model_tars=gcnv_model_tars,
            num_intervals_per_scatter=num_intervals_per_scatter,
            ref_fasta_dict=ref_fasta_dict,
            ref_fasta_fai=ref_fasta_fai,
            ref_fasta=ref_fasta,
            gatk_docker=gatk_docker,
            num_samples_per_scatter_block=num_samples_per_scatter_block,

            ##################################
            #### optional basic arguments ####
            ##################################
            gatk4_jar_override=gatk4_jar_override,
            preemptible_attempts=preemptible_attempts,

            # Required if BAM/CRAM is in a requester pays bucket
            gcs_project_for_requester_pays=gcs_project_for_requester_pays,

            ####################################################
            #### optional arguments for PreprocessIntervals ####
            ####################################################
            padding=padding,
            bin_length=bin_length,

            ##############################################
            #### optional arguments for CollectCounts ####
            ##############################################
            disabled_read_filters_for_collect_counts=disabled_read_filters_for_collect_counts,
            collect_counts_format=collect_counts_format,
            collect_counts_enable_indexing=collect_counts_enable_indexing,
            mem_gb_for_collect_counts=mem_gb_for_collect_counts,

            ######################################################################
            #### optional arguments for DetermineGermlineContigPloidyCaseMode ####
            ######################################################################
            ploidy_mapping_error_rate=ploidy_mapping_error_rate,
            ploidy_sample_psi_scale=ploidy_sample_psi_scale,
            mem_gb_for_determine_germline_contig_ploidy=mem_gb_for_determine_germline_contig_ploidy,
            cpu_for_determine_germline_contig_ploidy=cpu_for_determine_germline_contig_ploidy,
            disk_for_determine_germline_contig_ploidy=disk_for_determine_germline_contig_ploidy,

            ##########################################################
            #### optional arguments for GermlineCNVCallerCaseMode ####
            ##########################################################
            gcnv_p_alt=gcnv_p_alt,
            gcnv_cnv_coherence_length=gcnv_cnv_coherence_length,
            gcnv_max_copy_number=gcnv_max_copy_number,
            mem_gb_for_germline_cnv_caller=mem_gb_for_germline_cnv_caller,
            cpu_for_germline_cnv_caller=cpu_for_germline_cnv_caller,
            disk_for_germline_cnv_caller=disk_for_germline_cnv_caller,

            # optional arguments for germline CNV denoising model
            gcnv_mapping_error_rate=gcnv_mapping_error_rate,
            gcnv_sample_psi_scale=gcnv_sample_psi_scale,
            gcnv_depth_correction_tau=gcnv_depth_correction_tau,
            gcnv_copy_number_posterior_expectation_mode=gcnv_copy_number_posterior_expectation_mode,
            gcnv_active_class_padding_hybrid_mode=gcnv_active_class_padding_hybrid_mode,

            # optional arguments for Hybrid ADVI
            gcnv_learning_rate=gcnv_learning_rate,
            gcnv_adamax_beta_1=gcnv_adamax_beta_1,
            gcnv_adamax_beta_2=gcnv_adamax_beta_2,
            gcnv_log_emission_samples_per_round=gcnv_log_emission_samples_per_round,
            gcnv_log_emission_sampling_median_rel_error=gcnv_log_emission_sampling_median_rel_error,
            gcnv_log_emission_sampling_rounds=gcnv_log_emission_sampling_rounds,
            gcnv_max_advi_iter_first_epoch=gcnv_max_advi_iter_first_epoch,
            gcnv_max_advi_iter_subsequent_epochs=gcnv_max_advi_iter_subsequent_epochs,
            gcnv_min_training_epochs=gcnv_min_training_epochs,
            gcnv_max_training_epochs=gcnv_max_training_epochs,
            gcnv_initial_temperature=gcnv_initial_temperature,
            gcnv_num_thermal_advi_iters=gcnv_num_thermal_advi_iters,
            gcnv_convergence_snr_averaging_window=gcnv_convergence_snr_averaging_window,
            gcnv_convergence_snr_trigger_threshold=gcnv_convergence_snr_trigger_threshold,
            gcnv_convergence_snr_countdown_window=gcnv_convergence_snr_countdown_window,
            gcnv_max_calling_iters=gcnv_max_calling_iters,
            gcnv_caller_update_convergence_threshold=gcnv_caller_update_convergence_threshold,
            gcnv_caller_internal_admixing_rate=gcnv_caller_internal_admixing_rate,
            gcnv_caller_external_admixing_rate=gcnv_caller_external_admixing_rate,
            gcnv_disable_annealing=gcnv_disable_annealing,

            ###################################################
            #### arguments for PostprocessGermlineCNVCalls ####
            ###################################################
            ref_copy_number_autosomal_contigs=ref_copy_number_autosomal_contigs,
            allosomal_contigs=allosomal_contigs,

            ##########################
            #### arguments for QC ####
            ##########################
            maximum_number_events_per_sample=maximum_number_events_per_sample,
            maximum_number_pass_events_per_sample=maximum_number_pass_events_per_sample
    }

scatter (vcf in CNVGermlineCaseWorkflow.genotyped_segments_vcfs) {
        call ExtractSamplename {
            input:
                vcf_filename = basename(vcf),
                gatk_docker = gatk_docker
        }
        call FilterVCF{
            input:
                samplename = ExtractSamplename.samplename,
                vcf_file = vcf,
                filter_expression = filter_expression,
                ref_fasta = ref_fasta,
                gatk_docker = gatk_docker
        }
    }


    output {
        Array[File] filtered_vcfs = FilterVCF.filtered_vcf

        File preprocessed_intervals = CNVGermlineCaseWorkflow.preprocessed_intervals
        Array[File] read_counts_entity_id = CNVGermlineCaseWorkflow.read_counts_entity_id
        Array[File] read_counts = CNVGermlineCaseWorkflow.read_counts
        Array[File] sample_contig_ploidy_calls_tars = CNVGermlineCaseWorkflow.sample_contig_ploidy_calls_tars
        Array[Array[File]] gcnv_calls_tars = CNVGermlineCaseWorkflow.gcnv_calls_tars
        Array[File] gcnv_tracking_tars = CNVGermlineCaseWorkflow.gcnv_tracking_tars
        Array[File] genotyped_intervals_vcfs = CNVGermlineCaseWorkflow.genotyped_intervals_vcfs
        Array[File] genotyped_intervals_vcf_indexes = CNVGermlineCaseWorkflow.genotyped_intervals_vcf_indexes
        Array[File] genotyped_segments_vcfs = CNVGermlineCaseWorkflow.genotyped_segments_vcfs
        Array[File] genotyped_segments_vcf_indexes = CNVGermlineCaseWorkflow.genotyped_segments_vcf_indexes
        Array[File] qc_status_files = CNVGermlineCaseWorkflow.qc_status_files
        Array[String] qc_status_strings = CNVGermlineCaseWorkflow.qc_status_strings
        Array[File] denoised_copy_ratios = CNVGermlineCaseWorkflow.denoised_copy_ratios
    }
}


task ExtractSamplename {
    input {
        String vcf_filename
        String gatk_docker
    }
    command <<<
        echo ~{vcf_filename} | sed -E 's/genotyped-segments-(.*)\.cram\.vcf\.gz/\1/' > output.txt
        >>>
    output {
        String samplename = read_string("output.txt")
    }
    runtime {
        docker: gatk_docker
        memory: "2G"
        cpu: 2
    }
}


task FilterVCF {
    input {
        String samplename
        File vcf_file
        String filter_expression = "(FMT/CN>2 & QUAL<50) | (FMT/CN==1 & QUAL<100 ) | (FMT/CN==0 & QUAL<400)"
        String gatk_docker
        String filter_name = "LowQual"
    }

    command <<<
    	bcftools filter \
		-e '~{filter_expression}' \
		--soft-filter ~{filter_name} \
		-o ~{samplename}.filtered.genotyped-segments.vcf.gz \
		~{vcf_file}
    >>>

    output {
        File filtered_vcf = samplename + ".filtered.genotyped-segments.vcf.gz" #"filtered_${basename(vcf_file)}"
    }

    runtime {
        docker: gatk_docker #"biocontainers/bcftools:v1.10.2-1-deb_cv1"
        memory: "8G"
        cpu: 2
    }
}