# Workflow for running GATK GermlineCNVCaller on multiple case samples using a trained model (obtained from running 
# GATK GermlineCNVCaller in the cohort mode). Supports both WGS and WES.
#
# Notes:
#
# - The intervals argument is required for both WGS and WES workflows and accepts formats compatible with the
#   GATK -L argument (see https://gatkforums.broadinstitute.org/gatk/discussion/11009/intervals-and-interval-lists).
#   These intervals will be padded on both sides by the amount specified by padding (default 250)
#   and split into bins of length specified by bin_length (default 1000; specify 0 to skip binning,
#   e.g., for WES).  For WGS, the intervals should simply cover the chromosomes of interest.
#
# - Intervals can be blacklisted from coverage collection and all downstream steps by using the blacklist_intervals
#   argument, which accepts formats compatible with the GATK -XL argument
#   (see https://gatkforums.broadinstitute.org/gatk/discussion/11009/intervals-and-interval-lists).
#   This may be useful for excluding centromeric regions, etc. from analysis.  Alternatively, these regions may
#   be manually filtered from the final callset.
#
# - Example invocation:
#
#       java -jar cromwell.jar run cnv_germline_case_workflow.wdl -i my_parameters.json
#
#############

version 1.0

import "../cnv_common_tasks.wdl" as CNVTasks

workflow CNVGermlineCaseWorkflow {

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
      Int? disk_space_gb_for_postprocess_germline_cnv_calls
      Int? mem_gb_for_postprocess_germline_cnv_calls

      ##########################
      #### arguments for QC ####
      ##########################
      Int maximum_number_events_per_sample
      Int maximum_number_pass_events_per_sample
    }

    Array[Pair[String, String]] normal_bams_and_bais = zip(normal_bams, normal_bais)

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
            preemptible_attempts = preemptible_attempts
    }

    scatter (normal_bam_and_bai in normal_bams_and_bais) {
        call CNVTasks.CollectCounts {
            input:
                intervals = PreprocessIntervals.preprocessed_intervals,
                bam = normal_bam_and_bai.left,
                bam_idx = normal_bam_and_bai.right,
                ref_fasta = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,
                ref_fasta_dict = ref_fasta_dict,
                format = collect_counts_format,
                enable_indexing = collect_counts_enable_indexing,
                disabled_read_filters = disabled_read_filters_for_collect_counts,
                gatk4_jar_override = gatk4_jar_override,
                gatk_docker = gatk_docker,
                mem_gb = mem_gb_for_collect_counts,
                preemptible_attempts = preemptible_attempts,
                gcs_project_for_requester_pays = gcs_project_for_requester_pays
        }
    }

    call DetermineGermlineContigPloidyCaseMode {
        input:
            read_count_files = CollectCounts.counts,
            contig_ploidy_model_tar = contig_ploidy_model_tar,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            mem_gb = mem_gb_for_determine_germline_contig_ploidy,
            cpu = cpu_for_determine_germline_contig_ploidy,
            disk_space_gb = disk_for_determine_germline_contig_ploidy,
            mapping_error_rate = ploidy_mapping_error_rate,
            sample_psi_scale = ploidy_sample_psi_scale,
            preemptible_attempts = preemptible_attempts
    }

    call CNVTasks.ScatterIntervals {
        input:
            interval_list = filtered_intervals,
            num_intervals_per_scatter = num_intervals_per_scatter,
            gatk_docker = gatk_docker,
            preemptible_attempts = preemptible_attempts
    }

    scatter (scatter_index in range(length(ScatterIntervals.scattered_interval_lists))) {
        call GermlineCNVCallerCaseMode {
            input:
                scatter_index = scatter_index,
                read_count_files = CollectCounts.counts,
                contig_ploidy_calls_tar = DetermineGermlineContigPloidyCaseMode.contig_ploidy_calls_tar,
                gcnv_model_tar = gcnv_model_tars[scatter_index],
                gatk4_jar_override = gatk4_jar_override,
                gatk_docker = gatk_docker,
                mem_gb = mem_gb_for_germline_cnv_caller,
                cpu = cpu_for_germline_cnv_caller,
                p_alt = gcnv_p_alt,
                cnv_coherence_length = gcnv_cnv_coherence_length,
                max_copy_number = gcnv_max_copy_number,
                mapping_error_rate = gcnv_mapping_error_rate,
                sample_psi_scale = gcnv_sample_psi_scale,
                depth_correction_tau = gcnv_depth_correction_tau,
                copy_number_posterior_expectation_mode = gcnv_copy_number_posterior_expectation_mode,
                active_class_padding_hybrid_mode = gcnv_active_class_padding_hybrid_mode,
                learning_rate = gcnv_learning_rate,
                adamax_beta_1 = gcnv_adamax_beta_1,
                adamax_beta_2 = gcnv_adamax_beta_2,
                log_emission_samples_per_round = gcnv_log_emission_samples_per_round,
                log_emission_sampling_median_rel_error = gcnv_log_emission_sampling_median_rel_error,
                log_emission_sampling_rounds = gcnv_log_emission_sampling_rounds,
                max_advi_iter_first_epoch = gcnv_max_advi_iter_first_epoch,
                max_advi_iter_subsequent_epochs = gcnv_max_advi_iter_subsequent_epochs,
                min_training_epochs = gcnv_min_training_epochs,
                max_training_epochs = gcnv_max_training_epochs,
                initial_temperature = gcnv_initial_temperature,
                num_thermal_advi_iters = gcnv_num_thermal_advi_iters,
                convergence_snr_averaging_window = gcnv_convergence_snr_averaging_window,
                convergence_snr_trigger_threshold = gcnv_convergence_snr_trigger_threshold,
                convergence_snr_countdown_window = gcnv_convergence_snr_countdown_window,
                max_calling_iters = gcnv_max_calling_iters,
                caller_update_convergence_threshold = gcnv_caller_update_convergence_threshold,
                caller_internal_admixing_rate = gcnv_caller_internal_admixing_rate,
                caller_external_admixing_rate = gcnv_caller_external_admixing_rate,
                disable_annealing = gcnv_disable_annealing,
                preemptible_attempts = preemptible_attempts
        }
    }

    Array[Array[File]] call_tars_sample_by_shard = transpose(GermlineCNVCallerCaseMode.gcnv_call_tars)

    scatter (sample_index in range(length(normal_bams))) {
        call CNVTasks.PostprocessGermlineCNVCalls {
            input:
                entity_id = CollectCounts.entity_id[sample_index],
                gcnv_calls_tars = call_tars_sample_by_shard[sample_index],
                gcnv_model_tars = gcnv_model_tars,
                calling_configs = GermlineCNVCallerCaseMode.calling_config_json,
                denoising_configs = GermlineCNVCallerCaseMode.denoising_config_json,
                gcnvkernel_version = GermlineCNVCallerCaseMode.gcnvkernel_version_json,
                sharded_interval_lists = GermlineCNVCallerCaseMode.sharded_interval_list,
                allosomal_contigs = allosomal_contigs,
                ref_copy_number_autosomal_contigs = ref_copy_number_autosomal_contigs,
                contig_ploidy_calls_tar = DetermineGermlineContigPloidyCaseMode.contig_ploidy_calls_tar,
                sample_index = sample_index,
                maximum_number_events = maximum_number_events_per_sample,
                maximum_number_pass_events = maximum_number_pass_events_per_sample,
                gatk4_jar_override = gatk4_jar_override,
                gatk_docker = gatk_docker,
                preemptible_attempts = preemptible_attempts
        }
    }

    call CNVTasks.ScatterPloidyCallsBySample {
        input :
            contig_ploidy_calls_tar = DetermineGermlineContigPloidyCaseMode.contig_ploidy_calls_tar,
            samples = CollectCounts.entity_id,
            docker = gatk_docker,
            preemptible_attempts = preemptible_attempts
    }

    output {
        File preprocessed_intervals = PreprocessIntervals.preprocessed_intervals
        Array[File] read_counts_entity_id = CollectCounts.entity_id
        Array[File] read_counts = CollectCounts.counts
        Array[File] sample_contig_ploidy_calls_tars = ScatterPloidyCallsBySample.sample_contig_ploidy_calls_tar
        Array[Array[File]] gcnv_calls_tars = GermlineCNVCallerCaseMode.gcnv_call_tars
        Array[File] gcnv_tracking_tars = GermlineCNVCallerCaseMode.gcnv_tracking_tar
        Array[File] genotyped_intervals_vcfs = PostprocessGermlineCNVCalls.genotyped_intervals_vcf
        Array[File] genotyped_intervals_vcf_indexes = PostprocessGermlineCNVCalls.genotyped_intervals_vcf_index
        Array[File] genotyped_segments_vcfs = PostprocessGermlineCNVCalls.genotyped_segments_vcf
        Array[File] genotyped_segments_vcf_indexes = PostprocessGermlineCNVCalls.genotyped_segments_vcf_index
        Array[File] qc_status_files = PostprocessGermlineCNVCalls.qc_status_file
        Array[String] qc_status_strings = PostprocessGermlineCNVCalls.qc_status_string
        Array[File] denoised_copy_ratios = PostprocessGermlineCNVCalls.denoised_copy_ratios
    }
}

task DetermineGermlineContigPloidyCaseMode {
    input {
      Array[File] read_count_files
      File contig_ploidy_model_tar
      String? output_dir
      File? gatk4_jar_override

      # Runtime parameters
      String gatk_docker
      Int? mem_gb
      Int? disk_space_gb
      Boolean use_ssd = false
      Int? cpu
      Int? preemptible_attempts

      # Model parameters
      Float? mapping_error_rate
      Float? sample_psi_scale
    }

    # We do not expose Hybrid ADVI parameters -- the default values are decent

    Int machine_mem_mb = select_first([mem_gb, 7]) * 1000
    Int command_mem_mb = machine_mem_mb - 500

    # If optional output_dir not specified, use "out"
    String output_dir_ = select_first([output_dir, "out"])

    command <<<
        set -eu
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}
        export MKL_NUM_THREADS=~{default=8 cpu}
        export OMP_NUM_THREADS=~{default=8 cpu}

        mkdir contig-ploidy-model
        tar xzf ~{contig_ploidy_model_tar} -C contig-ploidy-model

        gatk --java-options "-Xmx~{command_mem_mb}m" DetermineGermlineContigPloidy \
            --input ~{sep=" --input " read_count_files} \
            --model contig-ploidy-model \
            --output ~{output_dir_} \
            --output-prefix case \
            --verbosity DEBUG \
            --mapping-error-rate ~{default="0.01" mapping_error_rate} \
            --sample-psi-scale ~{default="0.0001" sample_psi_scale}

        tar c -C ~{output_dir_}/case-calls . | gzip -1 > case-contig-ploidy-calls.tar.gz

        rm -rf contig-ploidy-model
    >>>

    runtime {
        docker: gatk_docker
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 150]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 8])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File contig_ploidy_calls_tar = "case-contig-ploidy-calls.tar.gz"
    }
}

task GermlineCNVCallerCaseMode {
    input {
      Int scatter_index
      Array[File] read_count_files
      File contig_ploidy_calls_tar
      File gcnv_model_tar
      String? output_dir
      File? gatk4_jar_override

      # Runtime parameters
      String gatk_docker
      Int? mem_gb
      Int? disk_space_gb
      Boolean use_ssd = false
      Int? cpu
      Int? preemptible_attempts

      # Caller parameters
      Float? p_alt
      Float? cnv_coherence_length
      Int? max_copy_number

      # Denoising model parameters
      Float? mapping_error_rate
      Float? sample_psi_scale
      Float? depth_correction_tau
      String? copy_number_posterior_expectation_mode
      Int? active_class_padding_hybrid_mode

      # Hybrid ADVI parameters
      Float? learning_rate
      Float? adamax_beta_1
      Float? adamax_beta_2
      Int? log_emission_samples_per_round
      Float? log_emission_sampling_median_rel_error
      Int? log_emission_sampling_rounds
      Int? max_advi_iter_first_epoch
      Int? max_advi_iter_subsequent_epochs
      Int? min_training_epochs
      Int? max_training_epochs
      Float? initial_temperature
      Int? num_thermal_advi_iters
      Int? convergence_snr_averaging_window
      Float? convergence_snr_trigger_threshold
      Int? convergence_snr_countdown_window
      Int? max_calling_iters
      Float? caller_update_convergence_threshold
      Float? caller_internal_admixing_rate
      Float? caller_external_admixing_rate
      Boolean? disable_annealing
    }

    Int machine_mem_mb = select_first([mem_gb, 7]) * 1000
    Int command_mem_mb = machine_mem_mb - 500

    # If optional output_dir not specified, use "out"
    String output_dir_ = select_first([output_dir, "out"])
    Int num_samples = length(read_count_files)

    command <<<
        set -eu
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}
        export MKL_NUM_THREADS=~{default=8 cpu}
        export OMP_NUM_THREADS=~{default=8 cpu}

        mkdir contig-ploidy-calls
        tar xzf ~{contig_ploidy_calls_tar} -C contig-ploidy-calls

        mkdir gcnv-model
        tar xzf ~{gcnv_model_tar} -C gcnv-model

        gatk --java-options "-Xmx~{command_mem_mb}m"  GermlineCNVCaller \
            --run-mode CASE \
            --input ~{sep=" --input " read_count_files} \
            --contig-ploidy-calls contig-ploidy-calls \
            --model gcnv-model \
            --output ~{output_dir_} \
            --output-prefix case \
            --verbosity DEBUG \
            --p-alt ~{default="1e-6" p_alt} \
            --cnv-coherence-length ~{default="10000.0" cnv_coherence_length} \
            --max-copy-number ~{default="5" max_copy_number} \
            --mapping-error-rate ~{default="0.01" mapping_error_rate} \
            --sample-psi-scale ~{default="0.0001" sample_psi_scale} \
            --depth-correction-tau ~{default="10000.0" depth_correction_tau} \
            --copy-number-posterior-expectation-mode ~{default="HYBRID" copy_number_posterior_expectation_mode} \
            --active-class-padding-hybrid-mode ~{default="50000" active_class_padding_hybrid_mode} \
            --learning-rate ~{default="0.05" learning_rate} \
            --adamax-beta-1 ~{default="0.9" adamax_beta_1} \
            --adamax-beta-2 ~{default="0.99" adamax_beta_2} \
            --log-emission-samples-per-round ~{default="50" log_emission_samples_per_round} \
            --log-emission-sampling-median-rel-error ~{default="0.005" log_emission_sampling_median_rel_error} \
            --log-emission-sampling-rounds ~{default="10" log_emission_sampling_rounds} \
            --max-advi-iter-first-epoch ~{default="5000" max_advi_iter_first_epoch} \
            --max-advi-iter-subsequent-epochs ~{default="100" max_advi_iter_subsequent_epochs} \
            --min-training-epochs ~{default="10" min_training_epochs} \
            --max-training-epochs ~{default="100" max_training_epochs} \
            --initial-temperature ~{default="2.0" initial_temperature} \
            --num-thermal-advi-iters ~{default="2500" num_thermal_advi_iters} \
            --convergence-snr-averaging-window ~{default="500" convergence_snr_averaging_window} \
            --convergence-snr-trigger-threshold ~{default="0.1" convergence_snr_trigger_threshold} \
            --convergence-snr-countdown-window ~{default="10" convergence_snr_countdown_window} \
            --max-calling-iters ~{default="10" max_calling_iters} \
            --caller-update-convergence-threshold ~{default="0.001" caller_update_convergence_threshold} \
            --caller-internal-admixing-rate ~{default="0.75" caller_internal_admixing_rate} \
            --caller-external-admixing-rate ~{default="1.00" caller_external_admixing_rate} \
            --disable-annealing ~{default="false" disable_annealing}

        tar czf case-gcnv-tracking-shard-~{scatter_index}.tar.gz -C ~{output_dir_}/case-tracking .

        CURRENT_SAMPLE=0
        NUM_SAMPLES=~{num_samples}
        NUM_DIGITS=${#NUM_SAMPLES}
        while [ $CURRENT_SAMPLE -lt $NUM_SAMPLES ]; do
            CURRENT_SAMPLE_WITH_LEADING_ZEROS=$(printf "%0${NUM_DIGITS}d" $CURRENT_SAMPLE)
            tar czf case-gcnv-calls-shard-~{scatter_index}-sample-$CURRENT_SAMPLE_WITH_LEADING_ZEROS.tar.gz -C ~{output_dir_}/case-calls/SAMPLE_$CURRENT_SAMPLE .
            let CURRENT_SAMPLE=CURRENT_SAMPLE+1
        done

        rm -rf contig-ploidy-calls
        rm -rf gcnv-model
    >>>

    runtime {
        docker: gatk_docker
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 150]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 8])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        Array[File] gcnv_call_tars = glob("case-gcnv-calls-shard-~{scatter_index}-sample-*.tar.gz")
        File gcnv_tracking_tar = "case-gcnv-tracking-shard-~{scatter_index}.tar.gz"
        File calling_config_json = "~{output_dir_}/case-calls/calling_config.json"
        File denoising_config_json = "~{output_dir_}/case-calls/denoising_config.json"
        File gcnvkernel_version_json = "~{output_dir_}/case-calls/gcnvkernel_version.json"
        File sharded_interval_list = "~{output_dir_}/case-calls/interval_list.tsv"
    }
}
