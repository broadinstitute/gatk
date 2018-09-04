# Workflow for running GATK GermlineCNVCaller on a single case sample. Supports both WGS and WES.
#
# Notes:
#
# - The intervals argument is required for both WGS and WES workflows and accepts formats compatible with the
#   GATK -L argument (see https://gatkforums.broadinstitute.org/gatk/discussion/11009/intervals-and-interval-lists).
#   These intervals will be padded on both sides by the amount specified by PreprocessIntervals.padding (default 250)
#   and split into bins of length specified by PreprocessIntervals.bin_length (default 1000; specify 0 to skip binning,
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

import "cnv_common_tasks.wdl" as CNVTasks

workflow CNVGermlineCaseWorkflow {

    ##################################
    #### required basic arguments ####
    ##################################
    File intervals
    File? blacklist_intervals
    File bam
    File bam_idx
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

    ####################################################
    #### optional arguments for PreprocessIntervals ####
    ####################################################
    Int? padding
    Int? bin_length

    ######################################################################
    #### optional arguments for DetermineGermlineContigPloidyCaseMode ####
    ######################################################################
    Float? ploidy_mapping_error_rate
    Float? ploidy_sample_psi_scale
    Int? mem_gb_for_determine_germline_contig_ploidy
    Int? cpu_for_determine_germline_contig_ploidy

    ##########################################################
    #### optional arguments for GermlineCNVCallerCaseMode ####
    ##########################################################
    Float? gcnv_p_alt
    Float? gcnv_cnv_coherence_length
    Int? gcnv_max_copy_number
    Int? mem_gb_for_germline_cnv_caller
    Int? cpu_for_germline_cnv_caller

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

    call CNVTasks.CollectCounts {
        input:
            intervals = PreprocessIntervals.preprocessed_intervals,
            bam = bam,
            bam_idx = bam_idx,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            ref_fasta_dict = ref_fasta_dict,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            preemptible_attempts = preemptible_attempts
    }

    call DetermineGermlineContigPloidyCaseMode {
        input:
            read_count_files = CollectCounts.counts,
            contig_ploidy_model_tar = contig_ploidy_model_tar,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            mem_gb = mem_gb_for_determine_germline_contig_ploidy,
            cpu = cpu_for_determine_germline_contig_ploidy,
            mapping_error_rate = ploidy_mapping_error_rate,
            sample_psi_scale = ploidy_sample_psi_scale,
            preemptible_attempts = preemptible_attempts
    }

    call CNVTasks.ScatterIntervals {
        input:
            interval_list = PreprocessIntervals.preprocessed_intervals,
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
                intervals = ScatterIntervals.scattered_interval_lists[scatter_index],
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

    call CNVTasks.PostprocessGermlineCNVCalls {
        input:
            entity_id = CollectCounts.entity_id,
            gcnv_calls_tars = GermlineCNVCallerCaseMode.gcnv_calls_tar,
            gcnv_model_tars = gcnv_model_tars,
            allosomal_contigs = allosomal_contigs,
            ref_copy_number_autosomal_contigs = ref_copy_number_autosomal_contigs,
            contig_ploidy_calls_tar = DetermineGermlineContigPloidyCaseMode.contig_ploidy_calls_tar,
            sample_index = 0,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            preemptible_attempts = preemptible_attempts
    }

    output {
        File preprocessed_intervals = PreprocessIntervals.preprocessed_intervals
        File read_counts_entity_id = CollectCounts.entity_id
        File read_counts = CollectCounts.counts
        File contig_ploidy_calls_tar = DetermineGermlineContigPloidyCaseMode.contig_ploidy_calls_tar
        Array[File] gcnv_calls_tars = GermlineCNVCallerCaseMode.gcnv_calls_tar
        Array[File] gcnv_tracking_tars = GermlineCNVCallerCaseMode.gcnv_tracking_tar
        File genotyped_intervals_vcf = PostprocessGermlineCNVCalls.genotyped_intervals_vcf
        File genotyped_segments_vcf = PostprocessGermlineCNVCalls.genotyped_segments_vcf
    }
}

task DetermineGermlineContigPloidyCaseMode {
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

    # We do not expose Hybrid ADVI parameters -- the default values are decent

    Int machine_mem_mb = select_first([mem_gb, 7]) * 1000
    Int command_mem_mb = machine_mem_mb - 500

    # If optional output_dir not specified, use "out"
    String output_dir_ = select_first([output_dir, "out"])

    command <<<
        set -e
        mkdir ${output_dir_}
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk4_jar_override}
        export MKL_NUM_THREADS=${default=8 cpu}
        export OMP_NUM_THREADS=${default=8 cpu}

        mkdir input-contig-ploidy-model
        tar xzf ${contig_ploidy_model_tar} -C input-contig-ploidy-model

        gatk --java-options "-Xmx${command_mem_mb}m" DetermineGermlineContigPloidy \
            --input ${sep=" --input " read_count_files} \
            --model input-contig-ploidy-model \
            --output ${output_dir_} \
            --output-prefix case \
            --verbosity DEBUG \
            --mapping-error-rate ${default="0.01" mapping_error_rate} \
            --sample-psi-scale ${default="0.0001" sample_psi_scale}

        tar czf case-contig-ploidy-calls.tar.gz -C ${output_dir_}/case-calls .
    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 150]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 8])
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        File contig_ploidy_calls_tar = "case-contig-ploidy-calls.tar.gz"
    }
}

task GermlineCNVCallerCaseMode {
    Int scatter_index
    Array[File] read_count_files
    File contig_ploidy_calls_tar
    File gcnv_model_tar
    File intervals
    File? annotated_intervals
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

    Int machine_mem_mb = select_first([mem_gb, 7]) * 1000
    Int command_mem_mb = machine_mem_mb - 500

    # If optional output_dir not specified, use "out"
    String output_dir_ = select_first([output_dir, "out"])

    command <<<
        set -e
        mkdir ${output_dir_}
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk4_jar_override}
        export MKL_NUM_THREADS=${default=8 cpu}
        export OMP_NUM_THREADS=${default=8 cpu}

        mkdir contig-ploidy-calls-dir
        tar xzf ${contig_ploidy_calls_tar} -C contig-ploidy-calls-dir

        mkdir gcnv-model
        tar xzf ${gcnv_model_tar} -C gcnv-model

        gatk --java-options "-Xmx${command_mem_mb}m"  GermlineCNVCaller \
            --run-mode CASE \
            --input ${sep=" --input " read_count_files} \
            --contig-ploidy-calls contig-ploidy-calls-dir \
            --model gcnv-model \
            --output ${output_dir_} \
            --output-prefix case \
            --verbosity DEBUG \
            --p-alt ${default="1e-6" p_alt} \
            --cnv-coherence-length ${default="10000.0" cnv_coherence_length} \
            --max-copy-number ${default="5" max_copy_number} \
            --mapping-error-rate ${default="0.01" mapping_error_rate} \
            --sample-psi-scale ${default="0.0001" sample_psi_scale} \
            --depth-correction-tau ${default="10000.0" depth_correction_tau} \
            --copy-number-posterior-expectation-mode ${default="HYBRID" copy_number_posterior_expectation_mode} \
            --active-class-padding-hybrid-mode ${default="50000" active_class_padding_hybrid_mode} \
            --learning-rate ${default="0.05" learning_rate} \
            --adamax-beta-1 ${default="0.9" adamax_beta_1} \
            --adamax-beta-2 ${default="0.99" adamax_beta_2} \
            --log-emission-samples-per-round ${default="50" log_emission_samples_per_round} \
            --log-emission-sampling-median-rel-error ${default="0.005" log_emission_sampling_median_rel_error} \
            --log-emission-sampling-rounds ${default="10" log_emission_sampling_rounds} \
            --max-advi-iter-first-epoch ${default="5000" max_advi_iter_first_epoch} \
            --max-advi-iter-subsequent-epochs ${default="100" max_advi_iter_subsequent_epochs} \
            --min-training-epochs ${default="10" min_training_epochs} \
            --max-training-epochs ${default="100" max_training_epochs} \
            --initial-temperature ${default="2.0" initial_temperature} \
            --num-thermal-advi-iters ${default="2500" num_thermal_advi_iters} \
            --convergence-snr-averaging-window ${default="500" convergence_snr_averaging_window} \
            --convergence-snr-trigger-threshold ${default="0.1" convergence_snr_trigger_threshold} \
            --convergence-snr-countdown-window ${default="10" convergence_snr_countdown_window} \
            --max-calling-iters ${default="10" max_calling_iters} \
            --caller-update-convergence-threshold ${default="0.001" caller_update_convergence_threshold} \
            --caller-internal-admixing-rate ${default="0.75" caller_internal_admixing_rate} \
            --caller-external-admixing-rate ${default="1.00" caller_external_admixing_rate} \
            --disable-annealing ${default="false" disable_annealing}

        tar czf case-gcnv-calls-${scatter_index}.tar.gz -C ${output_dir_}/case-calls .
        tar czf case-gcnv-tracking-${scatter_index}.tar.gz -C ${output_dir_}/case-tracking .
    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 150]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 8])
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        File gcnv_calls_tar = "case-gcnv-calls-${scatter_index}.tar.gz"
        File gcnv_tracking_tar = "case-gcnv-tracking-${scatter_index}.tar.gz"
    }
}
