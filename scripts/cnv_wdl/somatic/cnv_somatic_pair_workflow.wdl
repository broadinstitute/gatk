# Workflow for running the GATK CNV pipeline on a matched pair. Supports both WGS and WES.
#
# Notes:
#
# - The interval-list file is required for both WGS and WES workflows and should be a Picard or GATK-style interval list.
#   These intervals will be padded on both sides by the amount specified by PreprocessIntervals.padding (default 250)
#   and split into bins of length specified by PreprocessIntervals.bin_length (default 1000; specify 0 to skip binning).
#   For WGS, the intervals should simply cover the autosomal chromosomes (sex chromosomes may be included, but care
#   should be taken to 1) avoid creating panels of mixed sex, and 2) denoise case samples only with panels containing
#   individuals of the same sex as the case samples).
#
# - The sites file (common_sites) should be a Picard or GATK-style interval list.  This is a list of sites
#   of known variation at which allelic counts will be collected for use in modeling minor-allele fractions.
#
# - Example invocation:
#    java -jar cromwell.jar run cnv_somatic_pair_workflow.wdl myParameters.json
#   See cnv_somatic_pair_workflow_template.json for a template json file to modify with your own parameters (please save
#   your modified version with a different filename and do not commit to the gatk repository).
#
#############

import "cnv_common_tasks.wdl" as CNVTasks

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
        GATK_JAR=${default="/root/gatk.jar" gatk4_jar_override}

        java -Xmx${command_mem}m -jar $GATK_JAR DenoiseReadCounts \
            --input ${read_counts} \
            --readCountPanelOfNormals ${read_count_pon} \
            ${"--numberOfEigensamples " + number_of_eigensamples} \
            --standardizedCopyRatios ${entity_id}.standardizedCR.tsv \
            --denoisedCopyRatios ${entity_id}.denoisedCR.tsv
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
        GATK_JAR=${default="/root/gatk.jar" gatk4_jar_override}

        java -Xmx${command_mem}m -jar $GATK_JAR ModelSegments \
            --denoisedCopyRatios ${denoised_copy_ratios} \
            --allelicCounts ${allelic_counts} \
            ${"--normalAllelicCounts " + normal_allelic_counts} \
            --maxNumSegmentsPerChromosome ${default="500" max_num_segments_per_chromosome} \
            --minTotalAlleleCount ${default="30" min_total_allele_count} \
            --genotypingHomozygousLogRatioThreshold ${default="-10.0" genotyping_homozygous_log_ratio_threshold} \
            --genotypingBaseErrorRate ${default="0.05" genotyping_base_error_rate} \
            --kernelVarianceCopyRatio ${default="0.0" kernel_variance_copy_ratio} \
            --kernelVarianceAlleleFraction ${default="0.025" kernel_variance_allele_fraction} \
            --kernelScalingAlleleFraction ${default="1.0" kernel_scaling_allele_fraction} \
            --kernelApproximationDimension ${default="100" kernel_approximation_dimension} \
            --windowSize ${sep= " --windowSize " window_sizes} \
            --numChangepointsPenaltyFactor ${default="1.0" num_changepoints_penalty_factor} \
            --minorAlleleFractionPriorAlpha ${default="25.0" minor_allele_fraction_prior_alpha} \
            --numSamplesCopyRatio ${default=100 num_samples_copy_ratio} \
            --numBurnInCopyRatio ${default=50 num_burn_in_copy_ratio} \
            --numSamplesAlleleFraction ${default=100 num_samples_allele_fraction} \
            --numBurnInAlleleFraction ${default=50 num_burn_in_allele_fraction} \
            --smoothingThresholdCopyRatio ${default="2.0" smoothing_threshold_copy_ratio} \
            --smoothingThresholdAlleleFraction ${default="2.0" smoothing_threshold_allele_fraction} \
            --maxNumSmoothingIterations ${default=10 max_num_smoothing_iterations} \
            --numSmoothingIterationsPerFit ${default=0 num_smoothing_iterations_per_fit} \
            --output ${output_dir_} \
            --outputPrefix ${entity_id}

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
        GATK_JAR=${default="/root/gatk.jar" gatk4_jar_override}

        java -Xmx${command_mem}m -jar $GATK_JAR CallCopyRatioSegments \
            --input ${copy_ratio_segments} \
            --neutralSegmentCopyRatioThreshold ${default="0.1" neutral_segment_copy_ratio_threshold} \
            --outlierNeutralSegmentCopyRatioZScoreThreshold ${default="2.0" outlier_neutral_segment_copy_ratio_z_score_threshold} \
            --callingCopyRatioZScoreThreshold ${default="2.0" calling_copy_ratio_z_score_threshold} \
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
        GATK_JAR=${default="/root/gatk.jar" gatk4_jar_override}

        java -Xmx${command_mem}m -jar $GATK_JAR PlotDenoisedCopyRatios \
            --standardizedCopyRatios ${standardized_copy_ratios} \
            --denoisedCopyRatios ${denoised_copy_ratios} \
            -SD ${ref_fasta_dict} \
            --minimumContigLength ${default="1000000" minimum_contig_length} \
            --output ${output_dir_} \
            --outputPrefix ${entity_id}
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
        GATK_JAR=${default="/root/gatk.jar" gatk4_jar_override}

        java -Xmx${command_mem}m -jar $GATK_JAR PlotModeledSegments \
            --denoisedCopyRatios ${denoised_copy_ratios} \
            --allelicCounts ${het_allelic_counts} \
            --segments ${modeled_segments} \
            -SD ${ref_fasta_dict} \
            --minimumContigLength ${default="1000000" minimum_contig_length} \
            --output ${output_dir_} \
            --outputPrefix ${entity_id}
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
