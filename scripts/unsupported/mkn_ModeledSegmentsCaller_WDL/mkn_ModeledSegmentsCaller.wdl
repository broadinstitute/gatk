
workflow ModeledSegmentsCallerWorkflow {
    String entity_id
    File modeled_segments_input_file
    Boolean? load_copy_ratio
    Boolean? load_allele_fraction
    File? output_dir
    String? output_prefix
    Float? normal_minor_allele_fraction_threshold
    Float? copy_ratio_peak_min_relative_height
    Float? copy_ratio_kernel_density_bandwidth
    Float? min_fraction_of_points_in_normal_allele_fraction_region
    Float? min_weight_first_cr_peak_cr_data_only

    File? gatk4_jar_override

    # Runtime parameters
    String gatk_docker
    Int? mem_gb
    Int? disk_space_gb
    Int? cpu
    Int? preemptible_attempts

    call CallModeledSegmentsTest {
        input:
            entity_id = entity_id,
            modeled_segments_input_file = modeled_segments_input_file,
            load_copy_ratio = load_copy_ratio,
            load_allele_fraction = load_allele_fraction,
            normal_minor_allele_fraction_threshold = normal_minor_allele_fraction_threshold,
            copy_ratio_peak_min_relative_height = copy_ratio_peak_min_relative_height,
            copy_ratio_kernel_density_bandwidth = copy_ratio_kernel_density_bandwidth,
            min_weight_first_cr_peak_cr_data_only = min_weight_first_cr_peak_cr_data_only,
            min_fraction_of_points_in_normal_allele_fraction_region = min_fraction_of_points_in_normal_allele_fraction_region,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            mem_gb = mem_gb,
            disk_space_gb = disk_space_gb,
            preemptible_attempts = preemptible_attempts,
            output_dir = output_dir,
            output_prefix = output_prefix,
            cpu=cpu
    }

}

task CallModeledSegmentsTest {
    String entity_id
    File modeled_segments_input_file
    Boolean? load_copy_ratio
    Boolean? load_allele_fraction
    File? output_dir
    String? output_prefix
    Float? normal_minor_allele_fraction_threshold
    Float? copy_ratio_peak_min_relative_height
    Float? copy_ratio_kernel_density_bandwidth
    Float? min_fraction_of_points_in_normal_allele_fraction_region
    Float? min_weight_first_cr_peak_cr_data_only
    File? gatk4_jar_override

    # Runtime parameters
    String gatk_docker
    Int? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? cpu
    Int? preemptible_attempts

    Int machine_mem_mb = select_first([mem_gb, 7]) * 1000
    Int command_mem_mb = machine_mem_mb - 1000

    String output_dir_ = select_first([output_dir, "out/"])
    String output_prefix_ = select_first([output_prefix, entity_id])

    command <<<
        set -e
        mkdir ${output_dir_}
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk4_jar_override}

        gatk --java-options "-Xmx${command_mem_mb}m" CallModeledSegments \
            --input ${modeled_segments_input_file} \
            --load-copy-ratio ${default="true" load_copy_ratio} \
            --load-allele-fraction ${default="true" load_allele_fraction} \
            --output ${output_dir_} \
            --output-prefix ${output_prefix_} \
            --normal-minor-allele-fraction-threshold ${default="0.475" normal_minor_allele_fraction_threshold} \
            --copy-ratio-peak-min-relative-height ${default="0.03" copy_ratio_peak_min_relative_height} \
            --copy-ratio-kernel-density-bandwidth ${default="0.0" copy_ratio_kernel_density_bandwidth} \
            --min-weight-first-cr-peak-cr-data-only ${default="0.35" min_weight_first_cr_peak_cr_data_only} \
            --min-fraction-of-points-in-normal-allele-fraction-region ${default="0.15" min_fraction_of_points_in_normal_allele_fraction_region}
    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 100]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File called_modeled_segments_data = "${output_dir_}${output_prefix_}.called.seg"
    }
}