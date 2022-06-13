version 1.0

import "GvsUnified.wdl" as Unified

workflow GvsVariantCalling {
    input {
        String dataset_name
        String project_id
        Array[String] external_sample_names
        Array[File] input_vcfs
        Array[File] input_vcf_indexes
        String callset_identifier
    }

    Boolean samples_are_controls = false ## would setting this to True ever be useful? We wouldn't care about getting P&S on a filter model based on controls, right?
    File interval_list = "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list"
    Int load_data_batch_size = 5
    String filter_set_name = callset_identifier
    Array[String] indel_recalibration_annotation_values = ["AS_FS", "AS_ReadPosRankSum", "AS_MQRankSum", "AS_QD", "AS_SOR"]
    Array[String] snp_recalibration_annotation_values = ["AS_QD", "AS_MQRankSum", "AS_ReadPosRankSum", "AS_FS", "AS_MQ", "AS_SOR"]
    Int? INDEL_VQSR_max_gaussians_override = 4
    Int? SNP_VQSR_max_gaussians_override = 6
    String query_project = project_id
    String destination_project = project_id
    String destination_dataset = dataset_name
    String fq_temp_table_dataset = "~{destination_project}.~{destination_dataset}"
    File interval_weights_bed = "gs://broad-public-datasets/gvs/weights/gvs_vet_weights_1kb.bed"
    String extract_output_file_base_name = filter_set_name ## TODO make sure there are no spaces here!??!
    String extract_table_prefix = filter_set_name ## TODO make sure there are no spaces here!??!
    Boolean extract_do_not_filter_override = false

    call Unified.GvsUnified as Unified {
        input:
            dataset_name = dataset_name,
            project_id = project_id,
            external_sample_names = external_sample_names,
            samples_are_controls = samples_are_controls
            input_vcfs = input_vcfs,
            input_vcf_indexes = input_vcf_indexes,
            interval_list = interval_list,
            load_data_batch_size = load_data_batch_size
            filter_set_name = filter_set_name,
            indel_recalibration_annotation_values = indel_recalibration_annotation_values,
            snp_recalibration_annotation_values = snp_recalibration_annotation_values,
            extract_table_prefix = extract_table_prefix,
            query_project = query_project,
            destination_project = destination_project,
            destination_dataset = destination_dataset,
            fq_temp_table_dataset = fq_temp_table_dataset,
            interval_weights_bed = interval_weights_bed,
            extract_output_file_base_name = extract_output_file_base_name,
            extract_do_not_filter_override = extract_do_not_filter_override
    }

    output {
        Array[File] output_vcfs = Unified.output_vcfs
        Array[File] output_vcf_indexes = Unified.output_vcf_indexes
        Float total_vcfs_size_mb = Unified.total_vcfs_size_mb
        File manifest = Unified.manifest
    }
}
