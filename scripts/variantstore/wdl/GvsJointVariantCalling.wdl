version 1.0

import "GvsUnified.wdl" as GvsUnified

workflow GvsJointVariantCalling {
    input {
        String dataset_name
        String project_id
        Array[String] external_sample_names
        Array[File] input_vcfs
        Array[File] input_vcf_indexes
        String callset_identifier
        String? extract_output_gcs_dir
    }

    String extract_output_file_base_name = callset_identifier ## TODO make sure there are no spaces here!??!
    String extract_table_prefix = callset_identifier ## TODO make sure there are no spaces here!??!
    if(false) {
      Int extract_maxretries_override = ""
      Int extract_preemptible_override = ""
      Int extract_scatter_count = ""
      File gatk_override = ""
      Int load_data_preemptible_override = ""
      Int load_data_maxretries_override = ""
      Array[String] query_labels = []
      String service_account_json_path = ""
      File sample_names_to_extract = ""
      Int split_intervals_disk_size_override = ""
      Int split_intervals_mem_override = ""
      Int INDEL_VQSR_max_gaussians_override = 4
      Int INDEL_VQSR_mem_gb_override = ""
      Int SNP_VQSR_max_gaussians_override = 6
      Int SNP_VQSR_mem_gb_override = ""
    }
    File interval_list = "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list"
    Array[String] indel_recalibration_annotation_values = ["AS_FS", "AS_ReadPosRankSum", "AS_MQRankSum", "AS_QD", "AS_SOR"]
    Array[String] snp_recalibration_annotation_values = ["AS_QD", "AS_MQRankSum", "AS_ReadPosRankSum", "AS_FS", "AS_MQ", "AS_SOR"]
    File interval_weights_bed = "gs://broad-public-datasets/gvs/weights/gvs_vet_weights_1kb.bed"


    call GvsUnified.GvsUnified {
        input:
            dataset_name = dataset_name,
            project_id = project_id,
            external_sample_names = external_sample_names,
            input_vcf_indexes = input_vcf_indexes,
            input_vcfs = input_vcfs,
            filter_set_name = callset_identifier,
            extract_output_gcs_dir = extract_output_gcs_dir,
            destination_dataset = dataset_name,
            destination_project = project_id,
            extract_do_not_filter_override = false,
            extract_maxretries_override = extract_maxretries_override,
            extract_output_file_base_name = callset_identifier,
            extract_preemptible_override = extract_preemptible_override,
            extract_scatter_count = extract_scatter_count,
            extract_table_prefix = callset_identifier,
            fq_temp_table_dataset = "~{project_id}.~{dataset_name}",
            gatk_override = gatk_override,
            indel_recalibration_annotation_values = indel_recalibration_annotation_values,
            interval_list = interval_list,
            interval_weights_bed = interval_weights_bed,
            load_data_batch_size = 5,
            load_data_maxretries_override = load_data_maxretries_override,
            load_data_preemptible_override = load_data_preemptible_override,
            query_labels = query_labels,
            query_project = project_id,
            sample_names_to_extract = sample_names_to_extract,
            service_account_json_path = service_account_json_path,
            snp_recalibration_annotation_values = snp_recalibration_annotation_values,
            split_intervals_disk_size_override = split_intervals_disk_size_override,
            split_intervals_mem_override = split_intervals_mem_override,
            INDEL_VQSR_max_gaussians_override = INDEL_VQSR_max_gaussians_override,
            INDEL_VQSR_mem_gb_override = INDEL_VQSR_mem_gb_override,
            SNP_VQSR_max_gaussians_override = SNP_VQSR_max_gaussians_override,
            SNP_VQSR_mem_gb_override = SNP_VQSR_mem_gb_override,
    }

    output {
        Array[File] output_vcfs = GvsUnified.output_vcfs
        Array[File] output_vcf_indexes = GvsUnified.output_vcf_indexes
        Float total_vcfs_size_mb = GvsUnified.total_vcfs_size_mb
        File manifest = GvsUnified.manifest
        Boolean done = true
    }
}
