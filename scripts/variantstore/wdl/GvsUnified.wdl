version 1.0

import "GvsAssignIds.wdl" as AssignIds
import "GvsImportGenomes.wdl" as ImportGenomes
import "GvsCreateAltAllele.wdl" as CreateAltAllele
import "GvsCreateFilterSet.wdl" as CreateFilterSet
import "GvsPrepareRangesCallset.wdl" as PrepareRangesCallset
import "GvsExtractCallset.wdl" as ExtractCallset

workflow GvsUnified {
    input {
        # Begin GvsAssignIds
        String dataset_name
        String project_id

        Array[String] external_sample_names
        Boolean samples_are_controls = false

        File? gatk_override
        String? service_account_json_path
        # End GvsAssignIds

        # Begin GvsImportGenomes
        Array[File] input_vcfs
        Array[File] input_vcf_indexes
        File interval_list = "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list"

        Int? load_data_preemptible_override
        Int? load_data_maxretries_override
        # End GvsImportGenomes

        # Begin GvsCreateFilterSet
        String filter_set_name
        Array[String] indel_recalibration_annotation_values = ["AS_FS", "AS_ReadPosRankSum", "AS_MQRankSum", "AS_QD", "AS_SOR"]
        Int create_filter_set_scatter_count
        Array[String] snp_recalibration_annotation_values = ["AS_QD", "AS_MQRankSum", "AS_ReadPosRankSum", "AS_FS", "AS_MQ", "AS_SOR"]

        Int? INDEL_VQSR_max_gaussians_override = 4
        Int? INDEL_VQSR_mem_gb_override
        Int? SNP_VQSR_max_gaussians_override = 6
        Int? SNP_VQSR_mem_gb_override
        # End GvsCreateFilterSet

        # Begin GvsPrepareRangesCallset
        String extract_table_prefix

        String query_project = project_id
        String destination_project = project_id
        String destination_dataset = dataset_name
        String fq_temp_table_dataset = "~{destination_project}.~{destination_dataset}"

        Array[String]? query_labels
        File? sample_names_to_extract
        # End GvsPrepareRangesCallset

        # Begin GvsExtractCallset
        Int extract_scatter_count

        File interval_weights_bed = "gs://broad-public-datasets/gvs/weights/gvs_vet_weights_1kb.bed"

        String extract_output_file_base_name = filter_set_name

        Int? extract_maxretries_override
        Int? extract_preemptible_override
        String? extract_output_gcs_dir
        Int? split_intervals_disk_size_override
        Int? split_intervals_mem_override
        # End GvsExtractCallset
    }

    call AssignIds.GvsAssignIds as AssignIds {
        input:
            dataset_name = dataset_name,
            project_id = project_id,
            external_sample_names = external_sample_names,
            samples_are_controls = samples_are_controls,
            assign_ids_gatk_override = gatk_override,
            service_account_json_path = service_account_json_path
    }

    call ImportGenomes.GvsImportGenomes {
        input:
            go = AssignIds.done,
            dataset_name = dataset_name,
            project_id = project_id,
            external_sample_names = external_sample_names,
            input_vcfs = input_vcfs,
            input_vcf_indexes = input_vcf_indexes,
            interval_list = interval_list,
            load_data_preemptible_override = load_data_preemptible_override,
            load_data_maxretries_override = load_data_maxretries_override,
            load_data_gatk_override = gatk_override,
            service_account_json_path = service_account_json_path
    }

    call CreateAltAllele.GvsCreateAltAllele {
        input:
            go = GvsImportGenomes.done,
            dataset_name = dataset_name,
            project_id = project_id,
            service_account_json_path = service_account_json_path
    }

    call CreateFilterSet.GvsCreateFilterSet {
        input:
            go = GvsCreateAltAllele.done,
            dataset_name = dataset_name,
            project_id = project_id,
            filter_set_name = filter_set_name,
            indel_recalibration_annotation_values = indel_recalibration_annotation_values,
            scatter_count = create_filter_set_scatter_count,
            snp_recalibration_annotation_values = snp_recalibration_annotation_values,
            interval_list = interval_list,
            gatk_override = gatk_override,
            INDEL_VQSR_max_gaussians_override = INDEL_VQSR_max_gaussians_override,
            INDEL_VQSR_mem_gb_override = INDEL_VQSR_mem_gb_override,
            service_account_json_path = service_account_json_path,
            SNP_VQSR_max_gaussians_override = SNP_VQSR_max_gaussians_override,
            SNP_VQSR_mem_gb_override = SNP_VQSR_mem_gb_override
    }

    call PrepareRangesCallset.GvsPrepareCallset {
        input:
            go = GvsCreateFilterSet.done,
            project_id = project_id,
            dataset_name = dataset_name,
            control_samples = samples_are_controls,
            extract_table_prefix = extract_table_prefix,
            query_project = query_project,
            destination_project = destination_project,
            destination_dataset = destination_dataset,
            fq_temp_table_dataset = fq_temp_table_dataset,
            query_labels = query_labels,
            sample_names_to_extract = sample_names_to_extract,
            service_account_json_path = service_account_json_path
    }

    call ExtractCallset.GvsExtractCallset {
        input:
            go = GvsPrepareCallset.done,
            dataset_name = dataset_name,
            project_id = project_id,
            control_samples = samples_are_controls,
            extract_table_prefix = extract_table_prefix,
            filter_set_name = filter_set_name,
            query_project = query_project,
            scatter_count = extract_scatter_count,
            interval_list = interval_list,
            interval_weights_bed = interval_weights_bed,
            gatk_override = gatk_override,
            output_file_base_name = extract_output_file_base_name,
            extract_maxretries_override = extract_maxretries_override,
            extract_preemptible_override = extract_preemptible_override,
            output_gcs_dir = extract_output_gcs_dir,
            service_account_json_path = service_account_json_path,
            split_intervals_disk_size_override = split_intervals_disk_size_override,
            split_intervals_mem_override = split_intervals_mem_override
    }

    output {
        Array[File] output_vcfs = GvsExtractCallset.output_vcfs
        Array[File] output_vcf_indexes = GvsExtractCallset.output_vcf_indexes
        Float total_vcfs_size_mb = GvsExtractCallset.total_vcfs_size_mb
        File manifest = GvsExtractCallset.manifest
        Boolean done = true
    }
}
