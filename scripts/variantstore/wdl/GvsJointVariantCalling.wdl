version 1.0

import "GvsUnified.wdl" as GvsUnified
import "GvsBulkIngestGenomes.wdl" as BulkIngestGenomes
import "GvsPopulateAltAllele.wdl" as PopulateAltAllele
import "GvsCreateFilterSet.wdl" as CreateFilterSet
import "GvsPrepareRangesCallset.wdl" as PrepareRangesCallset
import "GvsExtractCallset.wdl" as ExtractCallset

workflow GvsJointVariantCalling {
    input {
        String dataset_name
        String project_id
        Array[String] external_sample_names
        Array[File] input_vcfs
        Array[File] input_vcf_indexes
        String call_set_identifier
        String samples_table_name = "sample"
        String sample_id_column_name = "sample_id"

        String? extract_output_gcs_dir
        Boolean is_beta_user = false
        String drop_state = "FORTY"
        Boolean use_classic_VQSR = true
        # Beta users have accounts with tighter quotas, and we must work around that
        Boolean tighter_gcp_quotas = true
    }

    # the call_set_identifier string is used to name many different things throughout this workflow (BQ tables, vcfs etc),
    # and so make sure nothing is broken by creative users, we replace spaces and underscores with hyphens
    String extract_output_file_base_name = sub(call_set_identifier, "\\s+|\_+", "-")
    String extract_table_prefix = sub(call_set_identifier, "\\s+|\_+", "-")
    String filter_set_name = sub(call_set_identifier, "\\s+|\_+", "-")
    if (false) {
        Int extract_maxretries_override = ""
        Int extract_preemptible_override = ""
        Int extract_scatter_count = ""
        Int load_data_batch_size = ""
        Int load_data_preemptible_override = ""
        Int load_data_maxretries_override = ""
        Array[String] query_labels = []
        File sample_names_to_extract = ""
        Int split_intervals_disk_size_override = ""
        Int split_intervals_mem_override = ""
        Int INDEL_VQSR_CLASSIC_max_gaussians_override = 4
        Int INDEL_VQSR_CLASSIC_mem_gb_override = ""
        Int SNP_VQSR_CLASSIC_max_gaussians_override = 6
        Int SNP_VQSR_CLASSIC_mem_gb_override = ""
    }
    # This is the most updated snapshot of the code as of June 2, 2023
    File gatk_override = "gs://gvs_quickstart_storage/jars/gatk-package-4.2.0.0-747-gbb08d81-SNAPSHOT-local.jar"
    File interval_list = "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list"

    File interval_weights_bed = "gs://broad-public-datasets/gvs/weights/gvs_vet_weights_1kb.bed"
    String fq_temp_table_dataset = "~{project_id}.~{dataset_name}"

    call BulkIngestGenomes.GvsBulkIngestGenomes as BulkIngestGenomes {
        input:
            dataset_name = dataset_name,
            project_id = project_id,
            call_set_identifier = call_set_identifier,
            samples_table_name = samples_table_name,
            sample_id_column_name = sample_id_column_name,
            interval_list = interval_list,
            drop_state = drop_state,
            gatk_override = gatk_override,
    }

    call PopulateAltAllele.GvsPopulateAltAllele {
        input:
            call_set_identifier = call_set_identifier,
            go = BulkIngestGenomes.done,
            dataset_name = dataset_name,
            project_id = project_id
    }

    call CreateFilterSet.GvsCreateFilterSet {
        input:
            go = GvsPopulateAltAllele.done,
            dataset_name = dataset_name,
            project_id = project_id,
            call_set_identifier = call_set_identifier,
            filter_set_name = filter_set_name,
            use_VQSR_lite = !use_classic_VQSR,
            interval_list = interval_list,
            gatk_override = gatk_override,
            INDEL_VQSR_CLASSIC_max_gaussians_override = INDEL_VQSR_CLASSIC_max_gaussians_override,
            INDEL_VQSR_CLASSIC_mem_gb_override = INDEL_VQSR_CLASSIC_mem_gb_override,
            SNP_VQSR_CLASSIC_max_gaussians_override = SNP_VQSR_CLASSIC_max_gaussians_override,
            SNP_VQSR_CLASSIC_mem_gb_override = SNP_VQSR_CLASSIC_mem_gb_override
    }

    call PrepareRangesCallset.GvsPrepareCallset {
        input:
            call_set_identifier = call_set_identifier,
            go = GvsCreateFilterSet.done,
            dataset_name = dataset_name,
            project_id = project_id,
            extract_table_prefix = extract_table_prefix,
            query_project = project_id,
            destination_project = project_id,
            destination_dataset = dataset_name,
            fq_temp_table_dataset = fq_temp_table_dataset,
            query_labels = query_labels,
            sample_names_to_extract = sample_names_to_extract
    }

    call ExtractCallset.GvsExtractCallset {
        input:
            go = GvsPrepareCallset.done,
            dataset_name = dataset_name,
            project_id = project_id,
            call_set_identifier = call_set_identifier,
            extract_table_prefix = extract_table_prefix,
            filter_set_name = filter_set_name,
            query_project = project_id,
            scatter_count = extract_scatter_count,
            interval_list = interval_list,
            interval_weights_bed = interval_weights_bed,
            gatk_override = gatk_override,
            output_file_base_name = extract_output_file_base_name,
            extract_maxretries_override = extract_maxretries_override,
            extract_preemptible_override = extract_preemptible_override,
            output_gcs_dir = extract_output_gcs_dir,
            split_intervals_disk_size_override = split_intervals_disk_size_override,
            split_intervals_mem_override = split_intervals_mem_override,
            do_not_filter_override = false,
            drop_state = drop_state
    }

    output {
        Array[File] output_vcfs = GvsExtractCallset.output_vcfs
        Array[File] output_vcf_indexes = GvsExtractCallset.output_vcf_indexes
        Array[File] output_vcf_interval_files = GvsExtractCallset.output_vcf_interval_files
        Float total_vcfs_size_mb = GvsExtractCallset.total_vcfs_size_mb
        File? sample_name_list = GvsExtractCallset.sample_name_list
        File manifest = GvsExtractCallset.manifest
    }
}
