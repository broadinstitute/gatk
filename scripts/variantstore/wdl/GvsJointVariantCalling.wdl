version 1.0

import "GvsBulkIngestGenomes.wdl" as BulkIngestGenomes
import "GvsPopulateAltAllele.wdl" as PopulateAltAllele
import "GvsCreateFilterSet.wdl" as CreateFilterSet
import "GvsPrepareRangesCallset.wdl" as PrepareRangesCallset
import "GvsExtractCallset.wdl" as ExtractCallset

workflow GvsJointVariantCalling {
    input {
        String dataset_name
        String project_id
        String call_set_identifier
        String? extract_output_gcs_dir
        String drop_state = "FORTY"
        Boolean use_classic_VQSR = true
        # Beta users have accounts with tighter quotas, and we must work around that
        Boolean tighter_gcp_quotas = true
        String? sample_id_column_name ## Note that a column WILL exist that is the <entity>_id from the table name. However, some users will want to specify an alternate column for the sample_name during ingest
        String? vcf_files_column_name
        String? vcf_index_files_column_name
        String? sample_set_name ## NOTE: currently we only allow the loading of one sample set at a time

        # Docker images in order of increasing size / complexity, with the intent of using the smallest image suitable
        # for most jobs.
        String basic_docker = "ubuntu:22.04"
        String cloud_sdk_docker = "gcr.io/google.com/cloudsdktool/cloud-sdk:435.0.0-alpine"
        String variants_docker = "us.gcr.io/broad-dsde-methods/variantstore:2023-08-03-alpine-d9f94010b"

        # NOTE: `gatk_override` is not intended for production runs of the GVS pipeline! If defined, `gatk_override`
        # should be at least as recent as `gatk_docker`. Legitimate uses of `gatk_override` include integration test
        # runs and feature development.
        File? gatk_override
        # This is the most updated snapshot of GATK code on `bulk_ingest_staging` (off `ah_var_store`) as of 2023-08-03.
        String gatk_docker = "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:varstore_bulk_ingest_staging_2023_08_03"

        File interval_list = "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list"
        Boolean use_interval_weights = true
        File interval_weights_bed = "gs://broad-public-datasets/gvs/weights/gvs_vet_weights_1kb.bed"
        Boolean extract_do_not_filter_override = false
        String? extract_output_file_base_name
        String? extract_table_prefix
        String? filter_set_name

        # Overrides to be passed to GvsCreateFilterSet
        Int? INDEL_VQSR_CLASSIC_max_gaussians_override = 4
        Int? INDEL_VQSR_CLASSIC_mem_gb_override
        Int? SNP_VQSR_CLASSIC_max_gaussians_override = 6
        Int? SNP_VQSR_CLASSIC_mem_gb_override
    }

    # the call_set_identifier string is used to name many different things throughout this workflow (BQ tables, vcfs etc),
    # and so make sure nothing is broken by creative users, we replace spaces and underscores with hyphens
    String effective_extract_output_file_base_name = select_first([extract_output_file_base_name, sub(call_set_identifier, "\\s+|\_+", "-")])
    String effective_extract_table_prefix = select_first([extract_table_prefix, sub(call_set_identifier, "\\s+|\_+", "-")])
    String effective_filter_set_name = select_first([filter_set_name, sub(call_set_identifier, "\\s+|\_+", "-")])

    String query_project = project_id
    String destination_project = project_id
    String destination_dataset = dataset_name
    String fq_temp_table_dataset = "~{destination_project}.~{destination_dataset}"

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
    }

    call BulkIngestGenomes.GvsBulkIngestGenomes as BulkIngestGenomes {
        input:
            dataset_name = dataset_name,
            project_id = project_id,
            basic_docker = basic_docker,
            cloud_sdk_docker = cloud_sdk_docker,
            variants_docker = variants_docker,
            gatk_docker = gatk_docker,
            gatk_override = gatk_override,
            interval_list = interval_list,
            drop_state = drop_state,
            sample_id_column_name = sample_id_column_name,
            vcf_files_column_name = vcf_files_column_name,
            vcf_index_files_column_name = vcf_index_files_column_name,
            sample_set_name = sample_set_name,
    }

    call PopulateAltAllele.GvsPopulateAltAllele {
        input:
            call_set_identifier = call_set_identifier,
            go = BulkIngestGenomes.done,
            dataset_name = dataset_name,
            project_id = project_id,
            variants_docker = variants_docker,
            cloud_sdk_docker = cloud_sdk_docker,
    }

    call CreateFilterSet.GvsCreateFilterSet {
        input:
            go = GvsPopulateAltAllele.done,
            dataset_name = dataset_name,
            project_id = project_id,
            call_set_identifier = call_set_identifier,
            filter_set_name = effective_filter_set_name,
            use_VQSR_lite = !use_classic_VQSR,
            interval_list = interval_list,
            variants_docker = variants_docker,
            gatk_docker = gatk_docker,
            gatk_override = gatk_override,
            INDEL_VQSR_CLASSIC_max_gaussians_override = INDEL_VQSR_CLASSIC_max_gaussians_override,
            INDEL_VQSR_CLASSIC_mem_gb_override = INDEL_VQSR_CLASSIC_mem_gb_override,
            SNP_VQSR_CLASSIC_max_gaussians_override = SNP_VQSR_CLASSIC_max_gaussians_override,
            SNP_VQSR_CLASSIC_mem_gb_override = SNP_VQSR_CLASSIC_mem_gb_override,
            cloud_sdk_docker = cloud_sdk_docker,
    }

    call PrepareRangesCallset.GvsPrepareCallset {
        input:
            call_set_identifier = call_set_identifier,
            go = GvsCreateFilterSet.done,
            dataset_name = dataset_name,
            project_id = project_id,
            extract_table_prefix = effective_extract_table_prefix,
            query_project = query_project,
            destination_project = destination_project,
            destination_dataset = destination_dataset,
            fq_temp_table_dataset = fq_temp_table_dataset,
            query_labels = query_labels,
            sample_names_to_extract = sample_names_to_extract,
            variants_docker = variants_docker,
    }

    call ExtractCallset.GvsExtractCallset {
        input:
            go = GvsPrepareCallset.done,
            dataset_name = dataset_name,
            project_id = project_id,
            call_set_identifier = call_set_identifier,
            extract_table_prefix = effective_extract_table_prefix,
            filter_set_name = effective_filter_set_name,
            query_project = query_project,
            scatter_count = extract_scatter_count,
            interval_list = interval_list,
            use_interval_weights = use_interval_weights,
            interval_weights_bed = interval_weights_bed,
            variants_docker = variants_docker,
            gatk_docker = gatk_docker,
            gatk_override = gatk_override,
            cloud_sdk_docker = cloud_sdk_docker,
            output_file_base_name = effective_extract_output_file_base_name,
            extract_maxretries_override = extract_maxretries_override,
            extract_preemptible_override = extract_preemptible_override,
            output_gcs_dir = extract_output_gcs_dir,
            split_intervals_disk_size_override = split_intervals_disk_size_override,
            split_intervals_mem_override = split_intervals_mem_override,
            do_not_filter_override = extract_do_not_filter_override,
            drop_state = drop_state,
    }

    output {
        Array[File] output_vcfs = GvsExtractCallset.output_vcfs
        Array[File] output_vcf_indexes = GvsExtractCallset.output_vcf_indexes
        Array[File] output_vcf_interval_files = GvsExtractCallset.output_vcf_interval_files
        Float total_vcfs_size_mb = GvsExtractCallset.total_vcfs_size_mb
        File? sample_name_list = GvsExtractCallset.sample_name_list
        File manifest = GvsExtractCallset.manifest
        Boolean done = true
    }
}
