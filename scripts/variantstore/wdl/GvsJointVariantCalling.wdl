version 1.0

import "GvsBulkIngestGenomes.wdl" as BulkIngestGenomes
import "GvsPopulateAltAllele.wdl" as PopulateAltAllele
import "GvsCreateFilterSet.wdl" as CreateFilterSet
import "GvsPrepareRangesCallset.wdl" as PrepareRangesCallset
import "GvsExtractCallset.wdl" as ExtractCallset

workflow GvsJointVariantCalling {
    input {
        String call_set_identifier
        String dataset_name
        String project_id

        String drop_state = "FORTY"
        Boolean extract_do_not_filter_override = false
        String? extract_output_gcs_dir
        File interval_list = "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list"
        Boolean is_beta_user = false
        String sample_id_column_name = "sample_id"
        String samples_table_name = "sample"
        Boolean tighter_gcp_quotas = true
        Boolean use_classic_VQSR = true
        String? vcf_files_column_name
        String? vcf_index_files_column_name

        # This is the most updated snapshot of the code as of June 22, 2023
        File gatk_override = "gs://gvs_quickstart_storage/jars/gatk-package-4.2.0.0-726-g63a0bea-SNAPSHOT-local.jar"
    }

    # the call_set_identifier string is used to name many different things throughout this workflow (BQ tables, vcfs etc),
    # and so make sure nothing is broken by creative users, we replace spaces and underscores with hyphens
    String extract_output_file_base_name = sub(call_set_identifier, "\\s+|\_+", "-")
    String extract_table_prefix = sub(call_set_identifier, "\\s+|\_+", "-")
    String filter_set_name = sub(call_set_identifier, "\\s+|\_+", "-")

    File interval_weights_bed = "gs://broad-public-datasets/gvs/weights/gvs_vet_weights_1kb.bed"
    String fq_temp_table_dataset = "~{project_id}.~{dataset_name}"

    call BulkIngestGenomes.GvsBulkIngestGenomes as BulkIngestGenomes {
        input:
            call_set_identifier = call_set_identifier,
            dataset_name = dataset_name,
            drop_state = drop_state,
            gatk_override = gatk_override,
            interval_list = interval_list,
            project_id = project_id,
            sample_id_column_name = sample_id_column_name,
            samples_table_name = samples_table_name,
            vcf_files_column_name = vcf_files_column_name,
            vcf_index_files_column_name = vcf_index_files_column_name,
    }

    call PopulateAltAllele.GvsPopulateAltAllele {
        input:
            call_set_identifier = call_set_identifier,
            dataset_name = dataset_name,
            go = BulkIngestGenomes.done,
            project_id = project_id,
    }

    call CreateFilterSet.GvsCreateFilterSet {
        input:
            call_set_identifier = call_set_identifier,
            dataset_name = dataset_name,
            filter_set_name = filter_set_name,
            gatk_override = gatk_override,
            go = GvsPopulateAltAllele.done,
            INDEL_VQSR_CLASSIC_max_gaussians_override = 4,
            interval_list = interval_list,
            project_id = project_id,
            SNP_VQSR_CLASSIC_max_gaussians_override = 6,
            use_VQSR_lite = !use_classic_VQSR,
    }

    call PrepareRangesCallset.GvsPrepareCallset {
        input:
            call_set_identifier = call_set_identifier,
            dataset_name = dataset_name,
            destination_dataset = dataset_name,
            destination_project = project_id,
            extract_table_prefix = extract_table_prefix,
            fq_temp_table_dataset = fq_temp_table_dataset,
            go = GvsPopulateAltAllele.done,
            project_id = project_id,
            query_project = project_id,
    }

    call ExtractCallset.GvsExtractCallset {
        input:
            call_set_identifier = call_set_identifier,
            dataset_name = dataset_name,
            do_not_filter_override = extract_do_not_filter_override,
            drop_state = drop_state,
            extract_table_prefix = extract_table_prefix,
            filter_set_name = GvsCreateFilterSet.loaded_filter_set_name,
            gatk_override = gatk_override,
            go = GvsPrepareCallset.done,
            interval_list = interval_list,
            interval_weights_bed = interval_weights_bed,
            output_file_base_name = extract_output_file_base_name,
            output_gcs_dir = extract_output_gcs_dir,
            project_id = project_id,
            query_project = project_id,
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
