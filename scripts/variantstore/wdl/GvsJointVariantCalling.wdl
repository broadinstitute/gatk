version 1.0

import "GvsBulkIngestGenomes.wdl" as BulkIngestGenomes
import "GvsPopulateAltAllele.wdl" as PopulateAltAllele
import "GvsCreateFilterSet.wdl" as CreateFilterSet
import "GvsPrepareRangesCallset.wdl" as PrepareRangesCallset
import "GvsExtractCallset.wdl" as ExtractCallset
import "GvsUtils.wdl" as Utils

workflow GvsJointVariantCalling {
    input {
        # If `git_branch_or_tag` is not specified by a caller (i.e. integration tests), default to the current
        # beta / early access tag. Note that setting this value does not "summon" this version of the workflow; only the
        # version selector in Terra can do that. The value of this `git_branch_or_tag` input should always be
        # consistent with the git tag assigned to its version.
        String git_branch_or_tag = "gvs_beta_0.2.0"
        # Potentially specified by a calling integration WDL.
        String? workflow_git_hash

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

        String? basic_docker
        String? cloud_sdk_docker
        String? variants_docker
        String? gatk_docker

        # NOTE: `gatk_override` is not intended for production runs of the GVS pipeline! If defined, `gatk_override`
        # should be at least as recent as `gatk_docker`. Legitimate uses of `gatk_override` include integration test
        # runs and feature development.
        File? gatk_override

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

    if (!defined(workflow_git_hash) || !defined(basic_docker) || !defined(cloud_sdk_docker) || !defined(variants_docker) || !defined(gatk_docker)) {
        call Utils.GetToolVersions {
            input:
                git_branch_or_tag = git_branch_or_tag,
        }
    }

    String effective_basic_docker = select_first([basic_docker, GetToolVersions.basic_docker])
    String effective_cloud_sdk_docker = select_first([cloud_sdk_docker, GetToolVersions.cloud_sdk_docker])
    String effective_variants_docker = select_first([variants_docker, GetToolVersions.variants_docker])
    String effective_gatk_docker = select_first([gatk_docker, GetToolVersions.gatk_docker])
    String effective_workflow_git_hash = select_first([workflow_git_hash, GetToolVersions.workflow_git_hash])

    call BulkIngestGenomes.GvsBulkIngestGenomes as BulkIngestGenomes {
        input:
            git_branch_or_tag = git_branch_or_tag,
            workflow_git_hash = effective_workflow_git_hash,
            dataset_name = dataset_name,
            project_id = project_id,
            basic_docker = effective_basic_docker,
            cloud_sdk_docker = effective_cloud_sdk_docker,
            variants_docker = effective_variants_docker,
            gatk_docker = effective_gatk_docker,
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
            git_branch_or_tag = git_branch_or_tag,
            workflow_git_hash = effective_workflow_git_hash,
            call_set_identifier = call_set_identifier,
            go = BulkIngestGenomes.done,
            dataset_name = dataset_name,
            project_id = project_id,
            variants_docker = effective_variants_docker,
            cloud_sdk_docker = effective_cloud_sdk_docker,
    }

    call CreateFilterSet.GvsCreateFilterSet {
        input:
            go = GvsPopulateAltAllele.done,
            git_branch_or_tag = git_branch_or_tag,
            workflow_git_hash = effective_workflow_git_hash,
            dataset_name = dataset_name,
            project_id = project_id,
            call_set_identifier = call_set_identifier,
            filter_set_name = effective_filter_set_name,
            use_VQSR_lite = !use_classic_VQSR,
            interval_list = interval_list,
            variants_docker = effective_variants_docker,
            gatk_docker = effective_gatk_docker,
            gatk_override = gatk_override,
            INDEL_VQSR_CLASSIC_max_gaussians_override = INDEL_VQSR_CLASSIC_max_gaussians_override,
            INDEL_VQSR_CLASSIC_mem_gb_override = INDEL_VQSR_CLASSIC_mem_gb_override,
            SNP_VQSR_CLASSIC_max_gaussians_override = SNP_VQSR_CLASSIC_max_gaussians_override,
            SNP_VQSR_CLASSIC_mem_gb_override = SNP_VQSR_CLASSIC_mem_gb_override,
            cloud_sdk_docker = effective_cloud_sdk_docker,
    }

    call PrepareRangesCallset.GvsPrepareCallset {
        input:
            call_set_identifier = call_set_identifier,
            go = GvsCreateFilterSet.done,
            git_branch_or_tag = git_branch_or_tag,
            workflow_git_hash = effective_workflow_git_hash,
            dataset_name = dataset_name,
            project_id = project_id,
            extract_table_prefix = effective_extract_table_prefix,
            query_project = query_project,
            destination_project = destination_project,
            destination_dataset = destination_dataset,
            fq_temp_table_dataset = fq_temp_table_dataset,
            query_labels = query_labels,
            sample_names_to_extract = sample_names_to_extract,
            variants_docker = effective_variants_docker,
    }

    call ExtractCallset.GvsExtractCallset {
        input:
            go = GvsPrepareCallset.done,
            git_branch_or_tag = git_branch_or_tag,
            workflow_git_hash = effective_workflow_git_hash,
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
            variants_docker = effective_variants_docker,
            gatk_docker = effective_gatk_docker,
            gatk_override = gatk_override,
            cloud_sdk_docker = effective_cloud_sdk_docker,
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
        String recorded_workflow_git_hash = effective_workflow_git_hash
        Boolean done = true
    }
}
