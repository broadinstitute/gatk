version 1.0

import "GvsBulkIngestGenomes.wdl" as BulkIngestGenomes
import "GvsPopulateAltAllele.wdl" as PopulateAltAllele
import "GvsCreateFilterSet.wdl" as CreateFilterSet
import "GvsPrepareRangesCallset.wdl" as PrepareRangesCallset
import "GvsExtractCallset.wdl" as ExtractCallset
import "GvsUtils.wdl" as Utils
import "PrepareReferenceFiles.wdl" as PrepareReferenceFiles

# Here we go again. 6
workflow GvsJointVariantCalling {
    input {
        Boolean go = true
        String call_set_identifier
        String dataset_name
        String extract_output_gcs_dir
        String project_id

        String sample_id_column_name
        String vcf_files_column_name
        String vcf_index_files_column_name

        Boolean bgzip_output_vcfs = false
        Boolean collect_variant_calling_metrics = false
        String drop_state = "FORTY"
        Boolean use_VQSR = false
        Boolean use_compressed_references = false
        Boolean load_vet_and_ref_ranges = true
        Boolean load_vcf_headers = false
        # Beta users have accounts with tighter quotas, and we must work around that
        Boolean tighter_gcp_quotas = true
        String? sample_set_name ## NOTE: currently we only allow the loading of one sample set at a time
        String? billing_project_id

        # If `git_branch_or_tag` is not specified by a caller (i.e. integration tests), default to the current branch or
        # tag as specified in `GvsUtils.GetToolVersions`.
        String? git_branch_or_tag
        # Potentially specified by a calling integration WDL.
        String? git_hash

        String? basic_docker
        String? cloud_sdk_docker
        String? variants_docker
        String? gatk_docker

        String? workspace_bucket
        String? workspace_id
        String? submission_id

        # NOTE: `gatk_override` is not intended for production runs of the GVS pipeline! If defined, `gatk_override`
        # should be at least as recent as `gatk_docker`. Legitimate uses of `gatk_override` include integration test
        # runs and feature development.
        File? gatk_override

        String reference_name = "hg38"

        # for supporting custom references... for now. Later map the references and use the reference_name above
        File? reference

        Boolean is_wgs = true
        File? interval_list

        Boolean extract_do_not_filter_override = false
        String? extract_output_file_base_name
        String? extract_table_prefix
        String? filter_set_name

        File? target_interval_list

        # Overrides to be passed to GvsCreateFilterSet
        Int? INDEL_VQSR_max_gaussians_override = 4
        Int? INDEL_VQSR_mem_gb_override
        Int? SNP_VQSR_max_gaussians_override = 6
        Int? SNP_VQSR_mem_gb_override

        File? training_python_script
        File? scoring_python_script

        Int? maximum_alternate_alleles
    }

    # The `call_set_identifier` string is used to name many different things throughout this workflow (BQ tables, vcfs etc).
    # To make sure nothing is broken by creative users we replace spaces and underscores with hyphens.
    String effective_extract_output_file_base_name = select_first([extract_output_file_base_name, sub(call_set_identifier, "\\s+|\_+", "-")])
    String effective_extract_table_prefix = select_first([extract_table_prefix, sub(call_set_identifier, "\\s+|\_+", "-")])
    String effective_filter_set_name = select_first([filter_set_name, sub(call_set_identifier, "\\s+|\_+", "-")])

    String query_project = project_id
    String destination_project = project_id
    String destination_dataset = dataset_name
    String fq_temp_table_dataset = "~{destination_project}.~{destination_dataset}"
    String ploidy_table_name = "sample_chromosome_ploidy"

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

    if (!defined(git_hash) ||
        !defined(basic_docker) || !defined(cloud_sdk_docker) || !defined(variants_docker) || !defined(gatk_docker) ||
        !defined(workspace_bucket) || !defined(submission_id)) {
        call Utils.GetToolVersions {
            input:
                git_branch_or_tag = git_branch_or_tag,
        }
    }

    String effective_basic_docker = select_first([basic_docker, GetToolVersions.basic_docker])
    String effective_cloud_sdk_docker = select_first([cloud_sdk_docker, GetToolVersions.cloud_sdk_docker])
    String effective_variants_docker = select_first([variants_docker, GetToolVersions.variants_docker])
    String effective_gatk_docker = select_first([gatk_docker, GetToolVersions.gatk_docker])
    String effective_git_hash = select_first([git_hash, GetToolVersions.git_hash])
    String effective_workspace_bucket = select_first([workspace_bucket, GetToolVersions.workspace_bucket])
    String effective_workspace_id = select_first([workspace_id, GetToolVersions.workspace_id])
    String effective_submission_id = select_first([submission_id, GetToolVersions.submission_id])

    call Utils.GetReference {
        input:
            reference_name = reference_name,
            basic_docker = effective_basic_docker,
    }

    # If `is_wgs` is true we'll use the WGS interval list else, otherwise we'll use the Exome interval list.
    # However if `interval_list` is defined, we'll use that instead of choosing based on `is_wgs`.
    File default_interval_list = if (is_wgs) then GetReference.reference.wgs_calling_interval_list
                                 else GetReference.reference.exome_calling_interval_list
    File interval_list_to_use = select_first([interval_list, default_interval_list])

    call PrepareReferenceFiles.GenerateBgzSequenceDictionaryAndIndex {
        input:
            reference_fasta = select_first([reference]),
            output_gcs_dir = effective_workspace_bucket + "/submissions/" + effective_submission_id,
            gatk_docker = effective_gatk_docker,
    }

    call PrepareReferenceFiles.GenerateContigMapping {
        input:
            sequence_dictionary = read_json(GenerateBgzSequenceDictionaryAndIndex.reference_files_json).sequence_dictionary,
            in_reference_json = GenerateBgzSequenceDictionaryAndIndex.reference_files_json,
            output_gcs_dir = effective_workspace_bucket + "/submissions/" + effective_submission_id,
            variants_docker = effective_variants_docker,
    }

    call BulkIngestGenomes.GvsBulkIngestGenomes as BulkIngestGenomes {
        input:
            git_branch_or_tag = git_branch_or_tag,
            git_hash = effective_git_hash,
            dataset_name = dataset_name,
            project_id = project_id,
            basic_docker = effective_basic_docker,
            cloud_sdk_docker = effective_cloud_sdk_docker,
            variants_docker = effective_variants_docker,
            gatk_docker = effective_gatk_docker,
            gatk_override = gatk_override,
            reference_name = reference_name,
            interval_list = interval_list_to_use,
            custom_ref_dictionary = read_json(GenerateBgzSequenceDictionaryAndIndex.reference_files_json).sequence_dictionary,
            custom_contig_mapping = read_json(GenerateContigMapping.reference_files_json).contig_mapping,
            drop_state = drop_state,
            sample_id_column_name = sample_id_column_name,
            vcf_files_column_name = vcf_files_column_name,
            vcf_index_files_column_name = vcf_index_files_column_name,
            sample_set_name = sample_set_name,
            billing_project_id = billing_project_id,
            use_compressed_references = use_compressed_references,
            load_vcf_headers = load_vcf_headers,
            load_vet_and_ref_ranges = load_vet_and_ref_ranges,
            workspace_bucket = effective_workspace_bucket,
            workspace_id = effective_workspace_id,
            tighter_gcp_quotas = tighter_gcp_quotas,
            is_wgs = is_wgs,
    }

    call PrepareReferenceFiles.CreateWeightedBedFile {
        input:
            project_id = project_id,
            dataset_name = dataset_name,
            reference_dictionary = read_json(GenerateBgzSequenceDictionaryAndIndex.reference_files_json).sequence_dictionary,
            contig_mapping = read_json(GenerateContigMapping.reference_files_json).contig_mapping,
            in_reference_json = GenerateContigMapping.reference_files_json,
            output_gcs_dir = effective_workspace_bucket + "/submissions/" + effective_submission_id,
            variants_docker = effective_variants_docker,
    }

    call PopulateAltAllele.GvsPopulateAltAllele {
        input:
            git_branch_or_tag = git_branch_or_tag,
            git_hash = effective_git_hash,
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
            git_hash = effective_git_hash,
            dataset_name = dataset_name,
            project_id = project_id,
            call_set_identifier = call_set_identifier,
            filter_set_name = effective_filter_set_name,
            use_VETS = !use_VQSR,
            reference_name = reference_name,
            interval_list = interval_list_to_use,
            custom_reference = reference,
            custom_contig_mapping = read_json(GenerateContigMapping.reference_files_json).contig_mapping,
            custom_training_resources = "gs://fc-59f20e4b-b37a-4506-8b09-1eae66d3b0e4/fake_sites_only.vcf",
            variants_docker = effective_variants_docker,
            gatk_docker = effective_gatk_docker,
            gatk_override = gatk_override,
            INDEL_VQSR_max_gaussians_override = INDEL_VQSR_max_gaussians_override,
            INDEL_VQSR_mem_gb_override = INDEL_VQSR_mem_gb_override,
            SNP_VQSR_max_gaussians_override = SNP_VQSR_max_gaussians_override,
            SNP_VQSR_mem_gb_override = SNP_VQSR_mem_gb_override,
            cloud_sdk_docker = effective_cloud_sdk_docker,
            training_python_script = training_python_script,
            scoring_python_script = scoring_python_script,
    }

    call PrepareRangesCallset.GvsPrepareCallset {
        input:
            call_set_identifier = call_set_identifier,
            go = GvsCreateFilterSet.done,
            git_branch_or_tag = git_branch_or_tag,
            git_hash = effective_git_hash,
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
            cloud_sdk_docker = effective_cloud_sdk_docker,
            enable_extract_table_ttl = true,
    }

    call ExtractCallset.GvsExtractCallset {
        input:
            go = GvsPrepareCallset.done,
            git_branch_or_tag = git_branch_or_tag,
            git_hash = effective_git_hash,
            dataset_name = dataset_name,
            project_id = project_id,
            call_set_identifier = call_set_identifier,
            extract_table_prefix = effective_extract_table_prefix,
            filter_set_name = effective_filter_set_name,
            query_project = query_project,
            scatter_count = extract_scatter_count,
            reference_name = reference_name,
            interval_list = interval_list_to_use,
            interval_weights_bed = read_json(CreateWeightedBedFile.updated_reference_files_json).weighted_bed,
            custom_reference = reference,
            custom_contig_mapping = read_json(GenerateContigMapping.reference_files_json).contig_mapping,
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
            bgzip_output_vcfs = bgzip_output_vcfs,
            collect_variant_calling_metrics = collect_variant_calling_metrics,
            is_wgs = is_wgs,
            maximum_alternate_alleles = maximum_alternate_alleles,
            target_interval_list = target_interval_list,
            ploidy_table_name = ploidy_table_name,
    }

    output {
        String output_gcs_path = extract_output_gcs_dir
        Array[File] output_vcfs = GvsExtractCallset.output_vcfs
        Array[File] output_vcf_indexes = GvsExtractCallset.output_vcf_indexes
        Array[File] output_vcf_interval_files = GvsExtractCallset.output_vcf_interval_files
        Float total_vcfs_size_mb = GvsExtractCallset.total_vcfs_size_mb
        File? sample_name_list = GvsExtractCallset.sample_name_list
        File manifest = GvsExtractCallset.manifest
        String recorded_git_hash = effective_git_hash
        Boolean done = true
        Boolean used_tighter_gcp_quotas = BulkIngestGenomes.used_tighter_gcp_quotas
    }
}
