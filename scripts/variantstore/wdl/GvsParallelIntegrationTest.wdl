version 1.0

import "GvsUtils.wdl" as Utils
import "GvsJointVariantCalling.wdl" as JointVariantCalling
import "GvsBulkIngestGenomes.wdl" as BulkIngestGenomes
import "GvsPopulateAltAllele.wdl" as PopulateAltAllele
import "GvsCreateFilterSet.wdl" as CreateFilterSet
import "GvsPrepareRangesCallset.wdl" as PrepareRangesCallset

workflow GvsQuickstartVcfIntegration {
    input {
        String git_branch_or_tag
        String? git_hash
        String expected_output_prefix
        Boolean use_VQSR_lite = true
        Boolean extract_do_not_filter_override = true
        Boolean use_compressed_references = true
        Boolean process_vcf_headers = false
        String drop_state = "FORTY"
        String dataset_suffix
        String snapshot_run_name
        Boolean is_wgs = true
        File? interval_list
        Boolean use_default_dockers = false
        String? basic_docker
        String? cloud_sdk_docker
        String? cloud_sdk_slim_docker
        String? variants_docker
        String? gatk_docker
        File? gatk_override
        String? sample_id_column_name ## Note that a column WILL exist that is the <entity>_id from the table name. However, some users will want to specify an alternate column for the sample_name during ingest
        String? vcf_files_column_name
        String? vcf_index_files_column_name
        String? sample_set_name ## NOTE: currently we only allow the loading of one sample set at a time

        String? workspace_bucket
        String? workspace_id
        String? submission_id
    }
    String project_id = "gvs-internal"
    String snapshot_dataset = "hatcher_vs_299_test_storage"

    # WDL 1.0 trick to set a variable ('none') to be undefined.
    if (false) {
        File? none = ""
    }

    if (!defined(workspace_bucket) || !defined(workspace_id) || !defined(submission_id) ||
        !defined(git_hash) || !defined(cloud_sdk_docker) || !defined(cloud_sdk_slim_docker) ||
        !defined(variants_docker) || !defined(basic_docker) || !defined(gatk_docker)) {
        call Utils.GetToolVersions {
            input:
                git_branch_or_tag = git_branch_or_tag,
        }
    }

    String effective_basic_docker = select_first([basic_docker, GetToolVersions.basic_docker])
    String effective_cloud_sdk_docker = select_first([cloud_sdk_docker, GetToolVersions.cloud_sdk_docker])
    String effective_cloud_sdk_slim_docker = select_first([cloud_sdk_slim_docker, GetToolVersions.cloud_sdk_slim_docker])
    String effective_variants_docker = select_first([variants_docker, GetToolVersions.variants_docker])
    String effective_gatk_docker = select_first([gatk_docker, GetToolVersions.gatk_docker])
    String effective_git_hash = select_first([git_hash, GetToolVersions.git_hash])

    String effective_workspace_bucket = select_first([workspace_bucket, GetToolVersions.workspace_bucket])
    String effective_workspace_id = select_first([workspace_id, GetToolVersions.workspace_id])
    String effective_submission_id = select_first([submission_id, GetToolVersions.submission_id])

    if (!use_default_dockers && !defined(gatk_override)) {
        call Utils.BuildGATKJar {
            input:
                git_branch_or_tag = git_branch_or_tag,
                cloud_sdk_slim_docker = effective_cloud_sdk_slim_docker,
        }
    }


    # Setting up the Bulkingest portion of the test.
    # Just make a new dataset, run that stage, then compare output.  It's the first stage, so no setup required
    call Utils.CreateDatasetForTest as CreateDatasetForIngestRun {
        input:
            git_branch_or_tag = git_branch_or_tag,
            dataset_prefix = "quickit_parallel_ingest",
            dataset_suffix = dataset_suffix,
            cloud_sdk_docker = effective_cloud_sdk_docker,
    }

    call BulkIngestGenomes.GvsBulkIngestGenomes as BulkIngestGenomes {
        input:
            git_branch_or_tag = git_branch_or_tag,
            git_hash = effective_git_hash,
            dataset_name = CreateDatasetForIngestRun.dataset_name,
            project_id = project_id,
            basic_docker = effective_basic_docker,
            cloud_sdk_docker = effective_cloud_sdk_docker,
            variants_docker = effective_variants_docker,
            gatk_docker = effective_gatk_docker,
            gatk_override = gatk_override,
            drop_state = drop_state,
            sample_id_column_name = sample_id_column_name,
            vcf_files_column_name = vcf_files_column_name,
            vcf_index_files_column_name = vcf_index_files_column_name,
            sample_set_name = sample_set_name,
            use_compressed_references = use_compressed_references,
            process_vcf_headers = process_vcf_headers,
            workspace_bucket = effective_workspace_bucket,
            workspace_id = effective_workspace_id,
            tighter_gcp_quotas = false,
            is_wgs = is_wgs,
    }

    call Utils.CheckResultsAgainstStoredState as CheckResultsOfBulkIngestRun {
        input:
            go = BulkIngestGenomes.done,
            project_id = project_id,
            run_dataset = CreateDatasetForIngestRun.dataset_name,
            snapshot_dataset = snapshot_dataset,
            run_name = snapshot_run_name,
            comparison_key = "after_bulk_ingest_before_alt_allele",
            cloud_sdk_docker = effective_cloud_sdk_docker,
    }
    # DONE WITH TESTING BULK INGEST



    # Setting up the alt allele portion of the test.
    # Make a new dataset, reestablish the tables, run that stage, then compare output.
    call Utils.CreateDatasetForTest as CreateDatasetForAltAlleleRun {
        input:
            git_branch_or_tag = git_branch_or_tag,
            dataset_prefix = "quickit_parallel_alt_allele",
            dataset_suffix = dataset_suffix,
            cloud_sdk_docker = effective_cloud_sdk_docker,
    }

    call Utils.RestoreSnapshotForRun as RestoreSnapshotForAltAllele {
        input:
            project_id = project_id,
            dest_dataset = CreateDatasetForAltAlleleRun.dataset_name,
            snapshot_dataset = snapshot_dataset,
            run_name = snapshot_run_name,
            start_retrieval_key = "after_bulk_ingest_before_alt_allele",
            end_retrieval_key = "after_alt_allele_before_filter_set",
            cloud_sdk_docker = effective_cloud_sdk_docker,
    }

    call PopulateAltAllele.GvsPopulateAltAllele {
        input:
            git_branch_or_tag = git_branch_or_tag,
            git_hash = effective_git_hash,
            call_set_identifier = git_branch_or_tag,
            go = RestoreSnapshotForAltAllele.done,
            dataset_name = CreateDatasetForAltAlleleRun.dataset_name,
            project_id = project_id,
            variants_docker = effective_variants_docker,
            cloud_sdk_docker = effective_cloud_sdk_docker,
    }

    call Utils.CheckResultsAgainstStoredState as CheckResultsOfAltAllele {
        input:
            go = GvsPopulateAltAllele.done,
            project_id = project_id,
            run_dataset = CreateDatasetForAltAlleleRun.dataset_name,
            snapshot_dataset = snapshot_dataset,
            run_name = snapshot_run_name,
            comparison_key = "after_alt_allele_before_filter_set",
            cloud_sdk_docker = effective_cloud_sdk_docker,
    }
    # DONE WITH TESTING Alt Allele


    # Setting up the alt allele portion of the test.
    # Make a new dataset, reestablish the tables, run that stage, then compare output.
    call Utils.CreateDatasetForTest as CreateDatasetForCreateFilterRun {
        input:
            git_branch_or_tag = git_branch_or_tag,
            dataset_prefix = "quickit_parallel_filter_creation",
            dataset_suffix = dataset_suffix,
            cloud_sdk_docker = effective_cloud_sdk_docker,
    }

    call Utils.RestoreSnapshotForRun as RestoreSnapshotForCreateFilter {
        input:
            project_id = project_id,
            dest_dataset = CreateDatasetForCreateFilterRun.dataset_name,
            snapshot_dataset = snapshot_dataset,
            run_name = snapshot_run_name,
            start_retrieval_key = "after_alt_allele_before_filter_set",
            end_retrieval_key = "after_filter_set_before_prepare",
            cloud_sdk_docker = effective_cloud_sdk_docker,
    }

    call CreateFilterSet.GvsCreateFilterSet {
        input:
            go = RestoreSnapshotForCreateFilter.done,
            git_branch_or_tag = git_branch_or_tag,
            git_hash = effective_git_hash,
            dataset_name = CreateDatasetForCreateFilterRun.dataset_name,
            project_id = project_id,
            call_set_identifier = git_branch_or_tag,
            filter_set_name = "quickit_parallel",
            use_VQSR_lite = use_VQSR_lite,
            variants_docker = effective_variants_docker,
            gatk_docker = effective_gatk_docker,
            gatk_override = gatk_override,
            cloud_sdk_docker = effective_cloud_sdk_docker,
    }

    call Utils.CheckResultsAgainstStoredState as CheckResultsOfCreateFilterSet {
        input:
            go = GvsCreateFilterSet.done,
            project_id = project_id,
            run_dataset = CreateDatasetForCreateFilterRun.dataset_name,
            snapshot_dataset = snapshot_dataset,
            run_name = snapshot_run_name,
            comparison_key = "after_filter_set_before_prepare",
            cloud_sdk_docker = effective_cloud_sdk_docker,
    }
    # DONE WITH TESTING Filter Set Creation


    output {
        String filter_set_name = "quickit"
        String recorded_git_hash = effective_git_hash
        Boolean done = true
        Boolean used_tighter_gcp_quotas = false
    }
}