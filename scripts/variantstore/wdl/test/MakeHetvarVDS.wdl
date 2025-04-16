version 1.0

import "../GvsUtils.wdl" as Utils
import "../GvsExtractAvroFilesForHail.wdl" as ExtractAvroFilesForHail
import "../GvsCreateVDS.wdl" as CreateVds
import "GvsQuickstartVcfIntegration.wdl" as QuickstartVcfIntegration
import "GvsTieOutVDS.wdl" as TieOutVDS

workflow MakeHetvarVDS {
    input {
        String git_branch_or_tag
        String? git_hash
        String reference_name = "hg38"
        Boolean is_wgs = true
        String dataset_name
        String filter_set_name
        File? interval_list
        Boolean use_compressed_references = false
        Boolean extract_do_not_filter_override = false
        String dataset_suffix = "hail"
        Boolean use_default_dockers = false
        Boolean bgzip_output_vcfs = false

        String? basic_docker
        String? cloud_sdk_docker
        String? cloud_sdk_slim_docker
        String? variants_docker
        String? gatk_docker

        File? gatk_override
        String? hail_version
        String? sample_set_name ## NOTE: currently we only allow the loading of one sample set at a time

        String? workspace_bucket
        String? workspace_id
        String? submission_id

        Int? maximum_alternate_alleles
        String ploidy_table_name = "sample_chromosome_ploidy"
    }

    String project_id = "gvs-internal"

    if (!defined(workspace_bucket) || !defined(workspace_id) || !defined(submission_id) ||
        !defined(git_hash) || !defined(basic_docker) || !defined(cloud_sdk_docker) || !defined(cloud_sdk_slim_docker) ||
        !defined(variants_docker) || !defined(gatk_docker) || !defined(hail_version)) {
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
    String effective_hail_version = select_first([hail_version, GetToolVersions.hail_version])

    String effective_workspace_bucket = select_first([workspace_bucket, GetToolVersions.workspace_bucket])
    String effective_workspace_id = select_first([workspace_id, GetToolVersions.workspace_id])
    String effective_submission_id = select_first([submission_id, GetToolVersions.submission_id])


    call ExtractAvroFilesForHail.GvsExtractAvroFilesForHail {
        input:
            git_branch_or_tag = git_branch_or_tag,
            git_hash = git_hash,
            project_id = project_id,
            dataset_name = dataset_name,
            filter_set_name = filter_set_name,
            ploidy_table_name = ploidy_table_name,
            scatter_width = 10,
            call_set_identifier = git_branch_or_tag,
            basic_docker = effective_basic_docker,
            cloud_sdk_docker = effective_cloud_sdk_docker,
            variants_docker = effective_variants_docker,
    }

    String vds_path = GvsExtractAvroFilesForHail.avro_path + "/gvs_export.vds"

    call CreateVds.GvsCreateVDS {
        input:
            git_branch_or_tag = git_branch_or_tag,
            hail_version = effective_hail_version,
            avro_path = GvsExtractAvroFilesForHail.avro_path,
            vds_destination_path = vds_path,
            cluster_prefix = "vds-cluster",
            gcs_subnetwork_name = "subnetwork",
            region = "us-central1",
            basic_docker = effective_basic_docker,
            variants_docker = effective_variants_docker,
            cloud_sdk_slim_docker = effective_cloud_sdk_slim_docker,
            cluster_max_age_minutes = 120,
            cluster_max_idle_minutes = 60,
            leave_cluster_running_at_end = false,
    }

    output {
        String vds_output_path = GvsExtractAvroFilesForHail.avro_path + "/gvs_export.vds"
        String recorded_git_hash = effective_git_hash
        Boolean done = true
    }
}
