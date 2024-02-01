version 1.0

import "GvsQuickstartVcfIntegration.wdl" as QuickstartVcfIntegration
import "GvsQuickstartHailIntegration.wdl" as QuickstartHailIntegration
import "GvsJointVariantCalling.wdl" as JointVariantCalling
import "GvsUtils.wdl" as Utils

workflow GvsQuickstartIntegration {
    input {
        String git_branch_or_tag
        Boolean use_default_dockers = false
        Boolean run_vcf_integration = true
        Boolean run_hail_integration = true
        Boolean run_exome_integration = true
        Boolean run_beta_integration = true
        String? wgs_sample_id_column_name ## Note that a column WILL exist that is the <entity>_id from the table name. However, some users will want to specify an alternate column for the sample_name during ingest
        String? wgs_vcf_files_column_name
        String? wgs_vcf_index_files_column_name
        String? wgs_sample_set_name ## NOTE: currently we only allow the loading of one sample set at a time
        String? exome_sample_id_column_name ## Note that a column WILL exist that is the <entity>_id from the table name. However, some users will want to specify an alternate column for the sample_name during ingest
        String? exome_vcf_files_column_name
        String? exome_vcf_index_files_column_name
        String? exome_sample_set_name ## NOTE: currently we only allow the loading of one sample set at a time
        String? basic_docker
        String? cloud_sdk_docker
        String? cloud_sdk_slim_docker
        String? variants_docker
        String? gatk_docker
        String? hail_version
        Boolean chr20_X_Y_only = true
    }

    File full_wgs_interval_list = "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list"
    File full_exome_interval_list = "gs://gcp-public-data--broad-references/hg38/v0/bge_exome_calling_regions.v1.1.interval_list"
    File expected_output_prefix = "gs://gvs-internal-quickstart/integration/2024-01-18/"

    Boolean process_vcf_headers = false
    Boolean load_headers_only = true

    # WDL 1.0 trick to set a variable ('none') to be undefined.
    if (false) {
        File? none = ""
    }

    # Always call `GetToolVersions` to get the git hash for this run as this is a top-level-only WDL (i.e. there are
    # no calling WDLs that might supply `git_hash`).
    call Utils.GetToolVersions {
        input:
            git_branch_or_tag = git_branch_or_tag,
    }

    String effective_basic_docker = select_first([basic_docker, GetToolVersions.basic_docker])
    String effective_cloud_sdk_docker = select_first([cloud_sdk_docker, GetToolVersions.cloud_sdk_docker])
    String effective_cloud_sdk_slim_docker = select_first([cloud_sdk_slim_docker, GetToolVersions.cloud_sdk_slim_docker])
    String effective_variants_docker = select_first([variants_docker, GetToolVersions.variants_docker])
    String effective_gatk_docker = select_first([gatk_docker, GetToolVersions.gatk_docker])
    String effective_hail_version = select_first([hail_version, GetToolVersions.hail_version])

    if (chr20_X_Y_only) {
        call FilterIntervalListChromosomes {
            input:
                full_interval_list = full_wgs_interval_list,
                chromosomes = ["chrX", "chrY", "chr20"],
                variants_docker = effective_variants_docker,
        }
    }

    if (!use_default_dockers) {
        call Utils.BuildGATKJar {
            input:
                git_branch_or_tag = git_branch_or_tag,
                cloud_sdk_slim_docker = effective_cloud_sdk_slim_docker,
        }
    }

    # Note for `GvsQuickstartIntegration` we use the git_branch_or_tag *input* and its corresponding git hash. This is not
    # necessarily the same as the branch name selected in Terra for the integration `GvsQuickstartIntegration` workflow,
    # though in practice likely they are the same.
    if (run_hail_integration) {
        call QuickstartHailIntegration.GvsQuickstartHailIntegration as GvsQuickstartHailVQSRLiteIntegration {
            input:
                git_branch_or_tag = git_branch_or_tag,
                git_hash = GetToolVersions.git_hash,
                use_VQSR_lite = true,
                extract_do_not_filter_override = false,
                dataset_suffix = "lite_hail",
                use_default_dockers = use_default_dockers,
                gatk_override = if (use_default_dockers) then none else BuildGATKJar.jar,
                is_wgs = true,
                interval_list = select_first([FilterIntervalListChromosomes.out, full_wgs_interval_list]),
                expected_output_prefix = expected_output_prefix,
                sample_id_column_name = wgs_sample_id_column_name,
                vcf_files_column_name = wgs_vcf_files_column_name,
                vcf_index_files_column_name = wgs_vcf_index_files_column_name,
                sample_set_name = select_first([wgs_sample_set_name, "wgs_integration_sample_set"]),
                basic_docker = effective_basic_docker,
                cloud_sdk_docker = effective_cloud_sdk_docker,
                cloud_sdk_slim_docker = effective_cloud_sdk_slim_docker,
                variants_docker = effective_variants_docker,
                gatk_docker = effective_gatk_docker,
                workspace_bucket = GetToolVersions.workspace_bucket,
                workspace_id = GetToolVersions.workspace_id,
                submission_id = GetToolVersions.submission_id,
                hail_version = effective_hail_version,
                process_vcf_headers = process_vcf_headers,
                load_headers_only = load_headers_only,
        }
        call QuickstartHailIntegration.GvsQuickstartHailIntegration as GvsQuickstartHailVQSRClassicIntegration {
            input:
                git_branch_or_tag = git_branch_or_tag,
                git_hash = GetToolVersions.git_hash,
                use_VQSR_lite = false,
                extract_do_not_filter_override = false,
                dataset_suffix = "classic_hail",
                use_default_dockers = use_default_dockers,
                gatk_override = if (use_default_dockers) then none else BuildGATKJar.jar,
                is_wgs = true,
                interval_list = FilterIntervalListChromosomes.out,
                expected_output_prefix = expected_output_prefix,
                sample_id_column_name = wgs_sample_id_column_name,
                vcf_files_column_name = wgs_vcf_files_column_name,
                vcf_index_files_column_name = wgs_vcf_index_files_column_name,
                sample_set_name = select_first([wgs_sample_set_name, "wgs_integration_sample_set"]),
                basic_docker = effective_basic_docker,
                cloud_sdk_docker = effective_cloud_sdk_docker,
                cloud_sdk_slim_docker = effective_cloud_sdk_slim_docker,
                variants_docker = effective_variants_docker,
                gatk_docker = effective_gatk_docker,
                workspace_bucket = GetToolVersions.workspace_bucket,
                workspace_id = GetToolVersions.workspace_id,
                submission_id = GetToolVersions.submission_id,
                hail_version = effective_hail_version,
                process_vcf_headers = process_vcf_headers,
                load_headers_only = load_headers_only,
        }

        if (GvsQuickstartHailVQSRLiteIntegration.used_tighter_gcp_quotas) {
            call Utils.TerminateWorkflow as HailVQSRLiteQuotaFail {
                input:
                    message = "GvsQuickstartHailVQSRLiteIntegration should not have used tighter GCP quotas but did!",
                    basic_docker = effective_basic_docker,
            }
        }

        if (GvsQuickstartHailVQSRClassicIntegration.used_tighter_gcp_quotas) {
            call Utils.TerminateWorkflow as HailVQSRClassicQuotaFail {
                input:
                    message = "GvsQuickstartHailVQSRClassicIntegration should not have used tighter GCP quotas but did!",
                    basic_docker = effective_basic_docker,
            }
        }
    }

    if (run_vcf_integration) {
        call QuickstartVcfIntegration.GvsQuickstartVcfIntegration as QuickstartVcfVQSRLiteIntegration {
            input:
                git_branch_or_tag = git_branch_or_tag,
                git_hash = GetToolVersions.git_hash,
                use_VQSR_lite = true,
                extract_do_not_filter_override = false,
                dataset_suffix = "lite_vcf",
                use_default_dockers = use_default_dockers,
                gatk_override = if (use_default_dockers) then none else BuildGATKJar.jar,
                is_wgs = true,
                interval_list = FilterIntervalListChromosomes.out,
                expected_output_prefix = expected_output_prefix,
                sample_id_column_name = wgs_sample_id_column_name,
                vcf_files_column_name = wgs_vcf_files_column_name,
                vcf_index_files_column_name = wgs_vcf_index_files_column_name,
                sample_set_name = select_first([wgs_sample_set_name, "wgs_integration_sample_set"]),
                drop_state = "FORTY",
                basic_docker = effective_basic_docker,
                cloud_sdk_docker = effective_cloud_sdk_docker,
                cloud_sdk_slim_docker = effective_cloud_sdk_slim_docker,
                variants_docker = effective_variants_docker,
                gatk_docker = effective_gatk_docker,
                workspace_bucket = GetToolVersions.workspace_bucket,
                workspace_id = GetToolVersions.workspace_id,
                submission_id = GetToolVersions.submission_id,
                process_vcf_headers = process_vcf_headers,
                load_headers_only = load_headers_only,
        }
        call QuickstartVcfIntegration.GvsQuickstartVcfIntegration as QuickstartVcfVQSRClassicIntegration {
            input:
                git_branch_or_tag = git_branch_or_tag,
                git_hash = GetToolVersions.git_hash,
                use_VQSR_lite = false,
                extract_do_not_filter_override = true,
                dataset_suffix = "classic_vcf",
                use_default_dockers = use_default_dockers,
                gatk_override = if (use_default_dockers) then none else BuildGATKJar.jar,
                is_wgs = true,
                interval_list = FilterIntervalListChromosomes.out,
                expected_output_prefix = expected_output_prefix,
                sample_id_column_name = wgs_sample_id_column_name,
                vcf_files_column_name = wgs_vcf_files_column_name,
                vcf_index_files_column_name = wgs_vcf_index_files_column_name,
                sample_set_name = select_first([wgs_sample_set_name, "wgs_integration_sample_set"]),
                drop_state = "FORTY",
                basic_docker = effective_basic_docker,
                cloud_sdk_docker = effective_cloud_sdk_docker,
                cloud_sdk_slim_docker = effective_cloud_sdk_slim_docker,
                variants_docker = effective_variants_docker,
                gatk_docker = effective_gatk_docker,
                workspace_bucket = GetToolVersions.workspace_bucket,
                workspace_id = GetToolVersions.workspace_id,
                submission_id = GetToolVersions.submission_id,
                process_vcf_headers = process_vcf_headers,
                load_headers_only = load_headers_only,
        }

        if (QuickstartVcfVQSRClassicIntegration.used_tighter_gcp_quotas) {
            call Utils.TerminateWorkflow as VcfVQSRClassicQuotaFail {
                input:
                    message = "QuickstartVcfVQSRLiteIntegration should not have used tighter GCP quotas but did!",
                    basic_docker = effective_basic_docker,
            }
        }

        if (QuickstartVcfVQSRLiteIntegration.used_tighter_gcp_quotas) {
            call Utils.TerminateWorkflow as VcfVQSRLiteQuotaFail {
                input:
                    message = "QuickstartVcfVQSRLiteIntegration should not have used tighter GCP quotas but did!",
                    basic_docker = effective_basic_docker,
            }
        }
    }

    if (run_exome_integration) {
        call QuickstartVcfIntegration.GvsQuickstartVcfIntegration as QuickstartVcfExomeIntegration {
            input:
                git_branch_or_tag = git_branch_or_tag,
                git_hash = GetToolVersions.git_hash,
                use_VQSR_lite = true,
                extract_do_not_filter_override = false,
                dataset_suffix = "exome",
                use_default_dockers = use_default_dockers,
                gatk_override = if (use_default_dockers) then none else BuildGATKJar.jar,
                is_wgs = false,
                interval_list = full_exome_interval_list,
                expected_output_prefix = expected_output_prefix,
                sample_id_column_name = exome_sample_id_column_name,
                vcf_files_column_name = exome_vcf_files_column_name,
                vcf_index_files_column_name = exome_vcf_index_files_column_name,
                sample_set_name = select_first([exome_sample_set_name, "exome_integration_sample_set"]),
                drop_state = "FORTY",
                basic_docker = effective_basic_docker,
                cloud_sdk_docker = effective_cloud_sdk_docker,
                cloud_sdk_slim_docker = effective_cloud_sdk_slim_docker,
                variants_docker = effective_variants_docker,
                gatk_docker = effective_gatk_docker,
                workspace_bucket = GetToolVersions.workspace_bucket,
                workspace_id = GetToolVersions.workspace_id,
                submission_id = GetToolVersions.submission_id,
        }

        if (QuickstartVcfExomeIntegration.used_tighter_gcp_quotas) {
            call Utils.TerminateWorkflow as ExomeQuotaFail {
                input:
                    message = "QuickstartVcfExomeIntegration should not have used tighter GCP quotas but did!",
                    basic_docker = effective_basic_docker,
            }
        }
    }

    if (run_beta_integration) {
        String project_id = "gvs-internal"

        call Utils.CreateDatasetForTest {
            input:
                git_branch_or_tag = git_branch_or_tag,
                dataset_prefix = "quickit",
                dataset_suffix = "beta",
                cloud_sdk_docker = effective_cloud_sdk_docker,
        }

        call JointVariantCalling.GvsJointVariantCalling as QuickstartBeta {
            input:
                call_set_identifier = git_branch_or_tag,
                dataset_name = CreateDatasetForTest.dataset_name,
                project_id = project_id,
                gatk_override = if (use_default_dockers) then none else BuildGATKJar.jar,
                extract_output_file_base_name = "quickit",
                filter_set_name = "quickit",
                extract_table_prefix = "quickit",
                basic_docker = effective_basic_docker,
                cloud_sdk_docker = effective_cloud_sdk_docker,
                variants_docker = effective_variants_docker,
                gatk_docker = effective_gatk_docker,
                workspace_bucket = GetToolVersions.workspace_bucket,
                workspace_id = GetToolVersions.workspace_id,
                submission_id = GetToolVersions.submission_id,
                git_branch_or_tag = git_branch_or_tag,
        }

        if (!QuickstartBeta.used_tighter_gcp_quotas) {
            call Utils.TerminateWorkflow as QuickstartBetaQuotaFail {
                input:
                    message = "QuickstartBeta should have used tighter GCP quotas but did not!",
                    basic_docker = effective_basic_docker,
            }
        }
    }

    output {
        String recorded_git_hash = GetToolVersions.git_hash
    }
}

task FilterIntervalListChromosomes {
    input {
        File full_interval_list
        Array[String]+ chromosomes
        String variants_docker
    }
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        python3 /app/filter_interval_list_chromosomes.py --input-interval-list ~{full_interval_list} \
          --output-interval-list "filtered.interval_list" --chromosome ~{sep=' --chromosome ' chromosomes}
    >>>
    runtime {
        docker: variants_docker
    }
    output {
        File out = "filtered.interval_list"
    }
}


