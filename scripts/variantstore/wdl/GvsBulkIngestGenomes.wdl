version 1.0

import "GvsUtils.wdl" as Utils
import "GvsAssignIds.wdl" as AssignIds
import "GvsImportGenomes.wdl" as ImportGenomes
import "GvsValidateVcfHeaders.wdl" as ValidateVcfHeaders

workflow GvsBulkIngestGenomes {
    input {
        # Begin GenerateImportFofnFromDataTable
        # for now set the entity type names with a default
        String data_table_name = "sample"
        # Note that a column WILL exist that is the <entity>_id from the table name. However, some users will want to
        # specify an alternate column for the sample_name during ingest
        String? sample_id_column_name
        String? vcf_files_column_name
        String? vcf_index_files_column_name
        String? sample_set_name
        # Optional FOFN of VCFs to ingest. If specified, the workflow will use this and not generate a FOFN from the data table.
        File? bulk_ingest_fofn
        # End GenerateImportFofnFromDataTable

        # Begin GvsAssignIds
        String dataset_name
        String project_id
        Boolean samples_are_controls = false
        String? basic_docker
        String? cloud_sdk_docker
        String? variants_docker
        String? gatk_docker
        String? git_branch_or_tag
        String? git_hash

        String? workspace_id
        String? workspace_bucket

        File? gatk_override
        # End GvsAssignIds

        # Begin GvsImportGenomes
        String reference_name = "hg38"
        File? interval_list

        # set to "NONE" to ingest all the reference data into GVS for VDS (instead of VCF) output
        String drop_state = "NONE"

        Int? load_data_scatter_width
        Int? load_data_preemptible_override
        Int? load_data_maxretries_override
        String? billing_project_id
        Boolean load_vet_and_ref_ranges = true
        # Controls whether VCF headers are loaded during data ingest when validate_vcf_headers is false.
        # Note: when validate_vcf_headers is true (default), VCF headers are automatically loaded during
        # the initial validation pass regardless of this setting. Set this to true only when bypassing
        # validation (validate_vcf_headers=false) while still wanting headers loaded in a single direct pass.
        Boolean load_vcf_headers = false
        Boolean tighter_gcp_quotas = false
        Boolean is_wgs = true
        # End GvsImportGenomes

        # Begin GvsValidateVcfHeaders (VS-1966 / VS-1995)
        # When true (default), validate the VCF headers before vet/ref data ingest. An initial headers-only
        # ingest pass populates header tables and runs GvsValidateVcfHeaders. If validation passes, the
        # workflow proceeds to vet/ref data ingest; if validation fails and fail_on_validation_errors is true,
        # it halts fast before any expensive vet/ref compute is spent. Set to false to bypass pre-ingest
        # header checks.
        Boolean validate_vcf_headers = true
        # Exact triplet ('3.7.8', AoU) or a range with optional interval notation ('3.4.12-3.7.8',
        # '[3.7.8-3.8)', '(3.7-3.8)'); see GvsValidateVcfHeaders.
        String? expected_dragen_version
        Boolean require_reblocking = true
        Boolean fail_on_validation_errors = true
        # End GvsValidateVcfHeaders

        Boolean use_parquet_ingest = true
        # `parquet_output_gcs_dir` must be defined if `use_parquet_ingest` is true.
        String? parquet_output_gcs_dir

        Boolean use_compressed_references = false
    }

    parameter_meta {
        data_table_name: "The name of the data table; This table holds the GVCFs to be ingested; `sample` is the default."
        sample_id_column_name: "The column that will be used for the sample name / id in GVS; the <data_table_name>_id will be used as the default"
        vcf_files_column_name: "The column that supplies the path for the GVCFs to be ingested. If not specified, the workflow will attempt to derive the column name."
        vcf_index_files_column_name: "The column that supplies the path for the GVCF index files to be ingested. If not specified, the workflow will attempt to derive the column name."
        sample_set_name: "The recommended way to load samples; Sample sets must be created by the user. If no sample_set_name is specified, all samples will be loaded into GVS"
        bulk_ingest_fofn: "An explicitly specified FOFN of VCFs to be ingested. If specified, the workflow will not generate a FOFN from the data table. This can be useful for avoiding the scale limitations of Terra data tables. The format is tab delimited with no header: sample_name<tab>gvcf_file_path<tab>gvcf_index_file_path. If this value is specified, none of the data table parameters should be specified."
    }

    if (!defined(git_hash) ||
        !defined(basic_docker) || !defined(cloud_sdk_docker) || !defined(variants_docker) || !defined(gatk_docker) ||
        !defined(workspace_id) || !defined(workspace_bucket)) {
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
    String effective_workspace_id = select_first([workspace_id, GetToolVersions.workspace_id])
    String effective_workspace_bucket = select_first([workspace_bucket, GetToolVersions.workspace_bucket])

    # WDL 1.0 idiom for an undefined optional string
    if (false) {
        String? none_string = ""
    }

    if (defined(bulk_ingest_fofn) && (defined(sample_id_column_name) || defined(vcf_files_column_name) || defined(vcf_index_files_column_name) || defined(sample_set_name))) {
        call Utils.TerminateWorkflow as MustNotSpecifyDataTableParamsWhenBulkIngestFofnSpecified {
            input:
                message = "GvsBulkIngestGenomes called with bulk_ingest_fofn specified, but also with at least one data table parameter specified (sample_id_column_name, vcf_files_column_name, vcf_index_files_column_name, sample_set_name).",
                basic_docker = effective_basic_docker,
        }
    }

    # Controls-only ingest skips pre-ingest header validation because control samples are reference materials
    # (e.g. NA12878) not subject to DRAGEN versioning or ReblockGVCF requirements.
    Boolean effective_validate_vcf_headers = validate_vcf_headers && !samples_are_controls

    if (!effective_validate_vcf_headers && !load_vcf_headers && !load_vet_and_ref_ranges) {
        call Utils.TerminateWorkflow as MustLoadAtLeastOneThing {
            input:
                message = "GvsBulkIngestGenomes called with validate_vcf_headers, load_vcf_headers, and load_vet_and_ref_ranges all set to false",
                basic_docker = effective_basic_docker,
        }
    }

    if (!defined(bulk_ingest_fofn)) {
        # If the user has not specified a FOFN, we will generate one from the data table.
        call GenerateImportFofnFromDataTable {
            input:
                variants_docker = effective_variants_docker,
                sample_set_name = sample_set_name,
                data_table_name = data_table_name,
                user_defined_sample_id_column_name = sample_id_column_name, ## NOTE: the user needs to define this, or it will default to the <entity>_id column
                vcf_files_column_name = vcf_files_column_name,
                vcf_index_files_column_name = vcf_index_files_column_name,
                workspace_id = effective_workspace_id,
                workspace_bucket = effective_workspace_bucket,
        }
    }

    call SplitBulkImportFofn {
        input:
            import_fofn = select_first([GenerateImportFofnFromDataTable.output_fofn, bulk_ingest_fofn]),
            basic_docker = effective_basic_docker,
    }

    call AssignIds.GvsAssignIds as AssignIds {
        input:
            git_branch_or_tag = git_branch_or_tag,
            git_hash = effective_git_hash,
            dataset_name = dataset_name,
            project_id = project_id,
            external_sample_names = SplitBulkImportFofn.sample_name_fofn,
            load_vcf_headers = (load_vcf_headers || effective_validate_vcf_headers),
            load_vet_and_ref_ranges = load_vet_and_ref_ranges,
            cloud_sdk_docker = effective_cloud_sdk_docker,
            use_compressed_references = use_compressed_references,
            samples_are_controls = samples_are_controls,
    }

    # Separate scratch Parquet directories for headers and data passes to avoid prefix collisions.
    # The data pass unconditionally uses /data so that recursive listing never picks up residual
    # files from a prior header pass, even if validate_vcf_headers is false.
    String? headers_parquet_output_gcs_dir = if (defined(parquet_output_gcs_dir) && effective_validate_vcf_headers) then select_first([parquet_output_gcs_dir]) + "/headers" else none_string
    String? data_parquet_output_gcs_dir = if defined(parquet_output_gcs_dir) then select_first([parquet_output_gcs_dir]) + "/data" else none_string

    # VS-1966 / VS-1995: If validate_vcf_headers is true, run an initial headers-only ingest pass,
    # validate the ingested headers, and generate a report. If validation fails, the workflow halts
    # fast before any expensive vet/ref data ingest.
    if (effective_validate_vcf_headers) {
        call ImportGenomes.GvsImportGenomes as ImportHeaders {
            input:
                go = AssignIds.done,
                git_branch_or_tag = git_branch_or_tag,
                git_hash = effective_git_hash,
                dataset_name = dataset_name,
                project_id = project_id,
                external_sample_names = SplitBulkImportFofn.sample_name_fofn,
                num_samples = SplitBulkImportFofn.sample_num,
                input_vcfs = SplitBulkImportFofn.vcf_file_name_fofn,
                input_vcf_indexes = SplitBulkImportFofn.vcf_index_file_name_fofn,
                reference_name = reference_name,
                interval_list = interval_list,
                load_data_scatter_width = load_data_scatter_width,
                load_data_maxretries_override = load_data_maxretries_override,
                load_data_preemptible_override = load_data_preemptible_override,
                basic_docker = effective_basic_docker,
                cloud_sdk_docker = effective_cloud_sdk_docker,
                variants_docker = effective_variants_docker,
                gatk_docker = effective_gatk_docker,
                load_data_gatk_override = gatk_override,
                drop_state = drop_state,
                billing_project_id = billing_project_id,
                use_compressed_references = use_compressed_references,
                load_vet_and_ref_ranges = false,
                load_vcf_headers = true,
                is_rate_limited_beta_customer = tighter_gcp_quotas,
                use_parquet_ingest = use_parquet_ingest,
                parquet_output_gcs_dir = headers_parquet_output_gcs_dir,
                is_wgs = is_wgs,
        }

        call ValidateVcfHeaders.GvsValidateVcfHeaders as ValidateHeaders {
            input:
                go = ImportHeaders.done,
                dataset_name = dataset_name,
                project_id = project_id,
                sample_names_file = SplitBulkImportFofn.sample_name_fofn,
                expected_dragen_version = expected_dragen_version,
                require_reblocking = require_reblocking,
                fail_on_validation_errors = fail_on_validation_errors,
                git_branch_or_tag = git_branch_or_tag,
                variants_docker = effective_variants_docker,
                basic_docker = effective_basic_docker,
        }
    }

    # Data ingest pass:
    # - If validate_vcf_headers is true, gated on ValidateHeaders.done to load vet/ref ranges. If
    #   load_vet_and_ref_ranges is false, this is intentionally skipped (headers-only pass).
    # - If validate_vcf_headers is false, runs directly gated on AssignIds.done (direct import).
    Boolean should_run_data_import = if (effective_validate_vcf_headers) then load_vet_and_ref_ranges else (load_vet_and_ref_ranges || load_vcf_headers)
    if (should_run_data_import) {
        call ImportGenomes.GvsImportGenomes as ImportGenomesData {
            input:
                go = select_first([ValidateHeaders.done, AssignIds.done]),
                git_branch_or_tag = git_branch_or_tag,
                git_hash = effective_git_hash,
                dataset_name = dataset_name,
                project_id = project_id,
                external_sample_names = SplitBulkImportFofn.sample_name_fofn,
                num_samples = SplitBulkImportFofn.sample_num,
                input_vcfs = SplitBulkImportFofn.vcf_file_name_fofn,
                input_vcf_indexes = SplitBulkImportFofn.vcf_index_file_name_fofn,
                reference_name = reference_name,
                interval_list = interval_list,
                load_data_scatter_width = load_data_scatter_width,
                load_data_maxretries_override = load_data_maxretries_override,
                load_data_preemptible_override = load_data_preemptible_override,
                basic_docker = effective_basic_docker,
                cloud_sdk_docker = effective_cloud_sdk_docker,
                variants_docker = effective_variants_docker,
                gatk_docker = effective_gatk_docker,
                load_data_gatk_override = gatk_override,
                drop_state = drop_state,
                billing_project_id = billing_project_id,
                use_compressed_references = use_compressed_references,
                load_vet_and_ref_ranges = load_vet_and_ref_ranges,
                load_vcf_headers = if (!effective_validate_vcf_headers) then load_vcf_headers else false,
                is_rate_limited_beta_customer = tighter_gcp_quotas,
                use_parquet_ingest = use_parquet_ingest,
                parquet_output_gcs_dir = data_parquet_output_gcs_dir,
                is_wgs = is_wgs,
        }
    }

    output {
        Boolean done = true
        String recorded_git_hash = effective_git_hash
        Boolean used_tighter_gcp_quotas = tighter_gcp_quotas
        # Optional because the ValidateHeaders call is conditional (validate_vcf_headers): Cromwell
        # surfaces an un-run conditional call's outputs as None, so consumers must treat these as
        # optional (defined()/select_first). This is the same pattern GvsValidateVDS uses.
        Boolean? vcf_headers_validation_passed = ValidateHeaders.validation_passed
        File? vcf_headers_validation_report = ValidateHeaders.validation_report
        String? vcf_headers_validation_report_contents = ValidateHeaders.validation_report_contents
    }
}


task GenerateImportFofnFromDataTable {
    ## In order to get the names of the columns with the GVCF and GVCF Index file paths, without requiring that the user input it manually, we apply heuristics
    input {
        String data_table_name ## NOTE: if not specified by the user, this has been set to "sample"
        String? sample_set_name
        String? user_defined_sample_id_column_name
        String? vcf_files_column_name
        String? vcf_index_files_column_name
        String workspace_id
        String workspace_bucket
        String variants_docker
    }
    meta {
        # Do not cache as this relies heavily on workspace state.
        volatile: true
    }

    ## set some output files
    String vcf_files_column_name_output_file = "vcf_files_column_name.txt"
    String vcf_index_files_column_name_output_file = "vcf_index_files_column_name.txt"

    String entity_id = data_table_name + "_id"

    String output_fofn_name = "output.tsv"
    String error_file_name = "errors.txt"

    String workspace_name_output = "workspace_name.txt"
    String workspace_namespace_output = "workspace_namespace.txt"

    String sample_name_column = if (defined(user_defined_sample_id_column_name)) then select_first([user_defined_sample_id_column_name]) else entity_id

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        export WORKSPACE_ID="~{workspace_id}"
        export WORKSPACE_BUCKET="~{workspace_bucket}"

        # Hit rawls with the workspace ID

        python3 /app/get_workspace_name_for_import.py \
            --workspace_id ${WORKSPACE_ID} \
            --workspace_name_output '~{workspace_name_output}' \
            --workspace_namespace_output '~{workspace_namespace_output}'

        export WORKSPACE_NAME="$(cat '~{workspace_name_output}')"
        export WORKSPACE_NAMESPACE="$(cat '~{workspace_namespace_output}')"

        # Get a list of all columns in the table. Apply basic heuristics to write the resulting vcf_files_column_name and vcf_index_files_column_name.

        python3 /app/get_columns_for_import.py \
             ~{"--user_defined_sample_id " + user_defined_sample_id_column_name} \
             ~{"--entity_set_name " + sample_set_name} \
             ~{"--user_defined_vcf " + vcf_files_column_name} \
             ~{"--user_defined_index " + vcf_index_files_column_name} \
            --entity_type ~{data_table_name} \
            --vcf_output ~{vcf_files_column_name_output_file} \
            --vcf_index_output ~{vcf_index_files_column_name_output_file}

        if [[ -z "~{vcf_files_column_name}" ]]
        then
            export VCF_COLUMN_NAME="$(cat ~{vcf_files_column_name_output_file})"
        else
            export VCF_COLUMN_NAME="~{vcf_files_column_name}"
        fi

        if [[ -z "~{vcf_index_files_column_name}" ]]
        then
            export VCF_INDEX_COLUMN_NAME="$(cat ~{vcf_index_files_column_name_output_file})"
        else
            export VCF_INDEX_COLUMN_NAME="~{vcf_index_files_column_name}"
        fi

        export GOOGLE_PROJECT="${WORKSPACE_NAMESPACE}"
        python3 /app/generate_fofn_for_import.py \
            --data-table-name ~{data_table_name} \
            --sample-id-column-name ~{sample_name_column} \
            --vcf-files-column-name "${VCF_COLUMN_NAME}" \
            --vcf-index-files-column-name "${VCF_INDEX_COLUMN_NAME}" \
            ~{"--sample-set-name " + sample_set_name} \
            --output-file-name ~{output_fofn_name} \
            --error-file-name ~{error_file_name}

        if [ -s ~{error_file_name} ]; then
            echo ""
            echo "-------- the following issues were found with the sample data, no samples were ingested in this run --------"
            cat ~{error_file_name}
            echo ""
            exit 1
        fi

    >>>

    runtime {
        docker: variants_docker
        memory: "3 GB"
        disks: "local-disk 200 HDD"
        cpu: 1
    }

    output {
        File output_fofn = output_fofn_name
    }
}


task SplitBulkImportFofn {
    input {
        File import_fofn
        String basic_docker
    }

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        cut -f 1 ~{import_fofn} > sample_names.txt
        cut -f 2 ~{import_fofn} > vcf_file_names.txt
        cut -f 3 ~{import_fofn} > vcf_index_file_names.txt
        wc -l < ~{import_fofn} > sample_num.txt
    >>>

    runtime {
        docker: basic_docker
        memory: "3 GB"
        disks: "local-disk 200 HDD"
        cpu: 1
    }

    output {
        File sample_name_fofn = "sample_names.txt"
        File vcf_file_name_fofn = "vcf_file_names.txt"
        File vcf_index_file_name_fofn = "vcf_index_file_names.txt"
        Int sample_num = read_int("sample_num.txt")
    }
}
