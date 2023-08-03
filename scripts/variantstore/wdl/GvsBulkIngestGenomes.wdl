version 1.0

import "GvsUtils.wdl" as Utils
import "GvsAssignIds.wdl" as AssignIds
import "GvsImportGenomes.wdl" as ImportGenomes


workflow GvsBulkIngestGenomes {
    input {
        # Begin GenerateImportFofnFromDataTable
        # for now set the entity type names with a default
        String data_table_name = "sample" ## Note that it is possible an advanced user has a different name for the table. We could glean some information from the sample_set name if that has been defined, but this has not, and try to use that information instead of simply using the default "sample"
        String? sample_id_column_name ## Note that a column WILL exist that is the <entity>_id from the table name. However, some users will want to specify an alternate column for the sample_name during ingest
        String? vcf_files_column_name
        String? vcf_index_files_column_name
        String? sample_set_name ## NOTE: currently we only allow the loading of one sample set at a time
        # End GenerateImportFofnFromDataTable

        # Begin GvsAssignIds
        String dataset_name
        String project_id

        File? gatk_override
        # End GvsAssignIds

        # Begin GvsImportGenomes
        File interval_list = "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list"

        # set to "NONE" to ingest all the reference data into GVS for VDS (instead of VCF) output
        String drop_state = "NONE"

        # The larger the `load_data_batch_size` the greater the probability of preemptions and non-retryable BigQuery errors,
        # so if specifying `load_data_batch_size`, adjust preemptible and maxretries accordingly. Or just take the defaults, as those should work fine in most cases.
        Int? load_data_batch_size
        Int? load_data_preemptible_override
        Int? load_data_maxretries_override
        # End GvsImportGenomes
    }

    parameter_meta {
        data_table_name: "The name of the data table; This table holds the GVCFs to be ingested; `sample` is the default."
        sample_id_column_name: "The column that will be used for the sample name / id in GVS; the <data_table_name>_id will be used as the default"
        vcf_files_column_name: "The column that supplies the path for the GVCFs to be ingested. If not specified, the workflow will attempt to derive the column name."
        vcf_index_files_column_name: "The column that supplies the path for the GVCF index files to be ingested. If not specified, the workflow will attempt to derive the column name."
        sample_set_name: "The recommended way to load samples; Sample sets must be created by the user. If no sample_set_name is specified, all samples will be loaded into GVS"
    }

    ## Start off by getting the Workspace ID to query for more information
    call GetWorkspaceId

    ## Next, use the workspace ID to get the Workspace Name
    call GetWorkspaceName {
        input:
            workspace_id = GetWorkspaceId.workspace_id,
            workspace_bucket = GetWorkspaceId.workspace_bucket,
    }


    call GetColumnNames {
        input:
            workspace_name = GetWorkspaceName.workspace_name,
            workspace_namespace = GetWorkspaceName.workspace_namespace,
            sample_set_name = sample_set_name,
            data_table_name = data_table_name,
            user_defined_sample_id_column_name = sample_id_column_name, ## NOTE: the user needs to define this, or it will default to the <entity>_id column
            vcf_files_column_name = vcf_files_column_name,
            vcf_index_files_column_name = vcf_index_files_column_name,
    }

    call GenerateImportFofnFromDataTable {
        input:
            google_project_id = GetWorkspaceName.workspace_namespace,
            workspace_name = GetWorkspaceName.workspace_name,
            workspace_namespace = GetWorkspaceName.workspace_namespace,
            workspace_bucket = GetWorkspaceId.workspace_bucket,
            samples_table_name = GetColumnNames.data_table,
            sample_id_column_name = GetColumnNames.sample_name_column,  ## NOTE: if no sample_id_column_name has been specified, this is now the <entity>_id column
            vcf_files_column_name = GetColumnNames.vcf_files_column_name_output,
            vcf_index_files_column_name = GetColumnNames.vcf_index_files_column_name_output,
            sample_set_name = sample_set_name,
    }

    call SplitBulkImportFofn {
        input:
            import_fofn = GenerateImportFofnFromDataTable.output_fofn,
    }

    call AssignIds.GvsAssignIds as AssignIds {
        input:
            dataset_name = dataset_name,
            project_id = project_id,
            external_sample_names = SplitBulkImportFofn.sample_name_fofn,
            samples_are_controls = false
    }

    call ImportGenomes.GvsImportGenomes as ImportGenomes {
        input:
            go = AssignIds.done,
            dataset_name = dataset_name,
            project_id = project_id,
            external_sample_names = SplitBulkImportFofn.sample_name_fofn,
            num_samples = SplitBulkImportFofn.sample_num,
            input_vcfs = SplitBulkImportFofn.vcf_file_name_fofn,
            input_vcf_indexes = SplitBulkImportFofn.vcf_index_file_name_fofn,
            interval_list = interval_list,

            # The larger the `load_data_batch_size` the greater the probability of preemptions and non-retryable
            # BigQuery errors so if specifying this adjust preemptible and maxretries accordingly. Or just take the defaults,
            # those should work fine in most cases.
            load_data_batch_size = load_data_batch_size,
            load_data_maxretries_override = load_data_maxretries_override,
            load_data_preemptible_override = load_data_preemptible_override,
            load_data_gatk_override = gatk_override,
            drop_state = drop_state,
    }

    output {
        Boolean done = true
    }
}

task GetWorkspaceId {
    ## In order to get the ID and bucket of the workspace, without requiring that the user input it manually, we check the delocalization script
    meta {
        volatile: true # always run this when asked otherwise you can get a previously run workspace
    }
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        # Sniff the workspace bucket out of the delocalization script and extract the workspace id from that.
        sed -n -E 's!.*gs://fc-(secure-)?([^\/]+).*!\2!p' /cromwell_root/gcs_delocalization.sh | sort -u > workspace_id.txt
        sed -n -E 's!.*gs://(fc-(secure-)?[^\/]+).*!\1!p' /cromwell_root/gcs_delocalization.sh | sort -u > workspace_bucket.txt
    >>>

    runtime {
        docker: "ubuntu:latest"
    }

    output {
        String workspace_id = read_string("workspace_id.txt")
        String workspace_bucket = read_string("workspace_bucket.txt")
    }
}

task GetWorkspaceName {
    ## In order to get the name and namespace of the workspace, without requiring that the user input it manually, we hit rawls directly for the exact name
    input {
        String workspace_id
        String workspace_bucket
    }

    String workspace_name_output = "workspace_name.txt"
    String workspace_namespace_output = "workspace_namespace.txt"

    command <<<
        # Hit rawls with the workspace ID

        export WORKSPACE_BUCKET='~{workspace_bucket}'

        python3 /app/get_workspace_name_for_import.py \
        --workspace_id ~{workspace_id} \
        --workspace_name_output ~{workspace_name_output} \
        --workspace_namespace_output ~{workspace_namespace_output} \

    >>>
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:2023-07-24-alpine-bd6b9d62e"
        memory: "3 GB"
        disks: "local-disk 10 HDD"
        cpu: 1
    }

    output {
        String workspace_name = read_string(workspace_name_output)
        String workspace_namespace = read_string(workspace_namespace_output)
    }
}

task GetColumnNames {
    ## In order to get the names of the columns with the GVCF and GVCF Index file paths, without requiring that the user input it manually, we apply heuristics
    input {
        String workspace_name
        String workspace_namespace
        String data_table_name ## NOTE: if not specified by the user, this has been set to "sample"
        String? sample_set_name
        String? user_defined_sample_id_column_name
        String? vcf_files_column_name
        String? vcf_index_files_column_name
    }

    ## set some output files
    String vcf_files_column_name_output_file = "vcf_files_column_name.txt"
    String vcf_index_files_column_name_output_file = "vcf_index_files_column_name.txt"

    String entity_id = data_table_name + "_id"

     command <<<
        # Get a list of all columns in the table. Apply basic heuristics to write the resulting vcf_files_column_name and vcf_index_files_column_name.

        export WORKSPACE_NAMESPACE='~{workspace_namespace}'
        export WORKSPACE_NAME='~{workspace_name}'

        python3 /app/get_columns_for_import.py \
           ~{"--user_defined_sample_id " + user_defined_sample_id_column_name} \
           ~{"--entity_set_name " + sample_set_name} \
           ~{"--user_defined_vcf " + vcf_files_column_name} \
           ~{"--user_defined_index " + vcf_index_files_column_name} \
          --entity_type ~{data_table_name} \
          --vcf_output ~{vcf_files_column_name_output_file} \
          --vcf_index_output ~{vcf_index_files_column_name_output_file}

    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:2023-07-24-alpine-bd6b9d62e"
        memory: "3 GB"
        disks: "local-disk 10 HDD"
        cpu: 1
    }

    output {
        String vcf_files_column_name_output = if (defined(vcf_files_column_name)) then select_first([vcf_files_column_name]) else read_string(vcf_files_column_name_output_file)
        String vcf_index_files_column_name_output = if (defined(vcf_index_files_column_name)) then select_first([vcf_index_files_column_name]) else read_string(vcf_index_files_column_name_output_file)
        String sample_name_column = if (defined(user_defined_sample_id_column_name)) then select_first([user_defined_sample_id_column_name]) else entity_id
        String data_table = data_table_name
    }
}

task SplitBulkImportFofn {
    input {
        File import_fofn
    }

    command <<<
        set -o errexit -o nounset -o xtrace -o pipefail
        PS4='\D{+%F %T} \w $ '

        cut -f 1 ~{import_fofn} > sample_names.txt
        cut -f 2 ~{import_fofn} > vcf_file_names.txt
        cut -f 3 ~{import_fofn} > vcf_index_file_names.txt
        wc -l < ~{import_fofn} > sample_num.txt
    >>>

    runtime {
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:441.0.0-alpine"
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


task GenerateImportFofnFromDataTable {
    ## In order to get the <entity>_ids in the sample_set for an inclusion list, we use Terra Notebook Utils
    ## This also allows us to validate that the requested sample_set exists
    input {
        String google_project_id
        String workspace_name
        String workspace_namespace
        String workspace_bucket
        String samples_table_name
        String sample_id_column_name ## NOTE: if the user has specified a different sample name column for GVS, it needs to be used independently of the sample_set info
        String vcf_files_column_name
        String vcf_index_files_column_name
        String? sample_set_name
    }

    String output_fofn_name = "output.tsv"
    String error_file_name = "errors.txt"

    command <<<
        set -o errexit -o nounset -o pipefail
        PS4='\D{+%F %T} \w $ '

        export GOOGLE_PROJECT='~{google_project_id}'
        export WORKSPACE_NAMESPACE='~{workspace_namespace}'
        export WORKSPACE_NAME='~{workspace_name}'
        export WORKSPACE_BUCKET='~{workspace_bucket}'

        python3 /app/generate_fofn_for_import.py \
            --data-table-name ~{samples_table_name} \
            --sample-id-column-name ~{sample_id_column_name} \
            --vcf-files-column-name ~{vcf_files_column_name} \
            --vcf-index-files-column-name ~{vcf_index_files_column_name} \
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
        docker: "us.gcr.io/broad-dsde-methods/variantstore:2023-07-24-alpine-bd6b9d62e"
        memory: "3 GB"
        disks: "local-disk 200 HDD"
        cpu: 1
    }

    output {
        File output_fofn = output_fofn_name
    }
}
