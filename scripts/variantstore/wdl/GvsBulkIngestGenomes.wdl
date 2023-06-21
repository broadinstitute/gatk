version 1.0

import "GvsAssignIds.wdl" as AssignIds
import "GvsImportGenomes.wdl" as ImportGenomes


workflow GvsBulkIngestGenomes {
    input {
        String dataset_name
        String project_id
        String call_set_identifier
        String samples_table_name
        String sample_id_column_name
        File interval_list
        String drop_state

        String? vcf_files_column_name
        String? vcf_index_files_column_name
        Int? load_data_batch_size
        Int? load_data_preemptible_override
        Int? load_data_maxretries_override

        File? gatk_override
    }

    call GetWorkspaceId

    call GetWorkspaceName {
        input:
            workspace_id = GetWorkspaceId.workspace_id,
            workspace_bucket = GetWorkspaceId.workspace_bucket,
    }

    if (!defined(vcf_files_column_name) || !defined(vcf_index_files_column_name)) {
        call GetColumnNames {
            input:
                workspace_id = GetWorkspaceId.workspace_id,
                workspace_name = GetWorkspaceName.workspace_name,
                workspace_namespace = GetWorkspaceName.workspace_namespace,
                samples_table_name = samples_table_name,
                sample_id_column_name = sample_id_column_name,
        }
    }


    call GenerateFOFNsFromDataTables {
        input:
            project_id = GetWorkspaceName.workspace_namespace,
            workspace_name = GetWorkspaceName.workspace_name,
            workspace_namespace = GetWorkspaceName.workspace_namespace,
            workspace_bucket  = GetWorkspaceId.workspace_bucket,
            samples_table_name = samples_table_name,
            sample_id_column_name = sample_id_column_name,
            vcf_files_column_name = select_first([vcf_files_column_name, GetColumnNames.vcf_files_column_name]),
            vcf_index_files_column_name = select_first([vcf_index_files_column_name, GetColumnNames.vcf_index_files_column_name]),
    }

    call AssignIds.GvsAssignIds as AssignIds {
        input:
            dataset_name = dataset_name,
            project_id = project_id,
            external_sample_names = read_lines(GenerateFOFNsFromDataTables.sampleFOFN),
            samples_are_controls = false
    }

    call ImportGenomes.GvsImportGenomes as ImportGenomes {
        input:
            go = AssignIds.done,
            dataset_name = dataset_name,
            project_id = project_id,
            external_sample_names = read_lines(GenerateFOFNsFromDataTables.sampleFOFN),
            input_vcfs = read_lines(GenerateFOFNsFromDataTables.vcfFOFN),
            input_vcf_indexes = read_lines(GenerateFOFNsFromDataTables.vcfIndexFOFN),

            interval_list = interval_list,

            # The larger the `load_data_batch_size` the greater the probability of preemptions and non-retryable
            # BigQuery errors so if specifying this adjust preemptible and maxretries accordingly. Or just take the defaults,
            # those should work fine in most cases.
            load_data_batch_size = load_data_batch_size,
            load_data_maxretries_override = load_data_maxretries_override,
            load_data_preemptible_override = load_data_preemptible_override,
            load_data_gatk_override = gatk_override,
    }

    output {
        Boolean done = ImportGenomes.done
    }
}

task GetColumnNames {
    input {
        String workspace_id
        String workspace_name
        String workspace_namespace
        String samples_table_name
        String sample_id_column_name
        String? vcf_files_column_name
        String? vcf_index_files_column_name
    }
    ## set some output files
    String vcf_files_column_name_output = "vcf_files_column_name.txt"
    String vcf_index_files_column_name_output = "vcf_index_files_column_name.txt"


    command <<<
        # Get a list of all columns in the table

        export WORKSPACE_NAMESPACE='~{workspace_namespace}'
        export WORKSPACE_NAME='~{workspace_name}'

        python3 /app/get_columns_for_import.py \
          --workspace_id ~{workspace_id} \
          --vcf_output ~{vcf_files_column_name_output} \
          --vcf_index_output ~{vcf_index_files_column_name_output} \


    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:2023-06-13-alpine"
        memory: "3 GB"
        disks: "local-disk 10 HDD"
        cpu: 1
    }

    output {
        String vcf_files_column_name = if (defined(vcf_files_column_name)) then select_first([vcf_files_column_name]) else read_string(vcf_files_column_name_output)
        String vcf_index_files_column_name = if (defined(vcf_index_files_column_name)) then select_first([vcf_index_files_column_name]) else read_string(vcf_index_files_column_name_output)
    }
}

task GetWorkspaceName {

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
        docker: "us.gcr.io/broad-dsde-methods/variantstore:2023-06-13-alpine"
        memory: "3 GB"
        disks: "local-disk 10 HDD"
        cpu: 1
    }

    output {
        String workspace_name = read_string(workspace_name_output)
        String workspace_namespace = read_string(workspace_namespace_output)
    }
}


# This copies the entirety of the data tables over.  We'll want a different task to make FOFNs for just a sample set
task GenerateFOFNsFromDataTables {
    input {
        String project_id
        String workspace_name
        String workspace_namespace
        String workspace_bucket
        String? samples_table_name
        String? sample_id_column_name
        String? vcf_files_column_name
        String? vcf_index_files_column_name
    }

    ## TODO I dont love that we are hardcoding them here and in the python--they need to be params!
    String sample_names_file_name = "sample_names.txt"
    String vcf_files_name = "vcf_files.txt"
    String vcf_index_files_name = "vcf_index_files.txt"
    String error_file_name = "errors.txt"

    command <<<
        set -o errexit -o nounset -o xtrace -o pipefail

        export GOOGLE_PROJECT='~{project_id}'
        export WORKSPACE_NAMESPACE='~{workspace_namespace}'
        export WORKSPACE_NAME='~{workspace_name}'
        export WORKSPACE_BUCKET='~{workspace_bucket}'

        python3 /app/generate_FOFNs_for_import.py \
        --data_table_name ~{samples_table_name} \
        --sample_id_column_name ~{sample_id_column_name} \
        --vcf_files_column_name ~{vcf_files_column_name} \
        --vcf_index_files_column_name ~{vcf_index_files_column_name}

    >>>
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:2023-06-13-alpine"
        memory: "3 GB"
        disks: "local-disk 200 HDD"
        cpu: 1
    }

    output {
        File sampleFOFN = sample_names_file_name
        File vcfFOFN = vcf_files_name
        File vcfIndexFOFN = vcf_index_files_name
        File errors = error_file_name
    }
}

task GetWorkspaceId {
    meta {
        volatile: true # always run this when asked otherwise you can get a previously run workspace!!!!
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

