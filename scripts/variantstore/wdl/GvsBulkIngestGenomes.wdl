version 1.0

import "GvsUtils.wdl" as Utils
import "GvsPrepareBulkImport.wdl" as PrepareBulkImport
import "GvsAssignIds.wdl" as AssignIds
import "GvsImportGenomes.wdl" as ImportGenomes


workflow GvsBulkIngestGenomes {
    input {
        # Begin GvsPrepareBulkImport
        # for now set the entity type names with a default
        String samples_table_name = "sample"
        String sample_id_column_name = "sample_id"
        String? vcf_files_column_name
        String? vcf_index_files_column_name
        String? sample_set_name ## currently we only allow for one sample set at a time
        # End GvsPrepareBulkImport

        # Begin GvsAssignIds
        String dataset_name
        String bq_project_id
        String call_set_identifier

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

    call GetWorkspaceId

    call GetWorkspaceName {
        input:
            workspace_id = GetWorkspaceId.workspace_id,
            workspace_bucket = GetWorkspaceId.workspace_bucket,
    }

    call GetColumnNames { ## TODO should we even run this at all if we have values for all 4?
        input:
            workspace_id = GetWorkspaceId.workspace_id,
            workspace_name = GetWorkspaceName.workspace_name,
            workspace_namespace = GetWorkspaceName.workspace_namespace,
            samples_table_name = samples_table_name,
            sample_id_column_name = sample_id_column_name,
            vcf_files_column_name = vcf_files_column_name,
            vcf_index_files_column_name = vcf_index_files_column_name,
    }

    call PrepareBulkImport.GvsPrepareBulkImport as PrepareBulkImport {
        input:
            project_id = GetWorkspaceName.workspace_namespace,
            workspace_name = GetWorkspaceName.workspace_name,
            workspace_namespace = GetWorkspaceName.workspace_namespace,
            workspace_bucket = GetWorkspaceId.workspace_bucket,
            samples_table_name = samples_table_name,
            sample_id_column_name = sample_id_column_name,
            vcf_files_column_name = GetColumnNames.vcf_files_column_name,
            vcf_index_files_column_name = GetColumnNames.vcf_index_files_column_name,
            sample_set_name = sample_set_name
    }

    call AssignIds.GvsAssignIds as AssignIds {
        input:
            dataset_name = dataset_name,
            project_id = bq_project_id,
            external_sample_names = read_lines(PrepareBulkImport.sampleFOFN),
            samples_are_controls = false
    }

    call ImportGenomes.GvsImportGenomes as ImportGenomes {
        input:
            go = AssignIds.done,
            dataset_name = dataset_name,
            project_id = bq_project_id,
            external_sample_names = read_lines(PrepareBulkImport.sampleFOFN),
            input_vcfs = read_lines(PrepareBulkImport.vcfFOFN),
            input_vcf_indexes = read_lines(PrepareBulkImport.vcfIndexFOFN),

            interval_list = interval_list,

            # The larger the `load_data_batch_size` the greater the probability of preemptions and non-retryable
            # BigQuery errors so if specifying this adjust preemptible and maxretries accordingly. Or just take the defaults,
            # those should work fine in most cases.
            load_data_batch_size = load_data_batch_size,
            load_data_maxretries_override = load_data_maxretries_override,
            load_data_preemptible_override = load_data_preemptible_override,
            load_data_gatk_override = gatk_override,
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

