version 1.0

import "GvsUtils.wdl" as Utils
import "GvsPrepareBulkImport.wdl" as PrepareBulkImport
import "GvsAssignIds.wdl" as AssignIds
import "GvsImportGenomes.wdl" as ImportGenomes


workflow GvsBulkIngestGenomes {
    input {
        # Begin GvsPrepareBulkImport
        String terra_project_id # TODO isn't this just the project namespace?!?!??!
        String? samples_table_name
        String? sample_id_column_name
        String? vcf_files_column_name
        String? vcf_index_files_column_name
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

    #call GetWorkspaceName {
    #    input:
    #        workspace_id = GetWorkspaceId.workspace_id,
    #        workspace_bucket = GetWorkspaceId.workspace_bucket,
    #        project_id = terra_project_id
    #}

    call GetColumnNames { ## TODO should we even run this at all if we have values for all 4?
        input:
            workspace_id = GetWorkspaceId.workspace_id,
            workspace_name = "gvs-dev", ## TODO hard-coding for testing only --- GetWorkspaceName.workspace_name,
            samples_table_name = samples_table_name,
            sample_id_column_name = sample_id_column_name,
            vcf_files_column_name = vcf_files_column_name,
            vcf_index_files_column_name = vcf_index_files_column_name,
    }
}

    task GetColumnNames {
        input {
            String workspace_id
            String workspace_name
            String? samples_table_name
            String? sample_id_column_name
            String? vcf_files_column_name
            String? vcf_index_files_column_name
        }
        command <<<
            # Get a list of all columns in the table

            # First we will check for the default named columns and make sure that each row has a value

            gsutil cp gs://fc-d5e319d4-b044-4376-afde-22ef0afc4088/get_columns_for_import.py  get_columns_for_import.py
            python get_columns_for_import.py --workspace_id ~{workspace_id} --workspace_name ~{workspace_name}
            ##s set some default vals
            # String samples_table_name = "sample"
            # String sample_id_column_name = "sample_id"
            # String vcf_files_column_name = "hg38_reblocked_v2_vcf"
            # String vcf_index_files_column_name = "hg38_reblocked_v2_vcf_index"


            #curl -X 'GET' \
            #'https://rawls.dsde-prod.broadinstitute.org/api/workspaces/{workspace_id}/{workspace_name}/entities?useCache=true' \
            #-H 'accept: application/json' \
            #-H 'Authorization: Bearer XXX

            #{
            #"participant": {
            #"attributeNames": [],
            #"count": 249045,
            #"idName": "participant_id"
            #},
            #"sample": {
            #"attributeNames": [
            #"blah blah these are the column names we want"
            #],
            #"count": 249049,
            #"idName": "sample_id"
            #},
            #"sample_set": {
            #"attributeNames": [
            #"samples"
            #],
            #"count": 12,
            #"idName": "sample_set_id"
            #}
            #}
        >>>

        runtime {
            docker: "us.gcr.io/broad-dsde-methods/variantstore:2023-03-23-bulk-ingest"
        }

        output {
            Array[String] bag_of_column_names = []
        }
    }

    task GetWorkspaceName {

        input {
            String workspace_id
            String workspace_bucket
            String project_id
        }
        command <<<
            # Hit rawls with the workspace ID <-- this is the optimized version that we need to figure out the auth on

            export GOOGLE_PROJECT='~{project_id}'
            export WORKSPACE_BUCKET='~{workspace_bucket}'

            gsutil cp gs://fc-d5e319d4-b044-4376-afde-22ef0afc4088/get_columns_for_import.py  get_columns_for_import.py
            python get_columns_for_import.py --workspace_id ~{workspace_id} > workspace_name.txt

        >>>
        runtime {
            docker: "us.gcr.io/broad-dsde-methods/variantstore:2023-04-13-alpine"
            memory: "3 GB"
            disks: "local-disk 10 HDD"
            cpu: 1
        }

    output {
        String workspace_name = read_string("workspace_name.txt")
        String workspace_namespace = read_string("workspace_name.txt")
        String google_project_id = read_string("workspace_name.txt")
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

