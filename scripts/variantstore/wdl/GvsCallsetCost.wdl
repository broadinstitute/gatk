version 1.0

workflow GvsCallsetCost {
    input {
        String project_id
        String dataset_name
        String workspace_namespace
        String workspace_name
        String call_set_identifier
        Array[String] excluded_submission_ids = []
    }

    call WorkflowComputeCosts {
        input:
            workspace_namespace = workspace_namespace,
            workspace_name = workspace_name,
            excluded_submission_ids = excluded_submission_ids
    }

    call CoreStorageModelSizes {
        input:
            project_id = project_id,
            dataset_name = dataset_name
    }

#    call BigQueryScannedCost {
#        input:
#            project_id = project_id,
#            dataset_name = dataset_name,
#            call_set_identifier = call_set_identifier
#    }
#
#    call BigQueryStorageAPIScannedCost {
#        input:
#            project_id = project_id,
#            dataset_name = dataset_name,
#            call_set_identifier = call_set_identifier
#    }

    output {
        File workflow_compute_costs = WorkflowComputeCosts.costs
        File workflow_compute_costs_log = WorkflowComputeCosts.log
        String vet_gib = CoreStorageModelSizes.vet_gib
        String ref_ranges_gib = CoreStorageModelSizes.ref_ranges_gib
        String alt_allele_gib = CoreStorageModelSizes.alt_allele_gib
    }
}


task WorkflowComputeCosts {
    meta {
        description: "Calculate workflow compute costs by calling Firecloud APIs for submissions in the specified workspace"
        volatile: true
    }
    input {
        String workspace_namespace
        String workspace_name
        Array[String] excluded_submission_ids
    }

    Array[String] excluded_ids = prefix('--exclude ', excluded_submission_ids)

    command <<<
        python3 /app/workflow_compute_costs.py \
            --workspace_namespace '~{workspace_namespace}' \
            --workspace_name '~{workspace_name}' \
            ~{sep=' ' excluded_ids} \
            > costs_by_workflow.json 2> workflow_compute_costs.log
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:vs_472_workflow_compute_costs-2022-06-21"
    }

    output {
        File costs = "costs_by_workflow.json"
        File log = "workflow_compute_costs.log"
    }
}

task CoreStorageModelSizes {
    meta {
        description: "Read sizes of vet_%, ref_ranges_%, and alt_allele tables from `INFORMATION_SCHEMA.PARTITIONS`."
        # Definitely don't cache this, the values will change while the inputs to this task will not!
        volatile: true
    }
    input {
        String project_id
        String dataset_name
    }
    command <<<

        get_billable_bytes_in_gib() {
            local table_pattern="$1"
            local output_file_name="$2"

            bq query --location=US --project_id='~{project_id}' --format=csv --use_legacy_sql=false \
                "SELECT round(sum(total_billable_bytes) / (1024*1024*1024),2) \
                    FROM \`~{project_id}.~{dataset_name}.INFORMATION_SCHEMA.PARTITIONS\` \
                    WHERE table_name LIKE '${table_pattern}'" | tail -1 > ${output_file_name}
        }

        get_billable_bytes_in_gib "vet_%"        vet_gib.txt
        get_billable_bytes_in_gib "ref_ranges_%" ref_ranges_gib.txt
        get_billable_bytes_in_gib "alt_allele"   alt_allele_gib.txt
    >>>
    runtime {
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:390.0.0"
    }
    output {
        Float vet_gib = read_float("vet_gib.txt")
        Float ref_ranges_gib = read_float("ref_ranges_gib.txt")
        Float alt_allele_gib = read_float("alt_allele_gib.txt")
    }
}

#task BigQueryScannedCost {
#    meta {
#        description: "Determine BigQuery scanned cost for GVSCreateAltAllele, GVSCreateFilterSet, and GVSPrepareRanges"
#        volatile: true
#    }
#
#    input {
#        String project_id
#        String dataset_name
#        String call_set_identifier
#    }
#
#    command <<<
#    >>>
#
#    runtime {
#        docker: ""
#    }
#
#    output {
#        Float create_alt_allele_gib = read_float("")
#        Float create_filter_set_gib = read_float("")
#        Float prepare_ranges_gib    = read_float("")
#        Float cost = 3
#    }
#}
#
#task BigQueryStorageAPIScannedCost {
#    meta {
#        description: "Determine BigQuery Storage API scanned cost for GvsCreateFilterSet and GvsExtractCallset"
#        volatile: true
#    }
#
#    input {
#        String project_id
#        String dataset_name
#        String call_set_identifier
#    }
#
#    command <<<
#    >>>
#
#    runtime {
#        docker: ""
#    }
#
#    output {
#        Float create_filter_set_gib = read_float("")
#        Float extract_callset_gib = read_float("")
#        Float cost = 3
#    }
#}
