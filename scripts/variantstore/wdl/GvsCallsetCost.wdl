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

    call ReadCostObservabilityTable {
        input:
            project_id = project_id,
            dataset_name = dataset_name,
            call_set_identifier = call_set_identifier,
    }

    output {
        File workflow_compute_costs = WorkflowComputeCosts.costs
        File workflow_compute_costs_log = WorkflowComputeCosts.log
        String vet_gib = CoreStorageModelSizes.vet_gib
        String ref_ranges_gib = CoreStorageModelSizes.ref_ranges_gib
        String alt_allele_gib = CoreStorageModelSizes.alt_allele_gib
        File cost_observability = ReadCostObservabilityTable.cost_observability
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
        docker: "us.gcr.io/broad-dsde-methods/variantstore:ah_var_store_2022_07_08"
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

task ReadCostObservabilityTable {
    meta {
        description: "Read data from cost_observability table for the specified parameters."
        # Definitely don't cache this, the values will change while the inputs to this task will not!
        volatile: true
    }
    input {
        String project_id
        String dataset_name
        String call_set_identifier
    }
    command <<<
        bq query --location=US --project_id='~{project_id}' --format=prettyjson --use_legacy_sql=false \
            "SELECT step, event_key, round(sum(event_bytes) / (1024*1024*1024), 2) AS sum_event_gibibytes \
                FROM \`~{project_id}.~{dataset_name}.cost_observability\` \
                WHERE call_set_identifier = '~{call_set_identifier}' GROUP BY step, event_key ORDER BY step" \
            > cost_observability.json
    >>>
    runtime {
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:390.0.0"
    }
    output {
        File cost_observability = "cost_observability.json"
    }
}
