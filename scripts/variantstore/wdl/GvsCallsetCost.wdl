version 1.0

import "GvsUtils.wdl" as Utils


workflow GvsCallsetCost {
    input {
        String project_id
        String dataset_name
        String workspace_namespace
        String workspace_name
        String call_set_identifier
        Array[String] excluded_submission_ids = []
        String? cloud_sdk_docker
        String? variants_docker
    }

    if (!defined(cloud_sdk_docker) || !defined(variants_docker)) {
        call Utils.GetToolVersions
    }

    String effective_cloud_sdk_docker = select_first([cloud_sdk_docker, GetToolVersions.cloud_sdk_docker])
    String effective_variants_docker = select_first([variants_docker, GetToolVersions.variants_docker])

    call WorkflowComputeCosts {
        input:
            workspace_namespace = workspace_namespace,
            workspace_name = workspace_name,
            excluded_submission_ids = excluded_submission_ids,
            variants_docker = effective_variants_docker,
    }

    call CoreStorageModelSizes {
        input:
            project_id = project_id,
            dataset_name = dataset_name,
            cloud_sdk_docker = effective_cloud_sdk_docker,
    }

    call ReadCostObservabilityTable {
        input:
            project_id = project_id,
            dataset_name = dataset_name,
            call_set_identifier = call_set_identifier,
            cloud_sdk_docker = effective_cloud_sdk_docker,
    }

    output {
        File workflow_compute_costs = WorkflowComputeCosts.costs
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
        String variants_docker
    }

    Array[String] excluded_ids = prefix('--exclude ', excluded_submission_ids)

    command <<<
        python3 /app/workflow_compute_costs.py \
            --workspace_namespace '~{workspace_namespace}' \
            --workspace_name '~{workspace_name}' \
            ~{sep=' ' excluded_ids} \
            > costs_by_workflow.json
    >>>

    runtime {
        docker: variants_docker
    }

    output {
        File costs = "costs_by_workflow.json"
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
        String cloud_sdk_docker
    }
    command <<<

        get_billable_bytes_in_gib() {
            local table_pattern="$1"
            local output_file_name="$2"

            bq --apilog=false query --project_id='~{project_id}' --format=csv --use_legacy_sql=false \
                "SELECT round(sum(total_billable_bytes) / (1024*1024*1024),2) \
                    FROM \`~{project_id}.~{dataset_name}.INFORMATION_SCHEMA.PARTITIONS\` \
                    WHERE table_name LIKE '${table_pattern}'" | tail -1 > ${output_file_name}
        }

        get_billable_bytes_in_gib "vet_%"        vet_gib.txt
        get_billable_bytes_in_gib "ref_ranges_%" ref_ranges_gib.txt
        get_billable_bytes_in_gib "alt_allele"   alt_allele_gib.txt
    >>>
    runtime {
        docker: cloud_sdk_docker
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
        String cloud_sdk_docker
    }
    command <<<
        bq --apilog=false query --project_id='~{project_id}' --format=prettyjson --use_legacy_sql=false \
            "SELECT step, event_key, round(sum(event_bytes) / (1024*1024*1024), 2) AS sum_event_gibibytes \
                FROM \`~{project_id}.~{dataset_name}.cost_observability\` \
                WHERE call_set_identifier = '~{call_set_identifier}' GROUP BY step, event_key ORDER BY step" \
            > cost_observability.json
    >>>
    runtime {
        docker: cloud_sdk_docker
    }
    output {
        File cost_observability = "cost_observability.json"
    }
}
