version 1.0

workflow GvsCallsetCost {
    input {
#        String project_id
#        String dataset_name
        String workspace_namespace
        String workspace_name
#        String callset_name
        Array[String] excluded_submission_ids = []
    }

    call WorkflowComputeCosts {
        input:
            workspace_namespace = workspace_namespace,
            workspace_name = workspace_name,
            excluded_submission_ids = excluded_submission_ids
    }

#    call BigQueryWriteAPICost {
#        input:
#            project_id = project_id,
#            dataset_name = dataset_name
#    }
#
#    call BigQueryScannedCost {
#        input:
#            project_id = project_id,
#            dataset_name = dataset_name,
#            callset_name = callset_name
#    }
#
#    call BigQueryStorageAPIScannedCost {
#        input:
#            project_id = project_id,
#            dataset_name = dataset_name,
#            callset_name = callset_name
#    }

    output {
        File workflow_compute_costs = WorkflowComputeCosts.costs
        File workflow_compute_costs_log = WorkflowComputeCosts.log
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

#task BigQueryWriteAPICost {
#    meta {
#        description: "Estimate GvsImportGenomes use of the BQ Write API via core storage costs from the sizes of vet_% and ref_ranges_% tables."
#        volatile: true
#    }
#
#    input {
#        String project_id
#        String dataset_name
#    }
#    command <<<
#    >>>
#
#    runtime {
#        docker: ""
#    }
#
#    output {
#        Float vet_gib = read_float("")
#        Float ref_ranges_gib = read_float("")
#        Float import_genomes_cost = 3
#    }
#}

#task BigQueryScannedCost {
#    meta {
#        description: "Determine BigQuery scanned cost for GVSCreateAltAllele, GVSCreateFilterSet, and GVSPrepareRanges"
#        volatile: true
#    }
#
#    input {
#        String project_id
#        String dataset_name
#        String callset_name
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
#        String callset_name
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
