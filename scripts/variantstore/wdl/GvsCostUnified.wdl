version 1.0

workflow GvsUnifiedCost {
    input {
        String project_id
        String dataset_name
        String workspace_namespace
        String workspace_name
        String callset_name
        Array[String] excluded_submission_ids = []
    }

    call WorkflowComputeCost {
        input:
            workspace_namespace = workspace_namespace,
            workspace_name = workspace_name,
            excluded_submission_ids = excluded_submission_ids
    }

    call BigQueryWriteAPICost {
        input:
            project_id = project_id,
            dataset_name = dataset_name
    }

    call BigQueryScannedCost {
        input:
            project_id = project_id,
            dataset_name = dataset_name,
            callset_name = callset_name
    }

    call BigQueryStorageAPIScannedCost {
        input:
            project_id = project_id,
            dataset_name = dataset_name,
            callset_name = callset_name
    }

    output {
        Float workflow_compute_cost = WorkflowComputeCost.cost
        Float core_storage_model_vet_gib = BigQueryWriteAPICost.vet_gib
        Float core_storage_model_ref_ranges_gib = BigQueryWriteAPICost.ref_ranges_gib
        Float write_api_cost = BigQueryWriteAPICost.import_genomes_cost
        Float scanned_create_alt_allele_gib = BigQueryScannedCost.create_alt_allele_gib
        Float scanned_create_filter_set_gib = BigQueryScannedCost.create_filter_set_gib
        Float scanned_prepare_ranges_gib = BigQueryScannedCost.prepare_ranges_gib
        Float scanned_cost = BigQueryScannedCost.cost
        Float storage_api_scanned_create_filter_set_gib = BigQueryStorageAPIScannedCost.create_filter_set_gib
        Float storage_api_scanned_extract_callset_gib = BigQueryStorageAPIScannedCost.extract_callset_gib
        Float storage_api_scanned_cost = BigQueryStorageAPIScannedCost.cost
        Float total_cost = workflow_compute_cost + write_api_cost + scanned_cost + storage_api_scanned_cost
    }
}

task WorkflowComputeCost {
    meta {
        description: "Calculate workflow compute cost by calling Firecloud APIs for submissions in the specified workspace"
        volatile: true
    }

    input {
        String workspace_namespace
        String workspace_name
        Array[String] excluded_submission_ids
    }

    command <<<
    >>>

    runtime {
        docker: ""
    }

    output {
        Float cost = read_float("")
    }
}

task BigQueryWriteAPICost {
    meta {
        description: "Estimate GvsImportGenomes use of the BQ Write API via core storage costs from the sizes of vet_% and ref_ranges_% tables."
        volatile: true
    }

    input {
        String project_id
        String dataset_name
    }
    command <<<
    >>>

    runtime {
        docker: ""
    }

    output {
        Float vet_gib = read_float("")
        Float ref_ranges_gib = read_float("")
        Float import_genomes_cost = 3
    }
}

task BigQueryScannedCost {
    meta {
        description: "Determine BigQuery scanned cost for GVSCreateAltAllele, GVSCreateFilterSet, and GVSPrepareRanges"
        volatile: true
    }

    input {
        String project_id
        String dataset_name
        String callset_name
    }

    command <<<
    >>>

    runtime {
        docker: ""
    }

    output {
        Float create_alt_allele_gib = read_float("")
        Float create_filter_set_gib = read_float("")
        Float prepare_ranges_gib    = read_float("")
        Float cost = 3
    }
}

task BigQueryStorageAPIScannedCost {
    meta {
        description: "Determine BigQuery Storage API scanned cost for GvsCreateFilterSet and GvsExtractCallset"
        volatile: true
    }

    input {
        String project_id
        String dataset_name
        String callset_name
    }

    command <<<
    >>>

    runtime {
        docker: ""
    }

    output {
        Float create_filter_set_gib = read_float("")
        Float extract_callset_gib = read_float("")
        Float cost = 3
    }
}
