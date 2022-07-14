version 1.0

import "GvsCallsetCost.wdl" as GvsCallsetCost

workflow GvsJointVariantCallsetCost {
    input {
        String project_id
        String dataset_name
        String workspace_namespace
        String workspace_name
        String call_set_identifier
    }

    # the call_set_identifier string is used to name many different things throughout this workflow (BQ tables, vcfs etc),
    # and so make sure nothing is broken by creative users, we replace spaces and underscores with hyphens

    if(false) {
      Array[String] excluded_submission_ids = []
    }
    File gatk_override = "gs://gvs_quickstart_storage/jars/gatk-package-4.2.0.0-552-g0f9780a-SNAPSHOT-local.jar"

    ## TB : GB; 1 : 1024

    Float write_API_cost = 0.025 ## BigQuery Storage Write API: $0.025 per 1 GB. The first 2 TB per month are free.
    Float query_cost = 0.0048828125 ## Queries (on-demand): $5 per TB. The first 1 TB per month is free.
    Float storage_api_cost = 0.00107421875 ## Storage API GB Scanned / Streaming reads (BigQuery Storage Read API):  $1.1 per TB read. Customers can read up to 300 TB of data per month at no charge.


    call GvsCallsetCost.GvsCallsetCost {
        input:
            call_set_identifier = call_set_identifier,
            dataset_name = dataset_name,
            project_id = project_id,
            workspace_namespace = workspace_namespace,
            workspace_name = workspace_name,
            excluded_submission_ids = excluded_submission_ids,
    }

    call FullCosts {
        input:
            vet_size = GvsCallsetCost.vet_gib,
            ref_size = GvsCallsetCost.ref_ranges_gib,
            cost_observability_json = GvsCallsetCost.cost_observability,
            workflow_compute_costs = GvsCallsetCost.workflow_compute_costs,
            storage_api_cost = storage_api_cost,
            query_cost = query_cost,
            write_API_cost = write_API_cost,
    }

        output {
            Float total_cost = FullCosts.total_cost
        }

  }

    task FullCosts {
        meta {
            description: "Calculate costs"
            volatile: true
        }
        input {
            Float vet_size
            Float ref_size
            File cost_observability_json
            File workflow_compute_costs
            Float storage_api_cost
            Float query_cost
            Float write_API_cost
        }

        ## Assign Ids
        # compute cost only

        ## Import Genomes
        # compute cost plus
        Float import_genomes_cost = write_API_cost * (vet_size + ref_size)

        ## AltAllele
        # compute cost plus
        # query_cost

        ## Filter
        # compute cost plus
        # query_cost
        # storage_api_cost

        ## Prepare
        # compute cost plus
        # query_cost

        ## Extract
        # compute cost plus
        # storage_api_cost

            command <<<
                set -e
                cat ~{cost_observability_json} | jq '.| map(select(.event_key=="BigQuery Query Scanned").sum_event_gibibytes | tonumber * 0.0048828125) |add'
                cat ~{cost_observability_json} | jq '.| map(select(.event_key=="Storage API Scanned").sum_event_gibibytes | tonumber * 0.00107421875) |add'
                echo ~{import_genomes_cost}

                cat ~{cost_observability_json} | jq '.| map((select(.event_key=="BigQuery Query Scanned").sum_event_gibibytes | tonumber * 0.0048828125), (select(.event_key=="Storage API Scanned").sum_event_gibibytes | tonumber * 0.00107421875)) | add '

                cat ~{cost_observability_json} | jq '.| map( \
                  (select(.event_key=="BigQuery Query Scanned").sum_event_gibibytes | tonumber * ~{query_cost}), \
                  (select(.event_key=="Storage API Scanned").sum_event_gibibytes | tonumber * ~{storage_api_cost}) \
                ) |add + ~{import_genomes_cost} ' > total_cost.txt

            >>>

        runtime {
            docker: "stedolan/jq:latest"
        }

        output {
            Float total_cost = read_float("total_cost.txt")
        }
    }

