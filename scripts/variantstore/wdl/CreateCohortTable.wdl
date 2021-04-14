version 1.0

workflow CreateCohortTable {
   input {
        String data_project
        String dataset
        String fq_cohort_sample_table = "${data_project}.${dataset}.sample_info"
        String fq_sample_mapping_table = "${data_project}.${dataset}.sample_info"

        String destination_cohort_table_name

        # TODO testing remove this ?
        String? query_project = data_project
        String destination_project = data_project
        String destination_dataset = dataset

        String? docker
    }

    # TODO update this docker source
    String docker_final = select_first([docker, "us.gcr.io/broad-dsde-methods/variantstore:latest"])

    call CreateCohortTableTask {
        input:
            data_project                    = data_project,
            dataset                         = dataset,
            fq_cohort_sample_table          = fq_cohort_sample_table,
            fq_sample_mapping_table         = fq_sample_mapping_table,

            destination_cohort_table_name   = destination_cohort_table_name,

            query_project                   = query_project,
            destination_project             = destination_project,
            destination_dataset             = destination_dataset,

            docker                          = docker_final
    }

}

task CreateCohortTableTask {
    # indicates that this task should NOT be call cached
    meta {
       volatile: true
    }

    input {
        String data_project
        String dataset
        String fq_cohort_sample_table
        String fq_sample_mapping_table

        String destination_cohort_table_name

        # TODO testing remove this ?
        String? query_project
        String destination_project
        String destination_dataset

        File? service_account_json
        String docker
    }

    command <<<
        set -e

        python3 /app/create_cohort_data_table.py \
            --fq_petvet_dataset ~{data_project}.~{dataset} \
            --fq_temp_table_dataset ~{destination_project}.temp_tables \
            --fq_destination_dataset ~{destination_project}.~{destination_dataset} \
            --destination_table ~{destination_cohort_table_name} \
            --fq_cohort_sample_names ~{fq_cohort_sample_table} \
            --query_project ~{query_project} \
            --fq_sample_mapping_table ~{fq_sample_mapping_table} \
            ~{"--sa_key_path " + service_account_json}
    >>>

    runtime {
        docker: docker
        memory: "10 GB"
        disks: "local-disk 100 HDD"
        bootDiskSizeGb: 15
        preemptible: 0
        cpu: 1
    }

 }



