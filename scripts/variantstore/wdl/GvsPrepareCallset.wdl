version 1.0

workflow GvsPrepareCallset {
   input {
        String data_project
        String default_dataset
        String destination_cohort_table_name

        # inputs with defaults
        String query_project = data_project
        Array[String]? query_labels
        String destination_project = data_project
        String destination_dataset = default_dataset

        String fq_petvet_dataset = "~{data_project}.~{default_dataset}"
        String fq_cohort_sample_table = "~{data_project}.~{default_dataset}.sample_info"
        String fq_sample_mapping_table = "~{data_project}.~{default_dataset}.sample_info"
        String fq_temp_table_dataset = "~{destination_project}.temp_tables"
        String fq_destination_dataset = "~{destination_project}.~{destination_dataset}"

        Int temp_table_ttl_in_hours = 72

        String? docker
    }

    String docker_final = select_first([docker, "us.gcr.io/broad-dsde-methods/variantstore:ah_var_store_20210507"])

    call PrepareCallsetTask {
        input:
            destination_cohort_table_name   = destination_cohort_table_name,
            query_project                   = query_project,
            query_labels                    = query_labels,
            fq_petvet_dataset               = fq_petvet_dataset,
            fq_cohort_sample_table          = fq_cohort_sample_table,
            fq_sample_mapping_table         = fq_sample_mapping_table,
            fq_temp_table_dataset           = fq_temp_table_dataset,
            fq_destination_dataset          = fq_destination_dataset,
            temp_table_ttl_in_hours         = temp_table_ttl_in_hours,

            docker                          = docker_final
    }

}

task PrepareCallsetTask {
    # indicates that this task should NOT be call cached
    meta {
       volatile: true
    }

    input {
        String destination_cohort_table_name
        String query_project
        Array[String]? query_labels

        String fq_petvet_dataset
        String fq_cohort_sample_table
        String fq_sample_mapping_table
        String fq_temp_table_dataset
        String fq_destination_dataset
        Int temp_table_ttl_in_hours

        File? service_account_json
        String docker
    }

    command <<<
        set -e

        python3 /app/create_cohort_extract_data_table.py \
            --fq_petvet_dataset ~{fq_petvet_dataset} \
            --fq_temp_table_dataset ~{fq_temp_table_dataset} \
            --fq_destination_dataset ~{fq_destination_dataset} \
            --destination_table ~{destination_cohort_table_name} \
            --fq_cohort_sample_names ~{fq_cohort_sample_table} \
            --query_project ~{query_project} \
            ~{"--query_labels " + query_labels} \
            --fq_sample_mapping_table ~{fq_sample_mapping_table} \
            --ttl ~{temp_table_ttl_in_hours} \
            ~{"--sa_key_path " + service_account_json}
    >>>

    runtime {
        docker: docker
        memory: "3 GB"
        disks: "local-disk 100 HDD"
        bootDiskSizeGb: 15
        preemptible: 0
        cpu: 1
    }

 }



