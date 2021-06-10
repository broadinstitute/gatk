version 1.0

workflow GvsPrepareCallset {
   input {
        String data_project
        String default_dataset
        String destination_cohort_table_prefix
        File sample_names_to_extract
        Boolean localize_sample_names_with_service_account = false

        # inputs with defaults
        String query_project = data_project
        Array[String]? query_labels
        String destination_project = data_project
        String destination_dataset = default_dataset

        String fq_petvet_dataset = "~{data_project}.~{default_dataset}"
        String fq_sample_mapping_table = "~{data_project}.~{default_dataset}.sample_info"
        String fq_temp_table_dataset = "~{destination_project}.temp_tables"
        String fq_destination_dataset = "~{destination_project}.~{destination_dataset}"

        Int temp_table_ttl_in_hours = 72
        File? service_account_json
        String? docker
    }

    String docker_final = select_first([docker, "us.gcr.io/broad-dsde-methods/variantstore:ah_var_store_20210604"])

    if (localize_sample_names_with_service_account && defined(service_account_json)) {
        call LocalizeFile {
            input:
              file = "~{sample_names_to_extract}",
              service_account_json = select_first([service_account_json])
        }
    }

    call PrepareCallsetTask {
        input:
            destination_cohort_table_prefix = destination_cohort_table_prefix,
            sample_names_to_extract         = select_first([LocalizeFile.localized_file, sample_names_to_extract]),
            query_project                   = query_project,
            query_labels                    = query_labels,
            fq_petvet_dataset               = fq_petvet_dataset,
            fq_sample_mapping_table         = fq_sample_mapping_table,
            fq_temp_table_dataset           = fq_temp_table_dataset,
            fq_destination_dataset          = fq_destination_dataset,
            temp_table_ttl_in_hours         = temp_table_ttl_in_hours,
            service_account_json            = service_account_json,

            docker                          = docker_final
    }

}

task PrepareCallsetTask {
    # indicates that this task should NOT be call cached
    meta {
       volatile: true
    }

    input {
        String destination_cohort_table_prefix
        File sample_names_to_extract
        String query_project
        Array[String]? query_labels

        String fq_petvet_dataset
        String fq_sample_mapping_table
        String fq_temp_table_dataset
        String fq_destination_dataset
        Int temp_table_ttl_in_hours

        File? service_account_json
        String docker
    }
    # Note the coercion of optional query_labels using select_first([expr, default])
    Array[String] query_label_args = if defined(query_labels) then prefix("--query_labels ", select_first([query_labels])) else []

    command <<<
        set -e

        python3 /app/create_cohort_extract_data_table.py \
            --fq_petvet_dataset ~{fq_petvet_dataset} \
            --fq_temp_table_dataset ~{fq_temp_table_dataset} \
            --fq_destination_dataset ~{fq_destination_dataset} \
            --destination_cohort_table_prefix ~{destination_cohort_table_prefix} \
            --sample_names_to_extract ~{sample_names_to_extract} \
            --query_project ~{query_project} \
            ~{sep=" " query_label_args} \
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

task LocalizeFile {
  input {
    String file
    File service_account_json
  }

  command {
    set -euo pipefail

    gcloud auth activate-service-account --key-file='~{service_account_json}'
    gsutil cp '~{file}' .
  }

  output {
    File localized_file = basename(file)
  }

  runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
    memory: "3.75 GiB"
    cpu: "1"
    disks: "local-disk 50 HDD"
  }
}


