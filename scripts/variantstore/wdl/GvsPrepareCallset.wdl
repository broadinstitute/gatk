version 1.0

workflow GvsPrepareCallset {
   input {
        String data_project
        String default_dataset
        String destination_cohort_table_prefix
        File? sample_names_to_extract

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
        Boolean skip_pet_insert = false
        String? service_account_json_path
        String? docker
    }

    String docker_final = select_first([docker, "us.gcr.io/broad-dsde-methods/variantstore:ah_var_store_20211014_2"])

    call PrepareCallsetTask {
        input:
            destination_cohort_table_prefix = destination_cohort_table_prefix,
            sample_names_to_extract         = sample_names_to_extract,
            query_project                   = query_project,
            query_labels                    = query_labels,
            fq_petvet_dataset               = fq_petvet_dataset,
            fq_sample_mapping_table         = fq_sample_mapping_table,
            fq_temp_table_dataset           = fq_temp_table_dataset,
            fq_destination_dataset          = fq_destination_dataset,
            temp_table_ttl_in_hours         = temp_table_ttl_in_hours,
            skip_pet_insert                 = skip_pet_insert,
            service_account_json_path       = service_account_json_path,
            docker                          = docker_final
    }

    output {
      String fq_cohort_extract_table_prefix = PrepareCallsetTask.fq_cohort_extract_table_prefix
    }

}

task PrepareCallsetTask {
    # indicates that this task should NOT be call cached
    meta {
       volatile: true
    }

    input {
        String destination_cohort_table_prefix
        File? sample_names_to_extract
        String query_project
        Array[String]? query_labels

        String fq_petvet_dataset
        String fq_sample_mapping_table
        String fq_temp_table_dataset
        String fq_destination_dataset
        Int temp_table_ttl_in_hours
        Boolean skip_pet_insert

        String? service_account_json_path
        String docker
    }
    # Note the coercion of optional query_labels using select_first([expr, default])
    Array[String] query_label_args = if defined(query_labels) then prefix("--query_labels ", select_first([query_labels])) else []

    String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'
    String use_sample_names_file = if (defined(sample_names_to_extract)) then 'true' else 'false'
    String sample_list_param = if (defined(sample_names_to_extract)) then '--sample_names_to_extract sample_names_file' else '--fq_cohort_sample_names ' + fq_sample_mapping_table

    parameter_meta {
      sample_names_to_extract: {
        localization_optional: true
      }
    }

    command <<<
        set -e

        echo ~{sample_list_param}

        if [ ~{has_service_account_file} = 'true' ]; then
            gsutil cp ~{service_account_json_path} local.service_account.json
            SERVICE_ACCOUNT_STANZA="--sa_key_path local.service_account.json "
        fi

        if [ ~{use_sample_names_file} = 'true' ]; then
            gsutil cp  ~{sample_names_to_extract} sample_names_file
        fi

        python3 /app/create_cohort_extract_data_table.py \
            --fq_petvet_dataset ~{fq_petvet_dataset} \
            --fq_temp_table_dataset ~{fq_temp_table_dataset} \
            --fq_destination_dataset ~{fq_destination_dataset} \
            --destination_cohort_table_prefix ~{destination_cohort_table_prefix} \
            ~{sample_list_param} \
            --query_project ~{query_project} \
            ~{sep=" " query_label_args} \
            --fq_sample_mapping_table ~{fq_sample_mapping_table} \
            --ttl ~{temp_table_ttl_in_hours} \
            --skip_pet_insert ~{skip_pet_insert} \
            $SERVICE_ACCOUNT_STANZA
    >>>

    output {
      String fq_cohort_extract_table_prefix = "~{fq_destination_dataset}.~{destination_cohort_table_prefix}" # implementation detail of create_cohort_extract_data_table.py
    }

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
    String service_account_json_path
  }

  command {
    set -euo pipefail

    gsutil cp ~{service_account_json_path} local.service_account.json
    gcloud auth activate-service-account --key-file=local.service_account.json
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
