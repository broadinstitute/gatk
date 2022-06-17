version 1.0

workflow GvsPrepareCallset {
  input {
    Boolean go = true
    String project_id
    String dataset_name

    # true for control samples only, false for participant samples only
    Boolean control_samples = false

    String call_set_identifier

    String extract_table_prefix = call_set_identifier
    String query_project = project_id
    String destination_project = project_id
    String destination_dataset = dataset_name
    String fq_temp_table_dataset = "~{destination_project}.~{destination_dataset}" # This only gets used if we have a file of sample names

    Array[String]? query_labels
    File? sample_names_to_extract
    String? service_account_json_path
  }

  String full_extract_prefix = if (control_samples) then "~{extract_table_prefix}_controls" else extract_table_prefix
  String fq_refvet_dataset = "~{project_id}.~{dataset_name}"
  String fq_sample_mapping_table = "~{project_id}.~{dataset_name}.sample_info"
  String fq_destination_dataset = "~{destination_project}.~{destination_dataset}"

  call PrepareRangesCallsetTask {
    input:
      call_set_identifier              = call_set_identifier,
      destination_cohort_table_prefix = full_extract_prefix,
      sample_names_to_extract         = sample_names_to_extract,
      query_project                   = query_project,
      query_labels                    = query_labels,
      fq_refvet_dataset               = fq_refvet_dataset,
      fq_sample_mapping_table         = fq_sample_mapping_table,
      fq_temp_table_dataset           = fq_temp_table_dataset,
      fq_destination_dataset          = fq_destination_dataset,
      temp_table_ttl_in_hours         = 72,
      control_samples                 = control_samples,
      service_account_json_path       = service_account_json_path,
  }

  output {
    String fq_cohort_extract_table_prefix = PrepareRangesCallsetTask.fq_cohort_extract_table_prefix
    Boolean done = true
  }
}

task PrepareRangesCallsetTask {
  input {
    String call_set_identifier
    String destination_cohort_table_prefix
    File? sample_names_to_extract
    String query_project

    Boolean control_samples
    String fq_refvet_dataset
    String fq_sample_mapping_table
    String fq_temp_table_dataset
    String fq_destination_dataset
    Array[String]? query_labels
    Int temp_table_ttl_in_hours = 24

    String? service_account_json_path
  }
  meta {
    # All kinds of BQ reading happening in the referenced Python script.
    volatile: true
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

      python3 /app/create_ranges_cohort_extract_data_table.py \
          --call_set_identifier ~{call_set_identifier} \
          --control_samples ~{control_samples} \
          --fq_ranges_dataset ~{fq_refvet_dataset} \
          --fq_temp_table_dataset ~{fq_temp_table_dataset} \
          --fq_destination_dataset ~{fq_destination_dataset} \
          --destination_cohort_table_prefix ~{destination_cohort_table_prefix} \
          ~{sample_list_param} \
          --query_project ~{query_project} \
          ~{sep=" " query_label_args} \
          --fq_sample_mapping_table ~{fq_sample_mapping_table} \
          --ttl ~{temp_table_ttl_in_hours} \
          $SERVICE_ACCOUNT_STANZA
  >>>
  output {
    String fq_cohort_extract_table_prefix = "~{fq_destination_dataset}.~{destination_cohort_table_prefix}" # implementation detail of create_ranges_cohort_extract_data_table.py
  }

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/variantstore:rsa_metadata_from_python_20220617"
    memory: "3 GB"
    disks: "local-disk 100 HDD"
    bootDiskSizeGb: 15
    preemptible: 0
    cpu: 1
  }
}
