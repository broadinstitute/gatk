version 1.0

import "GvsUtils.wdl" as Utils

workflow GvsPrepareCallset {
  input {
    Boolean go = true
    String project_id
    String dataset_name

    # true for control samples only, false for participant samples only
    Boolean control_samples = false

    String call_set_identifier

    String extract_table_prefix = call_set_identifier
    Boolean enable_extract_table_ttl = true
    String query_project = project_id
    String destination_project = project_id
    String destination_dataset = dataset_name
    String fq_temp_table_dataset = "~{destination_project}.~{destination_dataset}" # This only gets used if we have a file of sample names

    Array[String]? query_labels
    File? sample_names_to_extract
    Boolean only_output_vet_tables = false
    Boolean write_cost_to_db = true
    String? cloud_sdk_docker
    String? variants_docker
    String? git_branch_or_tag
    String? git_hash
  }

  String full_extract_prefix = if (control_samples) then "~{extract_table_prefix}_controls" else extract_table_prefix
  String fq_refvet_dataset = "~{project_id}.~{dataset_name}"
  String fq_sample_mapping_table = "~{project_id}.~{dataset_name}.sample_info"
  String fq_destination_dataset = "~{destination_project}.~{destination_dataset}"

  if (!defined(git_hash) || !defined(variants_docker) || !defined(cloud_sdk_docker)) {
    call Utils.GetToolVersions {
      input:
        git_branch_or_tag = git_branch_or_tag,
    }
  }

  String effective_variants_docker = select_first([variants_docker, GetToolVersions.variants_docker])
  String effective_cloud_sdk_docker = select_first([cloud_sdk_docker, GetToolVersions.cloud_sdk_docker])
  String effective_git_hash = select_first([git_hash, GetToolVersions.git_hash])

  call Utils.GetBQTableLastModifiedDatetime as RefTableDatetimeCheck {
    input:
      project_id = project_id,
      fq_table = "~{project_id}.~{dataset_name}.ref_ranges_001",
      cloud_sdk_docker = effective_cloud_sdk_docker,
  }

  call Utils.IsUsingCompressedReferences {
    input:
      project_id = query_project,
      dataset_name = dataset_name,
      ref_table_timestamp = RefTableDatetimeCheck.last_modified_timestamp,
      cloud_sdk_docker = effective_cloud_sdk_docker,
  }

  call PrepareRangesCallsetTask {
    input:
      call_set_identifier             = call_set_identifier,
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
      only_output_vet_tables          = only_output_vet_tables,
      write_cost_to_db                = write_cost_to_db,
      variants_docker                 = effective_variants_docker,
      use_compressed_references       = IsUsingCompressedReferences.is_using_compressed_references,
      enable_extract_table_ttl        = enable_extract_table_ttl
  }

  output {
    String fq_cohort_extract_table_prefix = PrepareRangesCallsetTask.fq_cohort_extract_table_prefix
    String recorded_git_hash = effective_git_hash
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
    Boolean only_output_vet_tables
    Boolean write_cost_to_db
    Boolean use_compressed_references
    String variants_docker
    Boolean enable_extract_table_ttl
  }
  meta {
    # All kinds of BQ reading happening in the referenced Python script.
    volatile: true
  }
  # Note the coercion of optional query_labels using select_first([expr, default])
  Array[String] query_label_args = if defined(query_labels) then prefix("--query_labels ", select_first([query_labels])) else []

  String use_sample_names_file = if (defined(sample_names_to_extract)) then 'true' else 'false'
  String sample_list_param = if (defined(sample_names_to_extract)) then '--sample_names_to_extract sample_names_file' else '--fq_cohort_sample_names ' + fq_sample_mapping_table

  parameter_meta {
    sample_names_to_extract: {
      localization_optional: true
    }
  }

  command <<<
    # Prepend date, time and pwd to xtrace log entries.
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o pipefail -o xtrace

      echo ~{sample_list_param}

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
          ~{true="--only_output_vet_tables True" false='' only_output_vet_tables} \
          ~{true="--write_cost_to_db True" false="--write_cost_to_db ''" write_cost_to_db} \
          ~{true="--use_compressed_references True" false='' use_compressed_references} \
          ~{true="--enable_extract_table_ttl True" false='' enable_extract_table_ttl}

  >>>
  output {
    String fq_cohort_extract_table_prefix = "~{fq_destination_dataset}.~{destination_cohort_table_prefix}" # implementation detail of create_ranges_cohort_extract_data_table.py
  }

  runtime {
    docker: variants_docker
    memory: "3 GB"
    disks: "local-disk 100 HDD"
    bootDiskSizeGb: 15
    preemptible: 0
    cpu: 1
  }
}
