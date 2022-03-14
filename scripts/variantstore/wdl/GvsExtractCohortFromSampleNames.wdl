version 1.0

import "GvsPrepareRangesCallset.wdl" as GvsPrepareCallset
import "GvsExtractCallset.wdl" as GvsExtractCallset

# Workflow used by AoU to extract variants for a given cohort of sample_names

workflow GvsExtractCohortFromSampleNames {

  input {
    File cohort_sample_names
    String query_project
    String gvs_project
    String gvs_dataset
    String cohort_table_prefix

    # not using the defaults in GvsPrepareCallset because we're using pre created datasets defined by the caller
    String destination_dataset_name
    String destination_project_id
    String extraction_uuid
    String filter_set_name
    String fq_gvs_extraction_temp_tables_dataset
    String output_file_base_name
    Int scatter_count

    String? output_gcs_dir
    String? service_account_json_path

    Int? extract_preemptible_override
    Int? extract_maxretries_override
    Int? split_intervals_disk_size_override
    Int? split_intervals_mem_override

    File? gatk_override
  }

  call GvsPrepareCallset.GvsPrepareCallset {
    input:
      extract_table_prefix            = cohort_table_prefix,
      sample_names_to_extract         = cohort_sample_names,
      project_id                      = gvs_project,
      query_labels                    = ["extraction_uuid=~{extraction_uuid}"],
      query_project                   = query_project,
      dataset_name                    = gvs_dataset, # unused if fq_* args are given
      destination_project             = destination_project_id,
      destination_dataset             = destination_dataset_name,
      fq_temp_table_dataset           = fq_gvs_extraction_temp_tables_dataset,
      service_account_json_path       = service_account_json_path
  }

  call GvsExtractCallset.GvsExtractCallset {
    input:
      project_id = gvs_project,
      query_project = query_project,
      dataset_name = gvs_dataset,
      extract_table_prefix = cohort_table_prefix,

      scatter_count = scatter_count,
      filter_set_name = filter_set_name,
      output_file_base_name = output_file_base_name,
      output_gcs_dir = output_gcs_dir,
      service_account_json_path = service_account_json_path,

      extract_preemptible_override = extract_preemptible_override,
      extract_maxretries_override = extract_maxretries_override,
      split_intervals_disk_size_override = split_intervals_disk_size_override,
      split_intervals_mem_override = split_intervals_mem_override
  }

  output {
    Float total_vcfs_size_mb = GvsExtractCallset.total_vcfs_size_mb
  }

}
