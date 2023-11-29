version 1.0

import "GvsPrepareRangesCallset.wdl" as GvsPrepareCallset
import "GvsExtractCallset.wdl" as GvsExtractCallset
import "GvsUtils.wdl" as Utils

# Workflow used by AoU to extract variants for a given cohort of sample_names

workflow GvsExtractCohortFromSampleNames {

  input {
    # cohort_sample_names_array will take precedence over cohort_sample_names if both are set
    Array[String]? cohort_sample_names_array
    File? cohort_sample_names

    String query_project
    String gvs_project
    String gvs_dataset
    String call_set_identifier
    String cohort_table_prefix

    # not using the defaults in GvsPrepareCallset because we're using pre created datasets defined by the caller
    String destination_dataset_name
    String destination_project_id
    String? fq_gvs_extraction_temp_tables_dataset
    String extraction_uuid
    String filter_set_name
    String output_file_base_name

    String? output_gcs_dir
    # set to "NONE" if all the reference data was loaded into GVS in GvsImportGenomes
    String drop_state = "NONE"

    Int? extract_preemptible_override
    Int? extract_maxretries_override
    Int? extract_scatter_count_override
    Int? extract_memory_override
    Int? extract_disk_override
    Int? split_intervals_disk_size_override
    Int? split_intervals_mem_override

    String? git_branch_or_tag
    File? gatk_override
    String? cloud_sdk_docker
  }

  Boolean write_cost_to_db = if ((gvs_project != destination_project_id) || (gvs_project != query_project)) then false else true

  # Always call `GetToolVersions` to get the git hash for this run as this is a top-level-only WDL (i.e. there are
  # no calling WDLs that might supply `git_hash`).
  call Utils.GetToolVersions {
    input:
      git_branch_or_tag = git_branch_or_tag,
  }

  String effective_cloud_sdk_docker = select_first([cloud_sdk_docker, GetToolVersions.cloud_sdk_docker])

  call Utils.GetBQTableLastModifiedDatetime as SamplesTableDatetimeCheck {
    input:
      project_id = query_project,
      fq_table = "~{gvs_project}.~{gvs_dataset}.sample_info",
      cloud_sdk_docker = effective_cloud_sdk_docker
  }

  call Utils.GetNumSamplesLoaded {
    input:
      fq_sample_table = "~{gvs_project}.~{gvs_dataset}.sample_info",
      project_id = gvs_project,
      sample_table_timestamp = SamplesTableDatetimeCheck.last_modified_timestamp,
      cloud_sdk_docker = effective_cloud_sdk_docker
  }

  Int effective_scatter_count = if defined(extract_scatter_count_override) then select_first([extract_scatter_count_override])
                                else if GetNumSamplesLoaded.num_samples < 100 then 50 # Quickstart
                                     else if GetNumSamplesLoaded.num_samples < 1000 then 250
                                          else if GetNumSamplesLoaded.num_samples < 5000 then 500
                                               else if GetNumSamplesLoaded.num_samples < 20000 then 1000 # Stroke Anderson
                                                    else if GetNumSamplesLoaded.num_samples < 50000 then 5000
                                                         else if GetNumSamplesLoaded.num_samples < 100000 then 10000 # Charlie
                                                              else 20000

  # writing the array to a file has to be done in a task
  # https://support.terra.bio/hc/en-us/community/posts/360071465631-write-lines-write-map-write-tsv-write-json-fail-when-run-in-a-workflow-rather-than-in-a-task
  if (defined(cohort_sample_names_array)) {
    call write_array_task {
      input:
        input_array = select_first([cohort_sample_names_array]),
        cloud_sdk_docker = effective_cloud_sdk_docker,
    }
  }

  File cohort_sample_names_file = select_first([write_array_task.output_file, cohort_sample_names])

  call GvsPrepareCallset.GvsPrepareCallset {
    input:
      call_set_identifier             = cohort_table_prefix,
      extract_table_prefix            = cohort_table_prefix,
      sample_names_to_extract         = cohort_sample_names_file,
      project_id                      = gvs_project,
      query_labels                    = ["extraction_uuid=~{extraction_uuid}"],
      query_project                   = query_project,
      dataset_name                    = gvs_dataset, # unused if fq_* args are given
      destination_project             = destination_project_id,
      destination_dataset             = destination_dataset_name,
      fq_temp_table_dataset           = fq_gvs_extraction_temp_tables_dataset,
      write_cost_to_db                = write_cost_to_db,
      cloud_sdk_docker                = effective_cloud_sdk_docker,
      extract_table_ttl               = false,
  }

  call GvsExtractCallset.GvsExtractCallset {
    input:
      go = GvsPrepareCallset.done,
      project_id = gvs_project,
      query_project = query_project,
      dataset_name = gvs_dataset,
      call_set_identifier = call_set_identifier,
      cohort_project_id = destination_project_id,
      cohort_dataset_name = destination_dataset_name,
      extract_table_prefix = cohort_table_prefix,

      scatter_count = effective_scatter_count,
      filter_set_name = filter_set_name,
      output_file_base_name = output_file_base_name,
      output_gcs_dir = output_gcs_dir,

      drop_state = drop_state,
      extract_preemptible_override = extract_preemptible_override,
      extract_maxretries_override = extract_maxretries_override,
      split_intervals_disk_size_override = split_intervals_disk_size_override,
      split_intervals_mem_override = split_intervals_mem_override,
      memory_override = extract_memory_override,
      disk_override = extract_disk_override,

      gatk_override = gatk_override,
      write_cost_to_db = write_cost_to_db
  }

  output {
    Float total_vcfs_size_mb = GvsExtractCallset.total_vcfs_size_mb
    String recorded_git_hash = GetToolVersions.git_hash
  }

}

task write_array_task {
  input {
    Array[String] input_array
    String cloud_sdk_docker
  }

  command <<<
  >>>

  output {
    File output_file = write_lines(input_array)
  }

  runtime {
    docker: cloud_sdk_docker
  }
}

