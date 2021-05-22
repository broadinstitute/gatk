version 1.0

import "GvsPrepareCallset.wdl" as GvsPrepareCallset
import "GvsExtractCallset.wdl" as GvsExtractCallset

# Workflow used by AoU to extract variants for a given cohort of sample_names

workflow GvsExtractCohortFromSampleNames {

  input {
    File cohort_sample_names
    String query_project
    String fq_gvs_dataset
    String gvs_extraction_cohorts_dataset
    String gvs_extraction_destination_dataset
    String gvs_extraction_temp_tables_dataset
    String extraction_uuid
    String? output_gcs_dir

    # Extract parameters
    File wgs_intervals
    Int scatter_count

    File reference
    File reference_index
    File reference_dict

    String output_file_base_name

    Boolean do_not_filter_override = false
    String? filter_set_name
    String fq_filter_set_info_table = "~{fq_gvs_dataset}.filter_set_info"
    String fq_filter_set_site_table = "~{fq_gvs_dataset}.filter_set_sites"
    String fq_filter_set_tranches_table = "~{fq_gvs_dataset}.filter_set_tranches"

    # if these are unset, default sensitivity levels will be used
    Float? snps_truth_sensitivity_filter_level_override
    Float? indels_truth_sensitivity_filter_level_override

    File? gatk_override
  }

  call CreateCohortSampleTable {
    input:
      cohort_sample_names = cohort_sample_names,
      fq_gvs_dataset = fq_gvs_dataset,
      gvs_extraction_cohorts_dataset = gvs_extraction_cohorts_dataset,
      query_project = query_project,
      extraction_uuid = extraction_uuid
  }

  call GvsPrepareCallset.GvsPrepareCallset {
    input:
      destination_cohort_table_name   = "destination_~{extraction_uuid}",
      data_project                    = query_project,
      default_dataset                 = "", # unused if fq_* args are given
      fq_petvet_dataset               = fq_gvs_dataset,
      fq_cohort_sample_table          = CreateCohortSampleTable.fq_cohort_sample_table,
      fq_sample_mapping_table         = "~{fq_gvs_dataset}.sample_info",
      fq_temp_table_dataset           = gvs_extraction_temp_tables_dataset,
      fq_destination_dataset          = gvs_extraction_destination_dataset
  }


  call GvsExtractCallset.GvsExtractCallset {
    input:
      data_project = "", # unused if fq_filter_set_* args are given or filtering is off
      query_project = query_project,
      default_dataset = "", # unused if fq_filter_set_* args are given or filtering is off

      wgs_intervals = wgs_intervals,
      scatter_count = scatter_count,

      reference = reference,
      reference_index = reference_index,
      reference_dict = reference_dict,

      fq_samples_to_extract_table = CreateCohortSampleTable.fq_cohort_sample_table,
      fq_cohort_extract_table = GvsPrepareCallset.fq_cohort_extract_table,

      do_not_filter_override = do_not_filter_override,
      filter_set_name = filter_set_name,
      fq_filter_set_info_table =  fq_filter_set_info_table,
      fq_filter_set_site_table =  fq_filter_set_site_table,
      fq_filter_set_tranches_table =  fq_filter_set_tranches_table,
      snps_truth_sensitivity_filter_level_override = snps_truth_sensitivity_filter_level_override,
      indels_truth_sensitivity_filter_level_override = indels_truth_sensitivity_filter_level_override,

      output_file_base_name = output_file_base_name,
      output_gcs_dir = output_gcs_dir,
      gatk_override = gatk_override
  }

  output {
    String fq_cohort_extract_table = GvsPrepareCallset.fq_cohort_extract_table
  }

}

task CreateCohortSampleTable {

  input {
    File cohort_sample_names
    String fq_gvs_dataset
    String gvs_extraction_cohorts_dataset
    String query_project
    String extraction_uuid
  }

  command <<<
    FQ_COHORT_SAMPLE_DATASET=$(echo ~{gvs_extraction_cohorts_dataset} | tr '.' ':')
    FQ_COHORT_SAMPLE_TABLE=${FQ_COHORT_SAMPLE_DATASET}.cohort_~{extraction_uuid}
    FQ_COHORT_SAMPLE_NAME_TABLE_COLON=${FQ_COHORT_SAMPLE_DATASET}.cohort_sample_names_~{extraction_uuid}
    FQ_COHORT_SAMPLE_NAME_TABLE_PERIOD=$(echo ${FQ_COHORT_SAMPLE_NAME_TABLE_COLON} | tr ':' '.')

    bq load --project_id ~{query_project} --format csv ${FQ_COHORT_SAMPLE_NAME_TABLE_COLON} ~{cohort_sample_names} sample_name:STRING

    bq query \
    --project_id ~{query_project} \
    --destination_table ${FQ_COHORT_SAMPLE_TABLE} \
    --use_legacy_sql=false \
    --max_rows=10000000 \
    --allow_large_results \
    "SELECT sample_info.sample_name,
            sample_info.sample_id
     FROM   `~{fq_gvs_dataset}.sample_info` sample_info
            JOIN `${FQ_COHORT_SAMPLE_NAME_TABLE_PERIOD}` cohort_sample_names
              ON sample_info.sample_name = cohort_sample_names.sample_name"

    echo ${FQ_COHORT_SAMPLE_TABLE} | tr ':' '.' > fq_cohort_sample_table.txt
  >>>

  output {
    String fq_cohort_sample_table = read_string("fq_cohort_sample_table.txt")
  }

  runtime {
    memory: "3.75 GB"
    bootDiskSizeGb: "15"
    disks: "local-disk 10 HDD"
    preemptible: 3
    docker: "us.gcr.io/broad-gatk/gatk:4.1.8.0"
  }
}

