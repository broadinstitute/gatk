version 1.0

import "GvsPrepareCallset.wdl" as GvsPrepareCallset
import "GvsExtractCallset.wdl" as GvsExtractCallset

# Workflow used by AoU to extract variants for a given cohort of sample_names

workflow GvsExtractCohortFromSampleNames {

  input {
    File cohort_sample_names
    String query_project
    String gvs_dataset
    String gvs_extraction_cohorts_dataset
    String gvs_extraction_destination_dataset
    String gvs_extraction_temp_tables_dataset
    String extraction_uuid
    String output_gcs_dir

    # Extract parameters
    File wgs_intervals
    Int scatter_count

    File reference
    File reference_index
    File reference_dict

    String output_file_base_name

    String? fq_filter_set_table
    String? filter_set_name
    File? gatk_override
  }

  call CreateCohortSampleTable {
    input:
      cohort_sample_names = cohort_sample_names,
      gvs_dataset = gvs_dataset,
      gvs_extraction_cohorts_dataset = gvs_extraction_cohorts_dataset,
      query_project = query_project,
      table_name = extraction_uuid
  }

  call GvsPrepareCallset.GvsPrepareCallset {
    input:
      destination_cohort_table_name   = extraction_uuid,
      data_project                    = query_project,
      default_dataset                 = "", # unused if fq_* args are given
      fq_petvet_dataset               = gvs_dataset,
      fq_cohort_sample_table          = CreateCohortSampleTable.fq_cohort_sample_table,
      fq_sample_mapping_table         = "~{gvs_dataset}.sample_info",
      fq_temp_table_dataset           = gvs_extraction_temp_tables_dataset,
      fq_destination_dataset          = gvs_extraction_destination_dataset
  }


  call GvsExtractCallset.GvsExtractCallset {
    input:
      data_project = "", # unused if fq_filter_set_* args are given or filtering is off
      query_project = query_project,
      default_dataset = "", # unused if fq_filter_set_* args are given or filtering is off
      filter_set_name = "", # unused if fq_filter_set_* args are given or filtering is off

      wgs_intervals = wgs_intervals,
      scatter_count = scatter_count,

      reference = reference,
      reference_index = reference_index,
      reference_dict = reference_dict,

      fq_samples_to_extract_table = CreateCohortSampleTable.fq_cohort_sample_table,
      fq_cohort_extract_table = GvsPrepareCallset.fq_cohort_extract_table,
      fq_filter_set_info_table =  "", # unused if filtering is off
      fq_filter_set_site_table =  "", # unused if filtering is off
      fq_filter_set_tranches_table =  "", # unused if filtering is off
      do_not_filter_override = true,

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
    String gvs_dataset
    String gvs_extraction_cohorts_dataset
    String query_project
    String table_name
  }

  command <<<
    echo "SELECT
    sample_id,
    sample_name
    FROM
    \`~{gvs_dataset}.sample_info\`
    WHERE
    sample_name IN " > create_cohort.sql

    PARTICIPANT_IDS=$(cat ~{cohort_sample_names} | awk '{print "\""$0"\""}' | paste -sd ",")
    echo "($PARTICIPANT_IDS)" >> create_cohort.sql

    DESTINATION_DATASET=$(echo ~{gvs_extraction_cohorts_dataset} | tr '.' ':')

    bq query \
    --project_id ~{query_project} \
    --destination_table ${DESTINATION_DATASET}.~{table_name} \
    --use_legacy_sql=false \
    --max_rows=10000000 \
    --allow_large_results < create_cohort.sql

    echo "~{gvs_extraction_cohorts_dataset}.~{table_name}" > fq_cohort_sample_table.txt
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

