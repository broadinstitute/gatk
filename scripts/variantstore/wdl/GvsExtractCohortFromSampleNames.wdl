version 1.0

import "GvsPrepareCallset.wdl"
import "GvsExtractCallset.wdl" as GvsExtractCallset

# Workflow used by AoU to extract variants for a given cohort of sample_names

workflow GvsExtractCohortFromSampleNames {

  input {
    File cohort_sample_names
    String query_project
    String gvs_project
    String gvs_dataset

    # not using the defaults in GvsPrepareCallset because we're using pre created datasets defined by the caller
    String fq_gvs_extraction_destination_dataset
    String fq_gvs_extraction_temp_tables_dataset

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
    String fq_filter_set_info_table = "~{gvs_project}.~{gvs_dataset}.filter_set_info"
    String fq_filter_set_site_table = "~{gvs_project}.~{gvs_dataset}.filter_set_sites"
    String fq_filter_set_tranches_table = "~{gvs_project}.~{gvs_dataset}.filter_set_tranches"

    # if these are unset, default sensitivity levels will be used
    Float? snps_truth_sensitivity_filter_level_override
    Float? indels_truth_sensitivity_filter_level_override

    File? gatk_override
  }

  call GvsPrepareCallset.GvsPrepareCallset {
    input:
      destination_cohort_table_prefix = extraction_uuid,
      sample_names_to_extract         = cohort_sample_names,
      data_project                    = query_project,
      default_dataset                 = gvs_dataset, # unused if fq_* args are given
      fq_petvet_dataset               = "~{gvs_project}.~{gvs_dataset}",
      fq_sample_mapping_table         = "~{gvs_project}.~{gvs_dataset}.sample_info",
      fq_temp_table_dataset           = fq_gvs_extraction_temp_tables_dataset,
      fq_destination_dataset          = fq_gvs_extraction_destination_dataset
  }

  call GvsExtractCallset.GvsExtractCallset {
    input:
      data_project = gvs_project,
      query_project = query_project,
      default_dataset = gvs_dataset,

      wgs_intervals = wgs_intervals,
      scatter_count = scatter_count,

      reference = reference,
      reference_index = reference_index,
      reference_dict = reference_dict,

      fq_cohort_extract_table_prefix = GvsPrepareCallset.fq_cohort_extract_table_prefix,

      do_not_filter_override = do_not_filter_override,
      filter_set_name = filter_set_name,
      snps_truth_sensitivity_filter_level_override = snps_truth_sensitivity_filter_level_override,
      indels_truth_sensitivity_filter_level_override = indels_truth_sensitivity_filter_level_override,

      output_file_base_name = output_file_base_name,
      output_gcs_dir = output_gcs_dir,
      gatk_override = gatk_override
  }

}
