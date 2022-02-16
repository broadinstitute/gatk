version 1.0

import "GvsExtractCallset.wdl" as CallsetInterval
import "https://raw.githubusercontent.com/broadinstitute/gatk/a08a18f709bf85e7bd429d41d6e8a36876f1098f/scripts/variantstore/wdl/GvsUtils.wdl" as Utils

workflow GvsRescatterCallsetInterval {
  input {
    String data_project
    String default_dataset
    String extract_table_prefix
    String output_file_base_name
    File reference
    File reference_dict
    File reference_index
    String? final_output_gcs_dir
    String? filter_set_name

    Int re_scatter_count
    String interval_file_dir
    Array[String] intervals_to_scatter    #e.g. ["0001", "0413", "9839"] match format of the interval file names

    Int? extract_preemptible_override
    Int? merge_disk_override
    File? gatk_override
    String? service_account_json_path
  }

  scatter(i in range(length(intervals_to_scatter))) {
    # take out leading 0s from interval file name number for VCF and index
    Int shard_num = intervals_to_scatter[i]
    String vcf_basename = "${output_file_base_name}_${shard_num}"

    call CallsetInterval.GvsExtractCallset as ExtractInterval {
      input:
        data_project = data_project,
        default_dataset = default_dataset,
        extract_table_prefix = extract_table_prefix,
        output_file_base_name = vcf_basename,
        reference = reference,
        reference_dict = reference_dict,
        reference_index = reference_index,
        scatter_count = re_scatter_count,
        wgs_intervals = sub(interval_file_dir, "/$", "") + '/' + intervals_to_scatter[i] + "-scattered.interval_list",
        extract_preemptible_override = extract_preemptible_override,
        filter_set_name = filter_set_name,
        gatk_override = gatk_override,
        service_account_json_path = service_account_json_path
    }

    call Utils.MergeVCFs as MergeVCFs {
      input:
        input_vcfs = ExtractInterval.output_vcfs,
        output_vcf_name = "${vcf_basename}.vcf.gz",
        output_directory = final_output_gcs_dir,
        merge_disk_override = merge_disk_override,
        service_account_json_path = service_account_json_path
    }
  }

  output {
    Array[File] output_vcf = MergeVCFs.output_vcf
    Array[File] output_vcf_index = MergeVCFs.output_vcf_index
  }
}
