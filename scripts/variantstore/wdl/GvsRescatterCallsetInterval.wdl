version 1.0

import "GvsExtractCallset.wdl" as CallsetInterval
import "GvsUtils.wdl" as Utils

workflow GvsRescatterCallsetInterval {
  input {
    String git_branch_or_tag
    String dataset_name
    String extract_table_prefix
    String filter_set_name
    String interval_file_dir
    Array[String] intervals_to_scatter    #e.g. ["0001", "0413", "9839"] match format of the interval file names
    String output_file_base_name
    String project_id
    Int re_scatter_count

    Int? extract_preemptible_override
    String? final_output_gcs_dir
    File? gatk_override
    Int? merge_disk_override
    String? variants_docker
    String? cloud_sdk_docker
    String? gatk_docker
  }

  # Always call `GetToolVersions` to get the git hash for this run as this is a top-level-only WDL (i.e. there are
  # no calling WDLs that might supply `git_hash`).
  call Utils.GetToolVersions {
    input:
      git_branch_or_tag = git_branch_or_tag,
  }

  String effective_gatk_docker = select_first([gatk_docker, GetToolVersions.gatk_docker])
  String effective_cloud_sdk_docker = select_first([cloud_sdk_docker, GetToolVersions.cloud_sdk_docker])
  String effective_variants_docker = select_first([variants_docker, GetToolVersions.variants_docker])

  scatter(i in range(length(intervals_to_scatter))) {
    # take out leading 0s from interval file name number for VCF and index
    Int shard_num = intervals_to_scatter[i]
    String vcf_basename = "${output_file_base_name}_${shard_num}"

    call CallsetInterval.GvsExtractCallset as ExtractInterval {
      input:
        project_id = project_id,
        dataset_name = dataset_name,
        extract_table_prefix = extract_table_prefix,
        output_file_base_name = vcf_basename,
        scatter_count = re_scatter_count,
        interval_list = sub(interval_file_dir, "/$", "") + '/' + intervals_to_scatter[i] + "-scattered.interval_list",
        extract_preemptible_override = extract_preemptible_override,
        filter_set_name = filter_set_name,
        gatk_docker = effective_gatk_docker,
        cloud_sdk_docker = effective_cloud_sdk_docker,
        variants_docker = effective_variants_docker,
    }

    call Utils.MergeVCFs as MergeVCFs {
      input:
        input_vcfs = ExtractInterval.output_vcfs,
        output_vcf_name = "${vcf_basename}.vcf.gz",
        output_directory = final_output_gcs_dir,
        merge_disk_override = merge_disk_override,
        gatk_docker = effective_gatk_docker,
    }
  }

  output {
    Array[File] output_vcf = MergeVCFs.output_vcf
    Array[File] output_vcf_index = MergeVCFs.output_vcf_index
    String recorded_git_hash = GetToolVersions.git_hash
  }
}
