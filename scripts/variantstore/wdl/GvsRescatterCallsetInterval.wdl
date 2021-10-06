version 1.0

import "GvsExtractCallset.wdl" as CallsetInterval

workflow GvsRescatterCallsetInterval {
  input {
    String data_project
    String default_dataset
    String fq_cohort_extract_table_prefix
    String output_file_base_name
    File reference
    File reference_dict
    File reference_index
    String final_output_gcs_dir
    String subshards_gcs_directory
    String? filter_set_name

    Int re_scatter_count
    String interval_file_dir
    Array[String] intervals_to_scatter    #e.g. [0001, 0413, 9839]

    Int? extract_preemptible_override
    Int? merge_disk_override
    File? gatk_override
    String? service_account_json_path
  }

  scatter(i in range(length(intervals_to_scatter))) {
    call CallsetInterval.GvsExtractCallset as ExtractInterval {
      input:
        data_project = data_project,
        default_dataset = default_dataset,
        fq_cohort_extract_table_prefix = fq_cohort_extract_table_prefix,
        output_file_base_name = output_file_base_name + '_' + intervals_to_scatter[i],
        reference = reference,
        reference_dict = reference_dict,
        reference_index = reference_index,
        scatter_count = re_scatter_count,
        wgs_intervals = sub(interval_file_dir, "/$", "") + '/' + intervals_to_scatter[i] + "-scattered.interval_list",
        extract_preemptible_override = extract_preemptible_override,
        filter_set_name = filter_set_name,
        gatk_override = gatk_override,
#        output_gcs_dir = subshards_gcs_directory,
        service_account_json_path = service_account_json_path
    }

#    call GenerateOrderedPaths as VCFpaths {
#      input:
#        root_path = sub(subshards_gcs_directory, "/$", "") + '/' + output_file_base_name + '_' + intervals_to_scatter[i] + '_',
#        num_files = re_scatter_count,
#        path_suffix = ".vcf.gz"
#    }

    call MergeVCFs {
      input:
        input_vcfs = ExtractInterval.output_vcfs,
        output_vcf_name = "${output_file_base_name}_${intervals_to_scatter[i]}.vcf.gz",
        output_directory = subshards_gcs_directory, # replace with final_output_gcs_dir after tested
        merge_disk_override = merge_disk_override,
        service_account_json_path = service_account_json_path
    }
  }
}

task MergeVCFs {
  input {
    Array[File] input_vcfs
    String gather_type = "BLOCK"
    String output_vcf_name
    String output_directory
    Int? merge_disk_override
    String? service_account_json_path
    File? gatk_override
  }

  Int disk_size = if (defined(merge_disk_override)) then merge_disk_override else ceil(size(input_vcfs, "GiB") * 2.5) + 10

  parameter_meta {
    input_vcfs: {
      localization_optional: true
    }
  }

  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

  command {
    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

    if [ ~{has_service_account_file} = 'true' ]; then
      gsutil cp ~{service_account_json_path} local.service_account.json
      gcloud auth activate-service-account --key-file=local.service_account.json
      export GOOGLE_APPLICATION_CREDENTIALS=local.service_account.json
    fi

    gatk --java-options -Xmx3g GatherVcfsCloud \
    --ignore-safety-checks --gather-type ~{gather_type} \
    --create-output-variant-index false \
    -I ~{sep=' -I ' input_vcfs} \
    --output ~{output_vcf_name}

    tabix ~{output_vcf_name}

    # Drop trailing slash if one exists
    OUTPUT_GCS_DIR=$(echo ~{output_directory} | sed 's/\/$//')

    gsutil cp ~{output_vcf_name} $OUTPUT_GCS_DIR/
    gsutil cp ~{output_vcf_name}.tbi $OUTPUT_GCS_DIR/
  }

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:varstore_d8a72b825eab2d979c8877448c0ca948fd9b34c7_change_to_hwe"
    preemptible: 1
    memory: "3 GiB"
    disks: "local-disk ~{disk_size} HDD"
  }

  output {
    File output_vcf = "~{output_vcf_name}"
    File output_vcf_index = "~{output_vcf_name}.tbi"
  }
}

task GenerateOrderedPaths {
  input {
    String root_path
    String num_files
    String path_suffix
  }

  command <<<
    set -e

    python3 /app/generate_ordered_paths.py \
    --root_path ~{root_path} \
    --path_suffix ~{path_suffix} \
    --number ~{num_files} > file_names.txt
  >>>

  output {
    Array[File] paths = read_lines("file_names.txt")
  }

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/variantstore:ah_var_store_20211005_2"
    memory: "3 GB"
    disks: "local-disk 10 HDD"
    cpu: 1
  }
}
