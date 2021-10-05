version 1.0

workflow GvsMergeScatteredCallsetShards {
  input {
    String input_vcfs_directory_plus_prefix
    Int num_shards
    String output_vcf_base_name
    String output_directory
    String? service_account_json_path
  }

  call GenerateOrderedPaths as VCFpaths {
    input:
      root_path = input_vcfs_directory_plus_prefix,
      num_files = num_shards,
      path_suffix = ".vcf.gz"
  }

  call MergeVCFs {
    input:
      input_vcfs = VCFpaths.paths,
      output_vcf_name = "${output_vcf_base_name}.vcf.gz",
      output_directory = output_directory,
      service_account_json_path = service_account_json_path
  }
}

task MergeVCFs {
  input {
    Array[File] input_vcfs
    String gather_type = "BLOCK"
    String output_vcf_name
    String output_directory
    String? service_account_json_path
    File? gatk_override
  }

  Int disk_size = ceil(size(input_vcfs, "GiB") * 2.5) + 10

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
