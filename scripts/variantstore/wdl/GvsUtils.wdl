version 1.0

task MergeVCFs {
  input {
    Array[File] input_vcfs
    String gather_type = "BLOCK"
    String output_vcf_name
    String? output_directory
    Int? merge_disk_override
    Int? preemptible_tries
    String? service_account_json_path
  }

  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'
  Int disk_size = if (defined(merge_disk_override)) then merge_disk_override else 100

  parameter_meta {
    input_vcfs: {
      localization_optional: true
    }
  }

  command {
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

    if [ -n "$OUTPUT_GCS_DIR" ]; then
      gsutil cp ~{output_vcf_name} $OUTPUT_GCS_DIR/
      gsutil cp ~{output_vcf_name}.tbi $OUTPUT_GCS_DIR/
    fi
  }

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:varstore_d8a72b825eab2d979c8877448c0ca948fd9b34c7_change_to_hwe"
    preemptible: select_first([preemptible_tries, 3])
    memory: "3 GiB"
    disks: "local-disk ~{disk_size} HDD"
  }

  output {
    File output_vcf = "~{output_vcf_name}"
    File output_vcf_index = "~{output_vcf_name}.tbi"
  }
}

task SplitIntervals {
  input {
    File intervals
    File ref_fasta
    File ref_fai
    File ref_dict
    Int scatter_count
    File? interval_weights_bed
    String? intervals_file_extension
    String? split_intervals_extra_args
    Int? split_intervals_disk_size_override
    Int? split_intervals_mem_override
    String? output_gcs_dir
    File? gatk_override = "gs://broad-dsp-spec-ops/scratch/bigquery-jointcalling/jars/rc-split-intervals-odd-02252022/gatk.jar"
    String? service_account_json_path
  }

  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

  Int disk_size = if (defined(split_intervals_disk_size_override)) then split_intervals_disk_size_override else 10
  Int disk_memory = if (defined(split_intervals_mem_override)) then split_intervals_mem_override else 16
  Int java_memory = disk_memory - 4

  String gatkTool = if (defined(interval_weights_bed)) then 'WeightedSplitIntervals' else 'SplitIntervals'

  parameter_meta {
    intervals: {
      localization_optional: true
    }
    ref_fasta: {
      localization_optional: true
    }
    ref_fai: {
      localization_optional: true
    }
    ref_dict: {
      localization_optional: true
    }
  }

  command {
    set -e
    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

    mkdir interval-files
    gatk --java-options "-Xmx~{java_memory}g" ~{gatkTool} \
      --dont-mix-contigs \
      -R ~{ref_fasta} \
      ~{"-L " + intervals} \
      ~{"--weight-bed-file " + interval_weights_bed} \
      -scatter ~{scatter_count} \
      -O interval-files \
      ~{"--extension " + intervals_file_extension} \
      --interval-file-num-digits 10 \
      ~{split_intervals_extra_args}
    cp interval-files/*.interval_list .

    # Drop trailing slash if one exists
    OUTPUT_GCS_DIR=$(echo ~{output_gcs_dir} | sed 's/\/$//')

    if [ ~{has_service_account_file} = 'true' ]; then
      gsutil cp ~{service_account_json_path} local.service_account.json
      gcloud auth activate-service-account --key-file=local.service_account.json
    fi

    if [ -n "$OUTPUT_GCS_DIR" ]; then
      gsutil -m cp *.interval_list $OUTPUT_GCS_DIR/
    fi
  }

  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:4.2.3.0"
    bootDiskSizeGb: 15
    memory: "~{disk_memory} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: 3
    cpu: 1
  }

  output {
    Array[File] interval_files = glob("*.interval_list")
  }
}

task GetBQTableLastModifiedDatetime {
  # because this is being used to determine if the data has changed, never use call cache
  meta {
    volatile: true
  }

  input {
    String query_project
    String fq_table
    String? service_account_json_path
  }

  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

  # ------------------------------------------------
  # try to get the last modified date for the table in question; fail if something comes back from BigQuery
  # that isn't in the right format (e.g. an error)
  command <<<
    set -e

    if [ ~{has_service_account_file} = 'true' ]; then
      gsutil cp ~{service_account_json_path} local.service_account.json
      gcloud auth activate-service-account --key-file=local.service_account.json
      gcloud config set project ~{query_project}
    fi

    echo "project_id = ~{query_project}" > ~/.bigqueryrc

    # bq needs the project name to be separate by a colon
    DATASET_TABLE_COLON=$(echo ~{fq_table} | sed 's/\./:/')

    LASTMODIFIED=$(bq --location=US --project_id=~{query_project} --format=json show ${DATASET_TABLE_COLON} | python3 -c "import sys, json; print(json.load(sys.stdin)['lastModifiedTime']);")
    if [[ $LASTMODIFIED =~ ^[0-9]+$ ]]; then
      echo $LASTMODIFIED
    else
      exit 1
    fi
  >>>

  output {
    String last_modified_timestamp = read_string(stdout())
  }

  runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
    memory: "3 GB"
    disks: "local-disk 10 HDD"
    preemptible: 3
    cpu: 1
  }
}

task GetBQTablesMaxLastModifiedTimestamp {
  # because this is being used to determine if the data has changed, never use call cache
  meta {
    volatile: true
  }

  input {
    String query_project
    String data_project
    String data_dataset
    Array[String] table_patterns
    String? service_account_json_path
  }

  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

  # ------------------------------------------------
  # try to get the latest last modified timestamp, in epoch microseconds, for all of the tables that match the provided prefixes
  command <<<
    set -e
    if [ ~{has_service_account_file} = 'true' ]; then
    gsutil cp ~{service_account_json_path} local.service_account.json
    gcloud auth activate-service-account --key-file=local.service_account.json
    fi

    echo "project_id = ~{query_project}" > ~/.bigqueryrc

    bq --location=US --project_id=~{query_project} query --format=csv --use_legacy_sql=false \
    "SELECT UNIX_MICROS(MAX(last_modified_time)) last_modified_time FROM \`~{data_project}\`.~{data_dataset}.INFORMATION_SCHEMA.PARTITIONS WHERE table_name like '~{sep="' OR table_name like '" table_patterns}'" > results.txt

    tail -1 results.txt | cut -d, -f1 > max_last_modified_timestamp.txt
  >>>

  output {
    String max_last_modified_timestamp = read_string("max_last_modified_timestamp.txt")
  }

  runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
    memory: "3 GB"
    disks: "local-disk 10 HDD"
    preemptible: 3
    cpu: 1
  }
}
