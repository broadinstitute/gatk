version 1.0

import "GvsUtils.wdl" as Utils

workflow GvsExtractCallset {
  input {
    String project_id
    String dataset_name

    String filter_set_name
    String extract_table_prefix
    String query_project = project_id
    Int scatter_count
    File interval_list = "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list"
    File interval_weights_bed = "gs://broad-public-datasets/gvs/weights/gvs_vet_weights_1kb.bed"
    String output_file_base_name = filter_set_name

    Int? extract_maxretries_override
    Int? extract_preemptible_override
    String? output_gcs_dir
    String? service_account_json_path
    Int? split_intervals_disk_size_override
    Int? split_intervals_mem_override
  }

  File reference = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
  File reference_dict = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict"
  File reference_index = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"

  File gatk_override = "gs://broad-dsp-spec-ops/scratch/bigquery-jointcalling/jars/rc-split-intervals-odd-02252022/gatk.jar"
  String fq_cohort_extract_table  = "~{project_id}.~{dataset_name}.~{extract_table_prefix}__DATA"
  String fq_filter_set_info_table = "~{project_id}.~{dataset_name}.filter_set_info"
  String fq_filter_set_site_table = "~{project_id}.~{dataset_name}.filter_set_sites"
  String fq_filter_set_tranches_table = "~{project_id}.~{dataset_name}.filter_set_tranches"
  String fq_ranges_cohort_ref_extract_table = "~{project_id}.~{dataset_name}.~{extract_table_prefix}__REF_DATA"
  String fq_ranges_cohort_vet_extract_table = "~{project_id}.~{dataset_name}.~{extract_table_prefix}__VET_DATA"
  String fq_samples_to_extract_table = "~{project_id}.~{dataset_name}.~{extract_table_prefix}__SAMPLES"
  String fq_ranges_dataset = "~{project_id}.~{dataset_name}"

  Array[String] tables_patterns_for_datetime_check = ["~{extract_table_prefix}__%"]

  call Utils.SplitIntervals {
    input:
      intervals = interval_list,
      ref_fasta = reference,
      ref_fai = reference_index,
      ref_dict = reference_dict,
      interval_weights_bed = interval_weights_bed,
      scatter_count = scatter_count,
      output_gcs_dir = output_gcs_dir,
      split_intervals_disk_size_override = split_intervals_disk_size_override,
      split_intervals_mem_override = split_intervals_mem_override,
      service_account_json_path = service_account_json_path,
  }

  call Utils.GetBQTablesMaxLastModifiedTimestamp {
    input:
      query_project = query_project,
      data_project = project_id,
      data_dataset = dataset_name,
      table_patterns = tables_patterns_for_datetime_check,
      service_account_json_path = service_account_json_path
  }

  call ValidateFilterSetName {
    input:
      query_project = query_project,
      filter_set_name = filter_set_name,
      data_project = project_id,
      data_dataset = dataset_name,
      service_account_json_path = service_account_json_path
  }

  scatter(i in range(length(SplitIntervals.interval_files))) {
    call ExtractTask {
      input:
        gatk_override                      = gatk_override,
        reference                          = reference,
        reference_index                    = reference_index,
        reference_dict                     = reference_dict,
        fq_samples_to_extract_table        = fq_samples_to_extract_table,
        interval_index                     = i,
        intervals                          = SplitIntervals.interval_files[i],
        fq_cohort_extract_table            = fq_cohort_extract_table,
        fq_ranges_cohort_ref_extract_table = fq_ranges_cohort_ref_extract_table,
        fq_ranges_cohort_vet_extract_table = fq_ranges_cohort_vet_extract_table,
        read_project_id                    = query_project,
        mode                               = "RANGES-PREPARED",
        do_not_filter_override             = false,
        fq_ranges_dataset                  = fq_ranges_dataset,
        fq_filter_set_info_table           = fq_filter_set_info_table,
        fq_filter_set_site_table           = fq_filter_set_site_table,
        fq_filter_set_tranches_table       = fq_filter_set_tranches_table,
        filter_set_name                    = filter_set_name,
        filter_set_name_verified           = ValidateFilterSetName.message,
        service_account_json_path          = service_account_json_path,
        drop_state                         = "FORTY",
        output_file                        = "${output_file_base_name}_${i}.vcf.gz",
        output_gcs_dir                     = output_gcs_dir,
        max_last_modified_timestamp        = GetBQTablesMaxLastModifiedTimestamp.max_last_modified_timestamp,
        extract_preemptible_override       = extract_preemptible_override,
        extract_maxretries_override        = extract_maxretries_override,
    }
  }

  call SumBytes {
    input:
      file_sizes_bytes = flatten([ExtractTask.output_vcf_bytes, ExtractTask.output_vcf_index_bytes])
  }

  call CreateManifest {
    input:
      manifest_lines = ExtractTask.manifest,
      output_gcs_dir = output_gcs_dir,
      service_account_json_path = service_account_json_path
  }

  output {
    Array[File] output_vcfs = ExtractTask.output_vcf
    Array[File] output_vcf_indexes = ExtractTask.output_vcf_index
    Float total_vcfs_size_mb = SumBytes.total_mb
    File manifest = CreateManifest.manifest
  }
}

task ValidateFilterSetName {
  input {
    String filter_set_name
    String data_project
    String data_dataset
    String query_project
    String? service_account_json_path
  }

  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

  command <<<
    set -e

    if [ ~{has_service_account_file} = 'true' ]; then
      gsutil cp ~{service_account_json_path} local.service_account.json
      gcloud auth activate-service-account --key-file=local.service_account.json
      gcloud config set project ~{query_project}
    fi

    echo "project_id = ~{query_project}" > ~/.bigqueryrc

    FILTERSETS=$(bq --location=US --project_id=~{query_project} --format=csv query --use_legacy_sql=false "SELECT filter_set_name as available_filter_set_names FROM ~{data_project}.~{data_dataset}.filter_set_info GROUP BY filter_set_name")

    echo "--------\n"
    echo $FILTERSETS
    echo "--------\n"

    if [[ $FILTERSETS == *"~{filter_set_name}"* ]]; then
      echo "Filter set name ~{filter_set_name} found."
    else
      echo "Error -- `~{filter_set_name}` is not an existing filter_set_name. The filter_set_names in ~{data_project}.~{data_dataset}. are: \n"
      echo $FILTERSETS
      exit 1
    fi
  >>>

  output {
    String message = read_string(stdout())
  }

  runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
    memory: "3 GB"
    disks: "local-disk 10 HDD"
    preemptible: 3
    cpu: 1
  }
}

task ExtractTask {
  input {
    File reference
    File reference_index
    File reference_dict

    String fq_samples_to_extract_table

    Int interval_index
    File intervals
    String? drop_state

    String fq_cohort_extract_table
    String fq_ranges_cohort_ref_extract_table
    String fq_ranges_cohort_vet_extract_table
    String read_project_id
    String output_file
    String? output_gcs_dir

    String mode

    Boolean do_not_filter_override
    String fq_ranges_dataset
    String fq_filter_set_info_table
    String fq_filter_set_site_table
    String fq_filter_set_tranches_table
    String? filter_set_name
    String filter_set_name_verified

    # Runtime Options:
    String? service_account_json_path
    File? gatk_override
    Int? extract_preemptible_override
    Int? extract_maxretries_override

    Int? local_sort_max_records_in_ram = 10000000

    # for call-caching -- check if DB tables haven't been updated since the last run
    String max_last_modified_timestamp
  }

  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

  command <<<
    set -e
    export GATK_LOCAL_JAR="~{default="/root/gatk.jar" gatk_override}"

    if [ ~{has_service_account_file} = 'true' ]; then
      gsutil cp ~{service_account_json_path} local.service_account.json
      export GOOGLE_APPLICATION_CREDENTIALS=local.service_account.json
      gcloud auth activate-service-account --key-file=local.service_account.json
      gcloud config set project ~{read_project_id}
    fi

    df -h

    if [ ~{do_not_filter_override} = 'true' ]; then
      FILTERING_ARGS=''
    else
      FILTERING_ARGS='--filter-set-info-table ~{fq_filter_set_info_table}
          --filter-set-site-table ~{fq_filter_set_site_table}
          --tranches-table ~{fq_filter_set_tranches_table}
          --filter-set-name ~{filter_set_name}'
    fi

    if [ ~{mode} = "RANGES-RAW" ]; then
      MODE_ARGS="--mode RANGES --vet-ranges-fq-dataset ~{fq_ranges_dataset} "
    elif [ ~{mode} = "RANGES-PREPARED" ]; then
      MODE_ARGS="--mode RANGES --vet-ranges-extract-fq-table ~{fq_ranges_cohort_vet_extract_table} --ref-ranges-extract-fq-table ~{fq_ranges_cohort_ref_extract_table} "
    else
      MODE_ARGS="--mode PET --cohort-extract-table ~{fq_cohort_extract_table} "
    fi

    gatk --java-options "-Xmx9g" \
      ExtractCohort \
        ${MODE_ARGS} \
        --ref-version 38 \
        -R ~{reference} \
        -O ~{output_file} \
        --local-sort-max-records-in-ram ~{local_sort_max_records_in_ram} \
        --sample-table ~{fq_samples_to_extract_table} \
        ~{"--inferred-reference-state " + drop_state} \
        -L ~{intervals} \
        --project-id ~{read_project_id} \
        ${FILTERING_ARGS}

    # Drop trailing slash if one exists
    OUTPUT_GCS_DIR=$(echo ~{output_gcs_dir} | sed 's/\/$//')

    OUTPUT_FILE_BYTES="$(du -b ~{output_file} | cut -f1)"
    echo ${OUTPUT_FILE_BYTES} > vcf_bytes.txt

    OUTPUT_FILE_INDEX_BYTES="$(du -b ~{output_file}.tbi | cut -f1)"
    echo ${OUTPUT_FILE_INDEX_BYTES} > vcf_index_bytes.txt

    if [ -n "${OUTPUT_GCS_DIR}" ]; then
      gsutil cp ~{output_file} ${OUTPUT_GCS_DIR}/
      gsutil cp ~{output_file}.tbi ${OUTPUT_GCS_DIR}/
      OUTPUT_FILE_DEST="${OUTPUT_GCS_DIR}/~{output_file}"
      OUTPUT_FILE_INDEX_DEST="${OUTPUT_GCS_DIR}/~{output_file}.tbi"
    else
      OUTPUT_FILE_DEST="~{output_file}"
      OUTPUT_FILE_INDEX_DEST="~{output_file}.tbi"
    fi

    # Parent Task will collect manifest lines and create a joined file
    # Currently, the schema is `[interval_number], [output_file_location], [output_file_size_bytes], [output_file_index_location], [output_file_size_bytes]`
    echo ~{interval_index},${OUTPUT_FILE_DEST},${OUTPUT_FILE_BYTES},${OUTPUT_FILE_INDEX_DEST},${OUTPUT_FILE_INDEX_BYTES} >> manifest.txt
  >>>
  runtime {
    docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:varstore_d8a72b825eab2d979c8877448c0ca948fd9b34c7_change_to_hwe"
    memory: "12 GB"
    disks: "local-disk 150 HDD"
    bootDiskSizeGb: 15
    preemptible: select_first([extract_preemptible_override, "2"])
    maxRetries: select_first([extract_maxretries_override, "3"])
    cpu: 2
  }

  # files sizes are floats instead of ints because they can be larger
  output {
    File output_vcf = "~{output_file}"
    Float output_vcf_bytes = read_float("vcf_bytes.txt")
    File output_vcf_index = "~{output_file}.tbi"
    Float output_vcf_index_bytes = read_float("vcf_index_bytes.txt")
    String manifest = read_string("manifest.txt")
  }
}

task SumBytes {
  input {
    Array[Float] file_sizes_bytes
  }

  command <<<
    set -e
    echo "~{sep=" " file_sizes_bytes}" | tr " " "\n" | python -c "
    import sys;
    total_bytes = sum(float(i.strip()) for i in sys.stdin);
    total_mb = total_bytes/10**6;
    print(total_mb);"
  >>>
  runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
    memory: "3 GB"
    disks: "local-disk 10 HDD"
    preemptible: 3
    cpu: 1
  }

  output {
    Float total_mb = read_float(stdout())
  }
}

task CreateManifest {
  input {
      Array[String] manifest_lines
      String? output_gcs_dir
      String? service_account_json_path
  }

  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

  command <<<
    set -e
    MANIFEST_LINES_TXT=~{write_lines(manifest_lines)}
    echo "vcf_file_location, vcf_file_bytes, vcf_index_location, vcf_index_bytes" >> manifest.txt
    sort -n ${MANIFEST_LINES_TXT} | cut -d',' -f 2- >> manifest.txt

    # Drop trailing slash if one exists
    OUTPUT_GCS_DIR=$(echo ~{output_gcs_dir} | sed 's/\/$//')

    if [ -n "$OUTPUT_GCS_DIR" ]; then
      if [ ~{has_service_account_file} = 'true' ]; then
        gsutil cp ~{service_account_json_path} local.service_account.json
        gcloud auth activate-service-account --key-file=local.service_account.json
      fi
      gsutil cp manifest.txt ${OUTPUT_GCS_DIR}/
    fi
  >>>
  output {
    File manifest = "manifest.txt"
  }

  runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
    memory: "3 GB"
    disks: "local-disk 10 HDD"
    preemptible: 3
    cpu: 1
  }
}
