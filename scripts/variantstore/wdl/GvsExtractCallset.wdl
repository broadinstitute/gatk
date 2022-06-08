version 1.0

import "GvsUtils.wdl" as Utils

workflow GvsExtractCallset {
  input {
    Boolean go = true
    String dataset_name
    String project_id

    String cohort_project_id = project_id
    String cohort_dataset_name = dataset_name
    Boolean do_not_filter_override = false
    Boolean control_samples = false
    String extract_table_prefix
    String filter_set_name
    String query_project = project_id
    # This is optional now since the workflow will choose an appropriate value below if this is unspecified.
    Int? scatter_count
    Boolean zero_pad_output_vcf_filenames = true

    File interval_list = "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list"
    File interval_weights_bed = "gs://broad-public-datasets/gvs/weights/gvs_vet_weights_1kb.bed"
    File gatk_override = "gs://broad-dsp-spec-ops/scratch/bigquery-jointcalling/jars/rc-add-withdrawn-back/gatk-package-4.2.0.0-527-ge5b9dd6-SNAPSHOT-local.jar"

    String output_file_base_name = filter_set_name

    Int? extract_maxretries_override
    Int? extract_preemptible_override
    String? output_gcs_dir
    String? service_account_json_path
    Int? split_intervals_disk_size_override
    Int? split_intervals_mem_override
    Float x_bed_weight_scaling = 4
    Float y_bed_weight_scaling = 4
  }

  File reference = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
  File reference_dict = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict"
  File reference_index = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"

  String fq_gvs_dataset = "~{project_id}.~{dataset_name}"
  String fq_cohort_dataset = "~{cohort_project_id}.~{cohort_dataset_name}"

  String full_extract_prefix = if (control_samples) then "~{extract_table_prefix}_controls" else extract_table_prefix
  String fq_filter_set_info_table = "~{fq_gvs_dataset}.filter_set_info"
  String fq_filter_set_site_table = "~{fq_gvs_dataset}.filter_set_sites"
  String fq_filter_set_tranches_table = "~{fq_gvs_dataset}.filter_set_tranches"
  String fq_sample_table = "~{fq_gvs_dataset}.sample_info"
  String fq_cohort_extract_table = "~{fq_cohort_dataset}.~{full_extract_prefix}__DATA"
  String fq_ranges_cohort_ref_extract_table = "~{fq_cohort_dataset}.~{full_extract_prefix}__REF_DATA"
  String fq_ranges_cohort_vet_extract_table = "~{fq_cohort_dataset}.~{full_extract_prefix}__VET_DATA"

  String fq_samples_to_extract_table = "~{fq_cohort_dataset}.~{full_extract_prefix}__SAMPLES"
  Array[String] tables_patterns_for_datetime_check = ["~{full_extract_prefix}__%"]

  Boolean emit_pls = false
  Boolean emit_ads = true

  String intervals_file_extension = if (zero_pad_output_vcf_filenames) then '-~{output_file_base_name}.vcf.gz.interval_list' else '-scattered.interval_list'

  call Utils.ScaleXYBedValues {
    input:
      interval_weights_bed = interval_weights_bed,
      x_bed_weight_scaling = x_bed_weight_scaling,
      y_bed_weight_scaling = y_bed_weight_scaling
  }

  call Utils.GetBQTableLastModifiedDatetime as SamplesTableDatetimeCheck {
    input:
      query_project = project_id,
      fq_table = fq_sample_table,
      service_account_json_path = service_account_json_path
  }

  call Utils.GetNumSamplesLoaded {
    input:
      fq_sample_table = fq_sample_table,
      fq_sample_table_lastmodified_timestamp = SamplesTableDatetimeCheck.last_modified_timestamp,
      service_account_json_path = service_account_json_path,
      project_id = project_id,
      control_samples = control_samples
  }

  Int effective_scatter_count = if defined(scatter_count) then select_first([scatter_count])
                                else if GetNumSamplesLoaded.num_samples < 100 then 100 # Quickstart
                                     else if GetNumSamplesLoaded.num_samples < 1000 then 500
                                          else if GetNumSamplesLoaded.num_samples < 5000 then 1000
                                               else if GetNumSamplesLoaded.num_samples < 20000 then 2000 # Stroke Anderson
                                                    else if GetNumSamplesLoaded.num_samples < 50000 then 10000
                                                         else if GetNumSamplesLoaded.num_samples < 100000 then 20000 # Charlie
                                                              else 40000

  call Utils.SplitIntervals {
    input:
      intervals = interval_list,
      ref_fasta = reference,
      ref_fai = reference_index,
      ref_dict = reference_dict,
      interval_weights_bed = ScaleXYBedValues.xy_scaled_bed,
      intervals_file_extension = intervals_file_extension,
      scatter_count = effective_scatter_count,
      output_gcs_dir = output_gcs_dir,
      split_intervals_disk_size_override = split_intervals_disk_size_override,
      split_intervals_mem_override = split_intervals_mem_override,
      gatk_override = gatk_override,
      service_account_json_path = service_account_json_path,
  }

  call Utils.GetBQTableLastModifiedDatetime as FilterSetInfoTimestamp {
       input:
       query_project = project_id,
       fq_table = "filter_set_info",
       service_account_json_path = service_account_json_path,
  }

  if ( !do_not_filter_override ) {
    call ValidateFilterSetName {
      input:
      query_project = query_project,
      filter_set_name = filter_set_name,
      filter_set_info_timestamp = FilterSetInfoTimestamp.last_modified_timestamp,
      data_project = project_id,
      data_dataset = dataset_name,
      service_account_json_path = service_account_json_path
    }
  }

  call Utils.GetBQTablesMaxLastModifiedTimestamp {
    input:
      query_project = query_project,
      data_project = project_id,
      data_dataset = dataset_name,
      table_patterns = tables_patterns_for_datetime_check,
      service_account_json_path = service_account_json_path
  }

  scatter(i in range(length(SplitIntervals.interval_files))) {
    String interval_filename = basename(SplitIntervals.interval_files[i])
    String vcf_filename = if (zero_pad_output_vcf_filenames) then sub(interval_filename, ".interval_list", "") else "~{output_file_base_name}_${i}.vcf.gz"

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
        do_not_filter_override             = do_not_filter_override,
        fq_filter_set_info_table           = fq_filter_set_info_table,
        fq_filter_set_site_table           = fq_filter_set_site_table,
        fq_filter_set_tranches_table       = fq_filter_set_tranches_table,
        filter_set_name                    = filter_set_name,
        filter_set_name_verified           = select_first([ValidateFilterSetName.done, "done"]),
        service_account_json_path          = service_account_json_path,
        drop_state                         = "FORTY",
        output_file                        = vcf_filename,
        output_gcs_dir                     = output_gcs_dir,
        max_last_modified_timestamp        = GetBQTablesMaxLastModifiedTimestamp.max_last_modified_timestamp,
        extract_preemptible_override       = extract_preemptible_override,
        extract_maxretries_override        = extract_maxretries_override,
        emit_pls                           = emit_pls,
        emit_ads                           = emit_ads,
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

  if (control_samples == false) {

    call Utils.GetBQTableLastModifiedDatetime {
      input:
        query_project = query_project,
        fq_table = fq_samples_to_extract_table,
        service_account_json_path = service_account_json_path
    }

    call GenerateSampleListFile {
      input:
        fq_samples_to_extract_table = fq_samples_to_extract_table,
        samples_to_extract_table_timestamp = GetBQTableLastModifiedDatetime.last_modified_timestamp,
        output_gcs_dir = output_gcs_dir,
        query_project = query_project,
        service_account_json_path = service_account_json_path
    }
  }

  output {
    Array[File] output_vcfs = ExtractTask.output_vcf
    Array[File] output_vcf_indexes = ExtractTask.output_vcf_index
    Float total_vcfs_size_mb = SumBytes.total_mb
    File manifest = CreateManifest.manifest
    Boolean done = true
  }
}

task ValidateFilterSetName {
  input {
    String filter_set_name
    String data_project
    String data_dataset
    String query_project
    String? service_account_json_path
    String filter_set_info_timestamp
  }
  meta {
    # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
  }

  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

  command <<<
    set -ex

    if [ ~{has_service_account_file} = 'true' ]; then
      gsutil cp ~{service_account_json_path} local.service_account.json
      gcloud auth activate-service-account --key-file=local.service_account.json
      gcloud config set project ~{query_project}
    fi

    echo "project_id = ~{query_project}" > ~/.bigqueryrc

    OUTPUT=$(bq --location=US --project_id=~{query_project} --format=csv query --use_legacy_sql=false "SELECT filter_set_name as available_filter_set_names FROM \`~{data_project}.~{data_dataset}.filter_set_info\` GROUP BY filter_set_name")
    FILTERSETS=${OUTPUT#"available_filter_set_names"}

    if [[ $FILTERSETS =~ "~{filter_set_name}" ]]; then
      echo "Filter set name '~{filter_set_name}' found."
    else
      echo "ERROR: '~{filter_set_name}' is not an existing filter_set_name. Available in ~{data_project}.~{data_dataset} are"
      echo $FILTERSETS
      exit 1
    fi
  >>>
  output {
    String done = read_string(stdout())
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

    Boolean emit_pls
    Boolean emit_ads

    Boolean do_not_filter_override
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
  meta {
    # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
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

    gatk --java-options "-Xmx9g" \
      ExtractCohort \
        --vet-ranges-extract-fq-table ~{fq_ranges_cohort_vet_extract_table} \
        --ref-ranges-extract-fq-table ~{fq_ranges_cohort_ref_extract_table} \
        --ref-version 38 \
        -R ~{reference} \
        -O ~{output_file} \
        --local-sort-max-records-in-ram ~{local_sort_max_records_in_ram} \
        --sample-table ~{fq_samples_to_extract_table} \
        ~{"--inferred-reference-state " + drop_state} \
        -L ~{intervals} \
        --project-id ~{read_project_id} \
        ~{true='--emit-pls' false='' emit_pls} \
        ~{true='--emit-ads' false='' emit_ads} \
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
  meta {
    # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
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
  meta {
    # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
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

task GenerateSampleListFile {
  input {
    String fq_samples_to_extract_table
    String samples_to_extract_table_timestamp
    String query_project

    String? output_gcs_dir
    String? service_account_json_path
  }
  meta {
    # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
  }

  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

  command <<<
    set -e

    if [ ~{has_service_account_file} = 'true' ]; then
      gsutil cp ~{service_account_json_path} local.service_account.json
      gcloud auth activate-service-account --key-file=local.service_account.json
      gcloud config set project ~{query_project}
    fi

    # Drop trailing slash if one exists
    OUTPUT_GCS_DIR=$(echo ~{output_gcs_dir} | sed 's/\/$//')

    echo "project_id = ~{query_project}" > ~/.bigqueryrc

    bq --location=US --project_id=~{query_project} --format=csv query --use_legacy_sql=false "SELECT sample_name FROM ~{fq_samples_to_extract_table}" | sed 1d > sample-name-list.txt

    if [ -n "$OUTPUT_GCS_DIR" ]; then
      gsutil cp sample-name-list.txt ${OUTPUT_GCS_DIR}/
    fi
  >>>
  output {
    File sample_name_list = "sample-name-list.txt"
  }

  runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
    memory: "3 GB"
    disks: "local-disk 10 HDD"
    preemptible: 3
    cpu: 1
  }
}
