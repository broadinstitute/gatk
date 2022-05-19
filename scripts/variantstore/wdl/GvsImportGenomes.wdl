version 1.0

workflow GvsImportGenomes {

  input {
    Boolean go = true
    String dataset_name
    String project_id

    Array[String] external_sample_names
    Array[File] input_vcfs
    Array[File] input_vcf_indexes

    File interval_list = "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list"
    Int? load_data_preemptible_override
    Int? load_data_maxretries_override
    File? load_data_gatk_override = "gs://broad-dsp-spec-ops/scratch/bigquery-jointcalling/jars/ah_var_store_20220415/gatk-package-4.2.0.0-492-g1387d47-SNAPSHOT-local.jar"
    String? service_account_json_path
  }

  # return an error if the lengths are not equal
  Int input_length = length(input_vcfs)
  Int input_indexes_length = length(input_vcf_indexes)
  if ((input_length != length(external_sample_names)) || (input_indexes_length != length(external_sample_names))) {
    call TerminateWorkflow {
      input:
        message = "The number of external_sample_names, sample input_vcfs and sample input_vcf_indexes are not the same."
    }
  }

  call GetUningestedSampleIds {
    input:
      dataset_name = dataset_name,
      project_id = project_id,
      external_sample_names = external_sample_names,
      table_name = "sample_info",
      service_account_json_path = service_account_json_path
  }

  call CurateInputLists {
    input:
      dataset_name = dataset_name,
      project_id = project_id,
      input_vcf_index_list = write_lines(input_vcf_indexes),
      input_vcf_list = write_lines(input_vcfs),
      input_sample_name_list = write_lines(external_sample_names),
      input_sample_map = GetUningestedSampleIds.sample_map,
      service_account_json_path = service_account_json_path
  }

  call CreateFOFNs {
    input:
      batch_size = 1,
      input_vcf_index_list = CurateInputLists.index_list,
      input_vcf_list = CurateInputLists.vcf_list,
      sample_name_list = CurateInputLists.sample_name_list
  }

  scatter (i in range(length(CreateFOFNs.vcf_batch_vcf_fofns))) {
    call LoadData {
      input:
        dataset_name = dataset_name,
        project_id = project_id,
        drop_state = "FORTY",
        drop_state_includes_greater_than = false,
        input_vcf_indexes = read_lines(CreateFOFNs.vcf_batch_vcf_index_fofns[i]),
        input_vcfs = read_lines(CreateFOFNs.vcf_batch_vcf_fofns[i]),
        interval_list = interval_list,
        gatk_override = load_data_gatk_override,
        load_data_preemptible_override = load_data_preemptible_override,
        load_data_maxretries_override = load_data_maxretries_override,
        sample_names = read_lines(CreateFOFNs.vcf_sample_name_fofns[i]),
        sample_map = GetUningestedSampleIds.sample_map,
        service_account_json_path = service_account_json_path,
    }
  }

  call SetIsLoadedColumn {
    input:
      load_done = LoadData.done,
      project_id = project_id,
      dataset_name = dataset_name,
      service_account_json_path = service_account_json_path,
  }

  output {
    Boolean done = true
    Array[File] load_data_stderrs = LoadData.stderr
  }
}

task CreateFOFNs {
  input {
    Int batch_size
    File input_vcf_index_list
    File input_vcf_list
    File sample_name_list
  }

  command <<<
    set -e

    split -d -a 5 -l ~{batch_size} ~{input_vcf_list} batched_vcfs.
    split -d -a 5 -l ~{batch_size} ~{input_vcf_index_list} batched_vcf_indexes.
    split -d -a 5 -l ~{batch_size} ~{sample_name_list} batched_sample_names.
  >>>
  runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
    bootDiskSizeGb: 15
    memory: "3 GB"
    disks: "local-disk 10 HDD"
    preemptible: 3
    cpu: 1
  }

  output {
    Array[File] vcf_batch_vcf_fofns = glob("batched_vcfs.*")
    Array[File] vcf_batch_vcf_index_fofns = glob("batched_vcf_indexes.*")
    Array[File] vcf_sample_name_fofns = glob("batched_sample_names.*")
  }
}

task LoadData {
  input {
    String dataset_name
    String project_id

    Array[File] input_vcf_indexes
    Array[File] input_vcfs
    File interval_list
    File sample_map
    Array[String] sample_names

    String? drop_state
    Boolean? drop_state_includes_greater_than = false
    Boolean force_loading_from_non_allele_specific = false

    File? gatk_override
    Int? load_data_preemptible_override
    Int? load_data_maxretries_override
    String? service_account_json_path
  }

  Boolean load_ref_ranges = true
  Boolean load_vet = true
  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

  meta {
    description: "Load data into BigQuery using the Write Api"
    volatile: true
  }

  parameter_meta {
    input_vcfs: {
      localization_optional: true
    }

    input_vcf_indexes: {
      localization_optional: true
    }
  }

  command <<<
    set -e

    # workaround for https://github.com/broadinstitute/cromwell/issues/3647
    export TMPDIR=/tmp

    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

    if [ ~{has_service_account_file} = 'true' ]; then
      gsutil cp ~{service_account_json_path} local.service_account.json
      export GOOGLE_APPLICATION_CREDENTIALS=local.service_account.json
      gcloud auth activate-service-account --key-file=local.service_account.json
    fi

    # translate WDL arrays into BASH arrays
    VCFS_ARRAY=(~{sep=" " input_vcfs})
    VCF_INDEXES_ARRAY=(~{sep=" " input_vcf_indexes})
    SAMPLE_NAMES_ARRAY=(~{sep=" " sample_names})

    # loop over the BASH arrays (See https://stackoverflow.com/questions/6723426/looping-over-arrays-printing-both-index-and-value)
    for i in "${!VCFS_ARRAY[@]}"; do
      input_vcf="${VCFS_ARRAY[$i]}"
      input_vcf_basename=$(basename $input_vcf)
      updated_input_vcf=$input_vcf
      input_vcf_index="${VCF_INDEXES_ARRAY[$i]}"
      sample_name="${SAMPLE_NAMES_ARRAY[$i]}"

      # we always do our own localization
      gsutil cp $input_vcf .
      gsutil cp $input_vcf_index .
      updated_input_vcf=$input_vcf_basename

      gatk --java-options "-Xmx2g" CreateVariantIngestFiles \
        -V ${updated_input_vcf} \
        -L ~{interval_list} \
        ~{"-IG " + drop_state} \
        --force-loading-from-non-allele-specific ~{force_loading_from_non_allele_specific} \
        --ignore-above-gq-threshold ~{drop_state_includes_greater_than} \
        --project-id ~{project_id} \
        --dataset-name ~{dataset_name} \
        --output-type BQ \
        --enable-reference-ranges ~{load_ref_ranges} \
        --enable-vet ~{load_vet} \
        -SN ${sample_name} \
        -SNM ~{sample_map} \
        --ref-version 38
    done
  >>>
  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:4.1.7.0"
    maxRetries: select_first([load_data_maxretries_override, 3])
    memory: "3.75 GB"
    disks: "local-disk 50 HDD"
    preemptible: select_first([load_data_preemptible_override, 5])
    cpu: 1
  }
  output {
    Boolean done = true
    File stderr = stderr()
  }
}

task TerminateWorkflow {
  input {
    String message
  }

  command <<<
    set -e
    echo ~{message}
    exit 1
  >>>
  runtime {
    docker: "python:3.8-slim-buster"
    memory: "1 GB"
    disks: "local-disk 10 HDD"
    preemptible: 3
    cpu: 1
  }
}

task SetIsLoadedColumn {
  meta {
    volatile: true
  }

  input {
    String dataset_name
    String project_id

    Array[String] load_done

    String? service_account_json_path
  }

  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

  command <<<
    set -ex

    if [ ~{has_service_account_file} = 'true' ]; then
      gsutil cp ~{service_account_json_path} local.service_account.json
      gcloud auth activate-service-account --key-file=local.service_account.json
      gcloud config set project ~{project_id}
    fi

    echo "project_id = ~{project_id}" > ~/.bigqueryrc

    # set is_loaded to true if there is a corresponding vet table partition with rows for that sample_id
    bq --location=US --project_id=~{project_id} query --format=csv --use_legacy_sql=false \
    'UPDATE `~{dataset_name}.sample_info` SET is_loaded = true WHERE sample_id IN (SELECT CAST(partition_id AS INT64) from `~{dataset_name}.INFORMATION_SCHEMA.PARTITIONS` WHERE partition_id NOT LIKE "__%" AND total_logical_bytes > 0 AND table_name LIKE "vet_%") OR sample_id IN (SELECT sample_id FROM `~{dataset_name}.sample_load_status` GROUP BY 1 HAVING COUNT(1) = 2)'
  >>>
  runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
    memory: "1 GB"
    disks: "local-disk 10 HDD"
    cpu: 1
  }

  output {
    String done = "done"
  }
}

task GetUningestedSampleIds {
  meta {
    volatile: true
  }

  input {
    String dataset_name
    String project_id

    Array[String] external_sample_names
    String table_name

    String? service_account_json_path
 }

  Int samples_per_table = 4000
  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'
  Int num_samples = length(external_sample_names)

  command <<<
    set -ex
    if [ ~{has_service_account_file} = 'true' ]; then
      gsutil cp ~{service_account_json_path} local.service_account.json
      gcloud auth activate-service-account --key-file=local.service_account.json
    fi

    echo "project_id = ~{project_id}" > ~/.bigqueryrc

    # create temp table with the sample_names and load external sample names into temp table -- make sure it doesn't exist already
    set +e
    TEMP_TABLE="~{dataset_name}.sample_names_to_load"
    bq show --project_id ~{project_id} ${TEMP_TABLE} > /dev/null
    BQ_SHOW_RC=$?
    set -e

    # if there is already a table of sample names or something else is wrong, bail
    if [ $BQ_SHOW_RC -eq 0 ]; then
      echo "There is already a list of sample names. This may need manual cleanup. Exiting"
      exit 1
    fi

    echo "Creating the external sample name list table ${TEMP_TABLE}"
    bq --project_id=~{project_id} mk ${TEMP_TABLE} "sample_name:STRING"
    NAMES_FILE=~{write_lines(external_sample_names)}
    bq load --project_id=~{project_id} ${TEMP_TABLE} $NAMES_FILE "sample_name:STRING"

    # get the current maximum id, or 0 if there are none
    bq --project_id=~{project_id} query --format=csv --use_legacy_sql=false \
      "SELECT IFNULL(MIN(sample_id),0) as min, IFNULL(MAX(sample_id),0) as max FROM \`~{dataset_name}.~{table_name}\` AS samples JOIN \`${TEMP_TABLE}\` AS temp ON samples.sample_name=temp.sample_name" > results

    # prep for being able to return min table id
    min_sample_id=$(tail -1 results | cut -d, -f1)
    max_sample_id=$(tail -1 results | cut -d, -f2)

    # no samples have been loaded or we don't have the right external_sample_names or something else is wrong, bail
    if [ $max_sample_id -eq 0 ]; then
      echo "Max id is 0. Exiting"
      exit 1
    fi

    python3 -c "from math import ceil; print(ceil($max_sample_id/~{samples_per_table}))" > max_sample_id
    python3 -c "from math import ceil; print(ceil($min_sample_id/~{samples_per_table}))" > min_sample_id

    # get sample map of samples that haven't been loaded yet
    bq --project_id=~{project_id} query --format=csv --use_legacy_sql=false -n ~{num_samples} \
      "SELECT sample_id, samples.sample_name FROM \`~{dataset_name}.~{table_name}\` AS samples JOIN \`${TEMP_TABLE}\` AS temp ON samples.sample_name=temp.sample_name WHERE samples.sample_id NOT IN (SELECT sample_id FROM \`~{dataset_name}.sample_load_status\` WHERE status='FINISHED')" > sample_map

    cut -d, -f1 sample_map > gvs_ids

    ## delete the table that was only needed for this ingest
    bq --project_id=~{project_id} rm -f=true ${TEMP_TABLE}
  >>>
  runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
    memory: "1 GB"
    disks: "local-disk 10 HDD"
    preemptible: 5
    cpu: 1
  }
  output {
    Int max_table_id = ceil(read_float("max_sample_id"))
    Int min_table_id = ceil(read_float("min_sample_id"))
    File sample_map = "sample_map"
    File gvs_ids = "gvs_ids"
  }
}

task CurateInputLists {
  input {
    String dataset_name
    String project_id
    File input_vcf_index_list
    File input_vcf_list
    File input_sample_map
    File input_sample_name_list

    String? service_account_json_path
  }

  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'
  command <<<
    set -ex
    if [ ~{has_service_account_file} = 'true' ]; then
      gsutil cp ~{service_account_json_path} local.service_account.json
      gcloud auth activate-service-account --key-file=local.service_account.json
    fi

    gsutil cp  ~{input_vcf_index_list} input_vcf_index_list
    gsutil cp  ~{input_vcf_list} input_vcf_list
    gsutil cp  ~{input_sample_map} input_sample_map
    gsutil cp  ~{input_sample_name_list} input_sample_name_list

    python3 /app/curate_input_array_files.py
  >>>
  runtime {
    docker: "us.gcr.io/broad-dsde-methods/variantstore:rsa_skip_samples_20220519"
    memory: "3 GB"
    disks: "local-disk 100 HDD"
    bootDiskSizeGb: 15
    preemptible: 0
    cpu: 1
  }

  output {
    File index_list = "output_vcf_index_list"
    File vcf_list = "output_vcf_list"
    File sample_name_list = "output_sample_name_list"
  }
}
