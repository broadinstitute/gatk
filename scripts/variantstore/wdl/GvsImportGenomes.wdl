version 1.0

workflow GvsImportGenomes {

  input {
    String dataset_name
    String project_id

    Array[String] external_sample_names
    Array[File] input_vcfs
    Array[File] input_vcf_indexes

    Int? load_data_preemptible_tries
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

  call GetSampleIds {
    input:
      external_sample_names = external_sample_names,
      project_id = project_id,
      dataset_name = dataset_name,
      table_name = "sample_info",
      service_account_json_path = service_account_json_path
  }

  call CheckForDuplicateData {
    input:
      project_id = project_id,
      dataset_name = dataset_name,
      sample_names = external_sample_names,
      service_account_json_path = service_account_json_path
  }

  call CreateFOFNs {
    input:
      input_vcf_list = write_lines(input_vcfs),
      input_vcf_index_list = write_lines(input_vcf_indexes),
      sample_name_list = write_lines(external_sample_names),
      batch_size = 1,
  }

  scatter (i in range(length(CreateFOFNs.vcf_batch_vcf_fofns))) {
    call LoadData {
      input:
        project_id = project_id,
        dataset_name = dataset_name,
        input_vcfs = read_lines(CreateFOFNs.vcf_batch_vcf_fofns[i]),
        input_vcf_indexes = read_lines(CreateFOFNs.vcf_batch_vcf_index_fofns[i]),
        sample_names = read_lines(CreateFOFNs.vcf_sample_name_fofns[i]),
        interval_list = "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list",
        service_account_json_path = service_account_json_path,
        sample_map = GetSampleIds.sample_map,
        drop_state = "FORTY",
        drop_state_includes_greater_than = false,
        preemptible_tries = load_data_preemptible_tries,
        duplicate_check_passed = CheckForDuplicateData.done
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
    Boolean loaded_in_gvs = true
  }
}

task CheckForDuplicateData {
  input {
    String project_id
    String dataset_name
    Array[String] sample_names
    String? service_account_json_path
  }

  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'
  Int num_samples = length(sample_names)

  meta {
    volatile: true
  }

  command <<<
    set -e

    if [ ~{has_service_account_file} = 'true' ]; then
      gsutil cp ~{service_account_json_path} local.service_account.json
      gcloud auth activate-service-account --key-file=local.service_account.json
      gcloud config set project ~{project_id}
    fi

    echo "project_id = ~{project_id}" > ~/.bigqueryrc

    INFO_SCHEMA_TABLE="~{dataset_name}.INFORMATION_SCHEMA.PARTITIONS"
    TEMP_TABLE="~{dataset_name}.sample_dupe_check"
    SAMPLE_INFO_TABLE="~{dataset_name}.sample_info"

    # create a temp table with the sample_names
    bq --project_id=~{project_id} mk ${TEMP_TABLE} "sample_name:STRING"
    NAMES_FILE=~{write_lines(sample_names)}
    bq load --project_id=~{project_id} ${TEMP_TABLE} $NAMES_FILE "sample_name:STRING"

    # check the INFORMATION_SCHEMA.PARTITIONS table to see if any of input sample names/ids have data loaded into their partitions
    # this returns the list of sample names that do already have data loaded
    echo "WITH items as (SELECT s.sample_id, s.sample_name, s.is_loaded FROM \`${TEMP_TABLE}\` t left outer join \`${SAMPLE_INFO_TABLE}\` s on (s.sample_name = t.sample_name)) " >> query.sql
    echo "SELECT i.sample_name FROM \`${INFO_SCHEMA_TABLE}\` p JOIN items i ON (p.partition_id = CAST(i.sample_id AS STRING)) WHERE p.total_logical_bytes > 0 AND (table_name like 'ref_ranges_%' OR table_name like 'vet_%' OR table_name like 'pet_%')" >> query.sql
    echo "UNION DISTINCT "  >> query.sql
    echo "SELECT i.sample_name FROM items i WHERE i.is_loaded = True "  >> query.sql
    echo "UNION DISTINCT "  >> query.sql
    echo "SELECT i.sample_name FROM items i WHERE i.sample_id IN (SELECT sample_id FROM \`~{dataset_name}.sample_load_status\`) "  >> query.sql


    cat query.sql | bq --location=US --project_id=~{project_id} query --format=csv -n ~{num_samples} --use_legacy_sql=false | sed -e '/sample_name/d' > duplicates

      # remove the temp table
      bq --project_id=~{project_id} rm -f -t ${TEMP_TABLE}

    # true if there is data in results
    if [ -s duplicates ]; then
      echo "ERROR: Trying to load samples that have already been loaded"
      cat duplicates
      exit 1
    fi
  >>>
  runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
    memory: "1 GB"
    disks: "local-disk 10 HDD"
    preemptible: 5
    cpu: 1
  }
  output {
    Boolean done = true
    File? duplicates = "duplicates"
  }
}

task CreateFOFNs {
  input {
    File input_vcf_list
    File input_vcf_index_list
    File sample_name_list
    Int batch_size
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
    Array[File] input_vcfs
    Array[File] input_vcf_indexes
    Array[String] sample_names
    File interval_list
    File sample_map
    String? service_account_json_path
    String? drop_state
    Boolean? drop_state_includes_greater_than = false

    Boolean duplicate_check_passed

    String project_id
    String dataset_name

    Boolean load_ref_ranges = true
    Boolean load_pet = false
    Boolean load_vet = true

    # runtime
    Int? preemptible_tries
  }

  String gatk_override = "gs://broad-dsp-spec-ops/scratch/bigquery-jointcalling/jars/kc_ranges_prepare_20220118/gatk-package-4.2.0.0-462-gc0e684c-SNAPSHOT-local.jar"
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
        --ignore-above-gq-threshold ~{drop_state_includes_greater_than} \
        --project-id ~{project_id} \
        --dataset-name ~{dataset_name} \
        --output-type BQ \
        --enable-reference-ranges ~{load_ref_ranges} \
        --enable-pet ~{load_pet} \
        --enable-vet ~{load_vet} \
        -SN ${sample_name} \
        -SNM ~{sample_map} \
        --ref-version 38
    done
  >>>
  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:4.1.7.0"
    maxRetries: 1
    memory: "3.75 GB"
    disks: "local-disk 50 HDD"
    preemptible: select_first([preemptible_tries, 5])
    cpu: 1
  }
  output {
    String done = "true"
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
    Array[String] load_done
    String dataset_name
    String project_id
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

task GetSampleIds {
  meta {
    volatile: true
  }

  input {
    Array[String] external_sample_names
    String project_id
    String dataset_name
    String table_name
    String? service_account_json_path
    Int samples_per_table = 4000

    # runtime
    Int? preemptible_tries
  }

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
    TEMP_TABLE="~{dataset_name}.sample_names_to_load"
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

    bq --project_id=~{project_id} query --format=csv --use_legacy_sql=false -n ~{num_samples} \
      "SELECT sample_id, samples.sample_name FROM \`~{dataset_name}.~{table_name}\` AS samples JOIN \`${TEMP_TABLE}\` AS temp ON samples.sample_name=temp.sample_name" > sample_map

    cut -d, -f1 sample_map > gvs_ids

    ## delete the table that was only needed for this ingest
    bq --project_id=~{project_id} rm -f=true ${TEMP_TABLE}
  >>>
  runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
    memory: "1 GB"
    disks: "local-disk 10 HDD"
    preemptible: select_first([preemptible_tries, 5])
    cpu: 1
  }
  output {
    Int max_table_id = ceil(read_float("max_sample_id"))
    Int min_table_id = ceil(read_float("min_sample_id"))
    File sample_map = "sample_map"
    File gvs_ids = "gvs_ids"
  }
}
