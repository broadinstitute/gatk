version 1.0

import "GvsUtils.wdl" as Utils

workflow GvsCreateAltAllele {
  input {
    Boolean go = true
    String dataset_name
    String project_id
    String call_set_identifier
    Boolean? skip_create_alt_allele_table = false

    String? service_account_json_path
  }

  String fq_alt_allele_table = "~{project_id}.~{dataset_name}.alt_allele"

  call GetMaxSampleId {
    input:
      dataset_name = dataset_name,
      project_id = project_id
  }

  call GetVetTableNames {
    input:
      dataset_name = dataset_name,
      project_id = project_id,
      service_account_json_path = service_account_json_path,
      max_sample_id = GetMaxSampleId.max_sample_id
  }

  if (!defined(skip_create_alt_allele_table)) {
    call CreateAltAlleleTable {
      input:
        dataset_name = dataset_name,
        project_id = project_id,
        service_account_json_path = service_account_json_path
    }
  }

  call Utils.GetBQTableLastModifiedDatetime {
    input:
      go = CreateAltAlleleTable.done,
      query_project = project_id,
      fq_table = fq_alt_allele_table,
      service_account_json_path = service_account_json_path
  }

  scatter (idx in range(length(GetVetTableNames.vet_tables))) {
    call PopulateAltAlleleTable {
      input:
        call_set_identifier = call_set_identifier,
        dataset_name = dataset_name,
        project_id = project_id,
        create_table_done = select_first(skip_create_alt_allele_table, CreateAltAlleleTable.done),
        vet_table_name = GetVetTableNames.vet_tables[idx],
        service_account_json_path = service_account_json_path,
        last_modified_timestamp = GetBQTableLastModifiedDatetime.last_modified_timestamp,
        max_sample_id = GetMaxSampleId.max_sample_id
    }
  }

  output {
    Array[String] vet_tables_loaded = PopulateAltAlleleTable.done
    Boolean done = true
  }
}

task GetMaxSampleId {
  input {
    String dataset_name
    String project_id
  }
  meta {
    # because this is being used to determine the current state of the GVS database, never use call cache
    volatile: true
  }

  command <<<
    set -e

    echo "project_id = ~{project_id}" > ~/.bigqueryrc
    bq query --location=US --project_id=~{project_id} --format=csv --use_legacy_sql=false \
    'SELECT MAX(sample_id) as max_sample_id FROM `~{dataset_name}.sample_info`' > num_rows.csv

    NUMROWS=$(python3 -c "csvObj=open('num_rows.csv','r');csvContents=csvObj.read();print(csvContents.split('\n')[1]);")

    [[ $NUMROWS =~ ^[0-9]+$ ]] && echo $NUMROWS || exit 1
  >>>

  output {
    Int max_sample_id = read_int(stdout())
  }

  runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
    memory: "3 GB"
    disks: "local-disk 10 HDD"
    preemptible: 3
    cpu: 1
  }
}

task GetVetTableNames {
  input {
    String dataset_name
    String project_id
    Int max_sample_id

    String? service_account_json_path
  }

  meta {
    # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
  }

  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'
  # add labels for DSP Cloud Cost Control Labeling and Reporting
  String bq_labels = "--label service:gvs --label team:variants --label managedby:create_alt_allele"

  command <<<
    set -e

    if [ ~{has_service_account_file} = 'true' ]; then
      gsutil cp ~{service_account_json_path} local.service_account.json
      gcloud auth activate-service-account --key-file=local.service_account.json
      gcloud config set project ~{project_id}
    fi

    if [ $(i % 4000 == 0) -eq 0 ]; then
      echo $(~{max_sample_id} / 4000) > min_vat_table_num.txt
    else
      echo $((~{max_sample_id} / 4000) + 1) > min_vat_table_num.txt
    fi

    echo "project_id = ~{project_id}" > ~/.bigqueryrc
    bq query --location=US --project_id=~{project_id} --format=csv --use_legacy_sql=false ~{bq_labels} \
    'SELECT table_name FROM `~{project_id}.~{dataset_name}.INFORMATION_SCHEMA.TABLES` WHERE table_name LIKE "vet_%" AND CAST(SUBSTRING(table_name, 5) AS INT64) >= $(cat min_vat_table_num.txt)' > vet_tables.csv

    # remove the header row from the CSV file
    sed -i 1d vet_tables.csv
  >>>
  runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
    memory: "3 GB"
    disks: "local-disk 10 HDD"
    preemptible: 3
    cpu: 1
  }

  output {
    Array[String] vet_tables = read_lines("vet_tables.csv")
  }
}

task CreateAltAlleleTable {
  input {
    Boolean go = true
    String dataset_name
    String project_id

    String? service_account_json_path
  }
  meta {
    # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
  }

  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'
  # add labels for DSP Cloud Cost Control Labeling and Reporting
  String bq_labels = "--label service:gvs --label team:variants --label managedby:create_alt_allele"

  command <<<
    set -e

    if [ ~{has_service_account_file} = 'true' ]; then
      gsutil cp ~{service_account_json_path} local.service_account.json
      gcloud auth activate-service-account --key-file=local.service_account.json
      gcloud config set project ~{project_id}
    fi

    echo "project_id = ~{project_id}" > ~/.bigqueryrc
    bq query --location=US --project_id=~{project_id} --format=csv --use_legacy_sql=false ~{bq_labels} \
    'CREATE OR REPLACE TABLE `~{project_id}.~{dataset_name}.alt_allele` (
      location INT64,
      sample_id INT64,
      ref STRING,
      allele STRING,
      allele_pos INT64,
      call_GT STRING,
      call_GQ INT64,
      as_raw_mq STRING,
      raw_mq INT64,
      as_raw_mqranksum STRING,
      raw_mqranksum_x_10 INT64,
      as_qualapprox STRING,
      qualapprox STRING,
      qual INT64,
      as_raw_readposranksum STRING,
      raw_readposranksum_x_10 INT64,
      as_sb_table STRING,
      sb_ref_plus INT64,
      sb_ref_minus INT64,
      sb_alt_plus INT64,
      sb_alt_minus INT64,
      call_AD STRING,
      ref_ad INT64,
      ad INT64
    ) PARTITION BY RANGE_BUCKET(location, GENERATE_ARRAY(0, 25000000000000, 1000000000000))
    CLUSTER BY location, sample_id;'

  >>>
  runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
    memory: "3 GB"
    disks: "local-disk 10 HDD"
    cpu: 1
  }

  output {
    Boolean done = true
  }
}

task PopulateAltAlleleTable {
  input {
    String dataset_name
    String project_id

    String create_table_done
    String vet_table_name
    String call_set_identifier
    Int max_sample_id

    String? service_account_json_path

    String last_modified_timestamp
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
      SERVICE_ACCOUNT_STANZA="--sa_key_path local.service_account.json "
    fi

    python3 /app/populate_alt_allele_table.py \
      --call_set_identifier ~{call_set_identifier} \
      --query_project ~{project_id} \
      --vet_table_name ~{vet_table_name} \
      --fq_dataset ~{project_id}.~{dataset_name} \
      --max_sample_id ~{max_sample_id} \
      $SERVICE_ACCOUNT_STANZA
  >>>
  runtime {
    docker: "us.gcr.io/broad-dsde-methods/variantstore:rsa_vs_52_incremental_alt_allele_2022_08_16"
    memory: "3 GB"
    disks: "local-disk 10 HDD"
    cpu: 1
  }

  output {
    String done = "~{vet_table_name}"
  }
}
