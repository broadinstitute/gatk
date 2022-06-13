version 1.0

import "GvsUtils.wdl" as Utils

workflow GvsCreateAltAllele {
  input {
    Boolean go = true
    String dataset_name
    String project_id

    String? service_account_json_path
  }

  String fq_alt_allele_table = "~{project_id}.~{dataset_name}.alt_allele"

  call GetVetTableNames {
    input:
      dataset_name = dataset_name,
      project_id = project_id,
      service_account_json_path = service_account_json_path
  }

  call CreateAltAlleleTable {
    input:
      dataset_name = dataset_name,
      project_id = project_id,
      service_account_json_path = service_account_json_path
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
        dataset_name = dataset_name,
        project_id = project_id,
        create_table_done = CreateAltAlleleTable.done,
        vet_table_name = GetVetTableNames.vet_tables[idx],
        service_account_json_path = service_account_json_path,
        last_modified_timestamp = GetBQTableLastModifiedDatetime.last_modified_timestamp
    }
  }

  output {
    Array[String] vet_tables_loaded = PopulateAltAlleleTable.done
    Boolean done = true
  }
}

task GetVetTableNames {
  input {
    String dataset_name
    String project_id

    String? service_account_json_path
  }

  meta {
    # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
  }

  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'
  # add labels for DSP Cloud Cost Control Labeling and Reporting
  String bq_labels = "--label service:gvs --label team:variants --label mangedby:createaltallele"

  command <<<
    set -e

    if [ ~{has_service_account_file} = 'true' ]; then
      gsutil cp ~{service_account_json_path} local.service_account.json
      gcloud auth activate-service-account --key-file=local.service_account.json
      gcloud config set project ~{project_id}
    fi

    echo "project_id = ~{project_id}" > ~/.bigqueryrc
    bq query --location=US --project_id=~{project_id} --format=csv --use_legacy_sql=false ~{bq_labels} \
    'SELECT table_name FROM `~{project_id}.~{dataset_name}.INFORMATION_SCHEMA.TABLES` WHERE table_name LIKE "vet_%" ORDER BY table_name' > vet_tables.csv

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
  String bq_labels = "--label service:gvs --label team:variants --label mangedby:createaltallele"

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
      --query_project ~{project_id} \
      --vet_table_name ~{vet_table_name} \
      --fq_dataset ~{project_id}.~{dataset_name} \
      $SERVICE_ACCOUNT_STANZA
  >>>
  runtime {
    docker: "us.gcr.io/broad-dsde-methods/variantstore:ah_var_store_2022_05_13"
    memory: "3 GB"
    disks: "local-disk 10 HDD"
    cpu: 1
  }

  output {
    String done = "~{vet_table_name}"
  }
}
