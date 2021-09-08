version 1.0

workflow GvsCreateAltAllele {
  input {
    String dataset_project
    String query_project_id = dataset_project
    String default_dataset

    String? service_account_json_path
  }

  call GetVetTableNames {
    input:
      query_project_id = query_project_id,
      dataset_project_id = dataset_project,
      dataset_name = default_dataset,
      service_account_json_path = service_account_json_path
  }

  call CreateAltAlleleTable {
    input:
      query_project_id = query_project_id,
      dataset_project_id = dataset_project,
      dataset_name = default_dataset,
      service_account_json_path = service_account_json_path
  }

  scatter (idx in range(length(GetVetTableNames.vet_tables))) {
    call PopulateAltAlleleTable {
      input:
        create_table_done = CreateAltAlleleTable.done,
        vet_table_name = GetVetTableNames.vet_tables[idx],
        query_project_id = query_project_id,
        dataset_project_id = dataset_project,
        dataset_name = default_dataset,
        service_account_json_path = service_account_json_path
    }
  }

  output {
    Array[String] vet_tables_loaded = PopulateAltAlleleTable.done
  }
}

task GetVetTableNames {
  meta {
    volatile: true
  }

  input {
    String query_project_id
    String dataset_project_id
    String dataset_name

    String? service_account_json_path
  }

  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

  command <<<
    set -e

    if [ ~{has_service_account_file} = 'true' ]; then
      gsutil cp ~{service_account_json_path} local.service_account.json
      gcloud auth activate-service-account --key-file=local.service_account.json
      gcloud config set project ~{query_project_id}
    fi

    echo "project_id = ~{query_project_id}" > ~/.bigqueryrc
    bq query --location=US --project_id=~{query_project_id} --format=csv --use_legacy_sql=false \
    "SELECT table_name FROM ~{dataset_project_id}.~{dataset_name}.INFORMATION_SCHEMA.TABLES WHERE table_name LIKE 'vet_%' ORDER BY table_name" > vet_tables.csv

    # remove the header row from the CSV file
    sed -i 1d vet_tables.csv
  >>>

  output {
    Array[String] vet_tables = read_lines("vet_tables.csv")
  }

  runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
    memory: "3 GB"
    disks: "local-disk 10 HDD"
    preemptible: 3
    cpu: 1
  }
}

task CreateAltAlleleTable {
  input {
    String query_project_id
    String dataset_project_id
    String dataset_name

    String? service_account_json_path
  }

  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

  command <<<
    set -e

    if [ ~{has_service_account_file} = 'true' ]; then
      gsutil cp ~{service_account_json_path} local.service_account.json
      gcloud auth activate-service-account --key-file=local.service_account.json
      gcloud config set project ~{query_project_id}
    fi

    echo "project_id = ~{query_project_id}" > ~/.bigqueryrc
    bq query --location=US --project_id=~{query_project_id} --format=csv --use_legacy_sql=false \
    "CREATE OR REPLACE TABLE ~{dataset_project_id}.~{dataset_name}.alt_allele (
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
    CLUSTER BY location, sample_id;"

  >>>

  output {
    String done = "done"
  }

  runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
    memory: "3 GB"
    disks: "local-disk 10 HDD"
    cpu: 1
  }
}

task PopulateAltAlleleTable {
  input {
    String create_table_done
    String vet_table_name
    String query_project_id
    String dataset_project_id
    String dataset_name

    String? service_account_json_path
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
      --query_project ~{query_project_id} \
      --vet_table_name ~{vet_table_name} \
      --fq_dataset ~{dataset_project_id}.~{dataset_name} \
      $SERVICE_ACCOUNT_STANZA
  >>>

  output {
    String done = "~{vet_table_name}"
  }

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/variantstore:ah_var_store_20210903"
    memory: "3 GB"
    disks: "local-disk 10 HDD"
    cpu: 1
  }
}
