version 1.0

workflow GvsCreateAltAllele {
  input {
    String query_project_id
    String dataset_project_id = query_project_id
    String dataset_name

    String? service_account_json_path
    Int? preemptible_tries
    File? gatk_override
    String? docker
  }

  String docker_final = select_first([docker, "us.gcr.io/broad-gatk/gatk:4.1.7.0"])

  call GetVetTableNames {
    input:
      query_project_id = query_project_id,
      dataset_project_id = dataset_project_id,
      dataset_name = dataset_name,
      service_account_json_path = service_account_json_path
  }

  scatter (vet_table in range(length(GetVetTableNames.vet_tables))) {
    call PopulateAltAlleleTable {
      input:
        vet_table_name = vet_table,
        query_project_id = query_project_id,
        dataset_project_id = dataset_project_id,
        dataset_name = dataset_name,
        service_account_json_path = service_account_json_path
    }
  }

  output {
    Array[String] vet_tables_loaded = PopulateAltAlleleTable.done
  }
}

task GetVetTableNames {
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

    echo "query_project_id = ~{query_project_id}" > ~/.bigqueryrc
    bq query --location=US --query_project_id=~{query_project_id} --format=csv --use_legacy_sql=false \
    "SELECT table_name FROM ~{dataset_project_id}.~{dataset_name}.INFORMATION_SCHEMA.TABLES WHERE table_name LIKE 'vet_%' ORDER BY table_name" > vet_tables.csv

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

task PopulateAltAlleleTable {
  input {
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
    docker: "us.gcr.io/broad-dsde-methods/variantstore:ah_var_store_20210818"
    memory: "3 GB"
    disks: "local-disk 10 HDD"
    preemptible: 3
    cpu: 1
  }
}
