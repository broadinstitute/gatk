version 1.0

workflow GvsWithdrawSamples {

  input {
    String dataset_name
    String project_id

    Array[String] sample_names

    String? service_account_json_path
  }

  call WithdrawSamples {
    input:
      project_id = project_id,
      dataset_name = dataset_name,
      sample_names = sample_names,
      service_account_json_path = service_account_json_path
  }

  output {
    Int num_rows_updated = WithdrawSamples.num_rows_updated
  }
}

task WithdrawSamples {
  input {
    String project_id
    String dataset_name

    Array[String] sample_names

    String? service_account_json_path
  }

  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

  meta {
    description: "Withdraw Samples from GVS by marking them as 'withdrawn' in the sample_info table"
    volatile: true
  }

  command <<<
    set -e
    set -x

    # make sure that sample names were actually passed, warn and exit if empty
    num_samples=~{length(sample_names)}
    if [ $num_samples -eq 0 ]; then
      echo "No sample names passed. Exiting"
      exit 0
    fi

    if [ ~{has_service_account_file} = 'true' ]; then
      gsutil cp ~{service_account_json_path} local.service_account.json
      gcloud auth activate-service-account --key-file=local.service_account.json
    fi

    echo "project_id = ~{project_id}" > ~/.bigqueryrc

    # perform actual update
    bq --project_id=~{project_id} query --format=csv --use_legacy_sql=false \
      'UPDATE `~{dataset_name}.sample_info` SET withdrawn = CURRENT_TIMESTAMP() WHERE sample_name IN ("~{sep='\", \"' sample_names}")' > log_message.txt;

    cat log_message.txt | sed -e 's/Number of affected rows: //' > rows_updated.txt
    typeset -i rows_updated=$(cat rows_updated.txt)

    if [ $num_samples -ne $rows_updated ]; then
      echo "Error: Expected to update $num_samples rows - but only updated $rows_updated."
      exit 1
    fi

  >>>
  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:4.2.6.0"
    memory: "3.75 GB"
    disks: "local-disk 10 HDD"
    cpu: 1
  }
  output {
    Int num_rows_updated = read_int("rows_updated.txt")
  }
}

