version 1.0

workflow GvsWithdrawSamples {

  input {
    String dataset_name
    String project_id

    Array[String] sample_names

    File? withdraw_samples_gatk_override
    String? service_account_json_path
  }

  String sample_info_table = "sample_info"

  call WithdrawSamples {
    input:
      project_id = project_id,
      dataset_name = dataset_name,
      sample_info_table = sample_info_table,
      sample_names = sample_names,
      gatk_override = withdraw_samples_gatk_override,
      service_account_json_path = service_account_json_path
  }

#  output {
#    Boolean gvs_ids_created = true
#    File gvs_ids_tsv = WithdrawSamples.gvs_ids_tsv
#  }
}

task WithdrawSamples {
  input {
    String project_id
    String dataset_name

    String sample_info_table
    Array[String] sample_names

    # runtime
    File? gatk_override
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

    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

    ## TODO - do I need to do this.
    if [ ~{has_service_account_file} = 'true' ]; then
      gsutil cp ~{service_account_json_path} local.service_account.json
      gcloud auth activate-service-account --key-file=local.service_account.json
    fi

    echo "project_id = ~{project_id}" > ~/.bigqueryrc

    NAMES_FILE=~{write_lines(sample_names)}

    echo "Hello!"
    echo ~{sep='"\', \'"' sample_names}
    echo "There!"

    # perform actual update
    ## TODO - retrieve the number updated and verify that it matches the number of samples?
    bq --project_id=~{project_id} query --format=csv --use_legacy_sql=false \
    'UPDATE `~{dataset_name}.~{sample_info_table}` SET withdrawn = CURRENT_TIMESTAMP() WHERE sample_name IN (~{sep='"\', \'"' sample_names})';

  >>>
  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:4.1.7.0"
    memory: "3.75 GB"
    disks: "local-disk " + 10 + " HDD"
    cpu: 1
  }
#  output {
#    File gvs_ids_tsv = "gvs_ids.tsv"
#    Int max_table_id = read_int("max_table_id")
#  }
}

