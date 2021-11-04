version 1.0

workflow CreateBQTables {

    input {
        Int max_table_id
        String project_id
        String dataset_name
        String? service_account_json_path
        String pet_schema = "location:INTEGER,sample_id:INTEGER,state:STRING"
        String vet_schema = "sample_id:INTEGER,location:INTEGER,ref:STRING,alt:STRING,AS_RAW_MQ:STRING,AS_RAW_MQRankSum:STRING,QUALapprox:STRING,AS_QUALapprox:STRING,AS_RAW_ReadPosRankSum:STRING,AS_SB_TABLE:STRING,AS_VarDP:STRING,call_GT:STRING,call_AD:STRING,call_GQ:INTEGER,call_PGT:STRING,call_PID:STRING,call_PL:STRING"
        String uuid = ""
        Int? preemptible_tries

    }

    call CreateTables as CreatePetTables {
        input:
        project_id = project_id,
        dataset_name = dataset_name,
        datatype = "pet",
        max_table_id = max_table_id,
        schema = pet_schema,
        superpartitioned = "true",
        partitioned = "true",
        uuid = uuid,
        service_account_json_path = service_account_json_path,
        preemptible_tries = preemptible_tries
    }

    call CreateTables as CreateVetTables {
        input:
        project_id = project_id,
        dataset_name = dataset_name,
        datatype = "vet",
        max_table_id = max_table_id,
        schema = vet_schema,
        superpartitioned = "true",
        partitioned = "true",
        uuid = uuid,
        service_account_json_path = service_account_json_path,
        preemptible_tries = preemptible_tries
    }
}


# Creates all the tables necessary for the LoadData operation
task CreateTables {
  meta {
    volatile: true
  }

	input {
      String project_id
      String dataset_name
      String datatype
      Int max_table_id
      String? schema
      String superpartitioned
      String partitioned
      String uuid
      String? service_account_json_path

      # runtime
      Int? preemptible_tries
    }

    String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

  command <<<
    set -x
    set -e

    if [ ~{has_service_account_file} = 'true' ]; then
      gsutil cp ~{service_account_json_path} local.service_account.json
      export GOOGLE_APPLICATION_CREDENTIALS=local.service_account.json
      gcloud auth activate-service-account --key-file=local.service_account.json
      gcloud config set project ~{project_id}
    fi

    echo "project_id = ~{project_id}" > ~/.bigqueryrc

    PREFIX=""
    if [ -n "~{uuid}" ]; then
      PREFIX="~{uuid}_"
    fi

    for TABLE_ID in $(seq 1 ~{max_table_id}); do
      PARTITION_STRING=""
      CLUSTERING_STRING=""
      if [ ~{partitioned} == "true" ]; then
        # assume clustering as well
        let "PARTITION_START=(${TABLE_ID}-1)*4000+1"
        let "PARTITION_END=$PARTITION_START+4000"
        let "PARTITION_STEP=1"
        PARTITION_FIELD="sample_id"
        CLUSTERING_FIELD="location"
        PARTITION_STRING="--range_partitioning=$PARTITION_FIELD,$PARTITION_START,$PARTITION_END,$PARTITION_STEP"
        CLUSTERING_STRING="--clustering_fields=$CLUSTERING_FIELD"
      fi

      if [ ~{superpartitioned} = "true" ]; then
        printf -v PADDED_TABLE_ID "%03d" ${TABLE_ID}
        TABLE="~{dataset_name}.${PREFIX}~{datatype}_${PADDED_TABLE_ID}"
      else
        TABLE="~{dataset_name}.${PREFIX}~{datatype}"
      fi

      # Check that the table has not been created yet
      set +e
      bq show --project_id ~{project_id} $TABLE > /dev/null
      BQ_SHOW_RC=$?
      set -e
      if [ $BQ_SHOW_RC -ne 0 ]; then
        echo "making table $TABLE"
        bq --location=US mk ${PARTITION_STRING} ${CLUSTERING_STRING} --project_id=~{project_id} $TABLE ~{schema}
      fi
    done
  >>>

  output {
    String done = "true"
  }

  runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:latest"
    memory: "3 GB"
    disks: "local-disk 10 HDD"
    preemptible: select_first([preemptible_tries, 5])
    cpu: 1
  }
}
