version 1.0

workflow CreateBQTables {
  input {
    String dataset_name
    String project_id

    Int max_table_id

    Int? preemptible_tries
  }

  String pet_schema_json = '[{"name": "location","type": "INTEGER","mode": "REQUIRED"},{"name": "sample_id","type": "INTEGER","mode": "REQUIRED"},{"name": "state","type": "STRING","mode": "REQUIRED"}]'
  String ref_ranges_schema_json = '[{"name": "location","type": "INTEGER","mode": "REQUIRED"},{"name": "sample_id","type": "INTEGER","mode": "REQUIRED"},{"name": "length","type": "INTEGER","mode": "REQUIRED"},{"name": "state","type": "STRING","mode": "REQUIRED"}]'
  String vet_schema_json = '[{"name": "sample_id", "type" :"INTEGER", "mode": "REQUIRED"},{"name": "location", "type" :"INTEGER", "mode": "REQUIRED"},{"name": "ref", "type" :"STRING", "mode": "REQUIRED"},{"name": "alt", "type" :"STRING", "mode": "REQUIRED"},{"name": "AS_RAW_MQ", "type" :"STRING", "mode": "NULLABLE"},{"name": "AS_RAW_MQRankSum", "type" :"STRING", "mode": "NULLABLE"},{"name": "QUALapprox", "type" :"STRING", "mode": "NULLABLE"},{"name": "AS_QUALapprox", "type" :"STRING", "mode": "NULLABLE"},{"name": "AS_RAW_ReadPosRankSum", "type" :"STRING", "mode": "NULLABLE"},{"name": "AS_SB_TABLE", "type" :"STRING", "mode": "NULLABLE"},{"name": "AS_VarDP", "type" :"STRING", "mode": "NULLABLE"},{"name": "call_GT", "type" :"STRING", "mode": "NULLABLE"},{"name": "call_AD", "type" :"STRING", "mode": "NULLABLE"},{"name": "call_GQ", "type" :"INTEGER", "mode": "NULLABLE"},{"name": "call_PGT", "type" :"STRING", "mode": "NULLABLE"},{"name": "call_PID", "type" :"STRING", "mode": "NULLABLE"},{"name": "call_PL", "type" :"STRING", "mode": "NULLABLE"}]'

  call CreateTables as CreateVetTables {
    input:
      project_id = project_id,
      dataset_name = dataset_name,
      datatype = "vet",
      max_table_id = max_table_id,
      schema_json = vet_schema_json,
      superpartitioned = "true",
      partitioned = "true"
  }

  call CreateTables as CreateRefRangesTables {
    input:
      project_id = project_id,
      dataset_name = dataset_name,
      datatype = "ref_ranges",
      max_table_id = max_table_id,
      schema_json = ref_ranges_schema_json,
      superpartitioned = "true",
      partitioned = "true"
  }

  output {
    String vetDone = CreateVetTables.done
    String refDone = CreateRefRangesTables.done
  }
}


# Creates all the tables necessary for the LoadData operation
task CreateTables {
  input {
    String project_id
    String dataset_name
    String datatype
    Int max_table_id
    String schema_json
    String superpartitioned
    String partitioned
  }
  meta {
    # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
  }

  command <<<
    set -x
    set -e

    echo "project_id = ~{project_id}" > ~/.bigqueryrc

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
        TABLE="~{dataset_name}.~{datatype}_${PADDED_TABLE_ID}"
      else
        TABLE="~{dataset_name}.~{datatype}"
      fi

      # Check that the table has not been created yet
      set +e
      bq show --project_id ~{project_id} $TABLE > /dev/null
      BQ_SHOW_RC=$?
      set -e
      if [ $BQ_SHOW_RC -ne 0 ]; then
        echo "making table $TABLE"
        echo '~{schema_json}' > schema.json
        bq --location=US mk ${PARTITION_STRING} ${CLUSTERING_STRING} --project_id=~{project_id} $TABLE schema.json
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
    cpu: 1
  }
}
