version 1.0

import "GvsCreateTables.wdl" as GvsCreateTables

workflow GvsAssignIds {

  input {
    String dataset_name
    String project_id

    Array[String] external_sample_names

    File? assign_ids_gatk_override
    Int? create_tables_preemptible_override
    String? service_account_json_path
  }

  String sample_info_table = "sample_info"
  String sample_info_schema_json = '[{"name": "sample_name","type": "STRING","mode": "REQUIRED"},{"name": "sample_id","type": "INTEGER","mode": "NULLABLE"},{"name":"is_loaded","type":"BOOLEAN","mode":"NULLABLE"}]'
  String sample_load_status_json = '[{"name": "sample_id","type": "INTEGER","mode": "REQUIRED"},{"name":"status","type":"STRING","mode":"REQUIRED"}, {"name":"event_timestamp","type":"TIMESTAMP","mode":"REQUIRED"}]'


  call GvsCreateTables.CreateTables as CreateSampleInfoTable {
  	input:
      project_id = project_id,
      dataset_name = dataset_name,
      datatype = "sample_info",
      schema_json = sample_info_schema_json,
      max_table_id = 1,
      superpartitioned = "false",
      partitioned = "false",
      service_account_json_path = service_account_json_path,
      preemptible_tries = create_tables_preemptible_override,
  }

  call GvsCreateTables.CreateTables as CreateSampleLoadStatusTable {
    input:
      project_id = project_id,
      dataset_name = dataset_name,
      datatype = "sample_load_status",
      schema_json = sample_load_status_json,
      max_table_id = 1,
      superpartitioned = "false",
      partitioned = "false",
      service_account_json_path = service_account_json_path,
      preemptible_tries = create_tables_preemptible_override,
  }

  call AssignIds {
    input:
      sample_names = external_sample_names,
      project_id = project_id,
      dataset_name = dataset_name,
      sample_info_table = sample_info_table,
      table_creation_done = CreateSampleInfoTable.done,
      gatk_override = assign_ids_gatk_override,
      service_account_json_path = service_account_json_path
  }

  call GvsCreateTables.CreateBQTables as CreateTablesForMaxId {
    input:
      project_id = project_id,
      dataset_name = dataset_name,
      max_table_id = AssignIds.max_table_id,
      service_account_json_path = service_account_json_path,
      preemptible_tries = create_tables_preemptible_override
  }

  output {
    Boolean gvs_ids_created = true
    File gvs_ids_tsv = AssignIds.gvs_ids_tsv
  }
}

task AssignIds {
  input {
    String dataset_name
    String project_id

    String sample_info_table
    Array[String] sample_names
    String table_creation_done

    # runtime
    File? gatk_override
    String? service_account_json_path
  }

  Int samples_per_table = 4000
  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

  meta {
    description: "Assigns Ids to samples"
    volatile: true
  }

  command <<<
    set -e
    set -x

    # make sure that sample names were actually passed, fail if empty
    num_samples=~{length(sample_names)}
    if [ $num_samples -eq 0 ]; then
      echo "No sample names passed. Exiting"
      exit 1
    fi

    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

    if [ ~{has_service_account_file} = 'true' ]; then
      gsutil cp ~{service_account_json_path} local.service_account.json
      gcloud auth activate-service-account --key-file=local.service_account.json
    fi

    echo "project_id = ~{project_id}" > ~/.bigqueryrc

    # Check that lock table has not been created yet
    set +e
    bq --project_id=~{project_id} show ~{dataset_name}.sample_id_assignment_lock > /dev/null
    BQ_SHOW_RC=$?
    set -e
    if [ $BQ_SHOW_RC -eq 0 ]; then
      echo "lock table already exists. Exiting"
      exit 1
    fi

    # create the lock table
    bq --project_id=~{project_id} mk ~{dataset_name}.sample_id_assignment_lock "sample_name:STRING"

    NAMES_FILE=~{write_lines(sample_names)}

    # first load name into the lock table - will check for dupes when adding to sample_info table
    bq load --project_id=~{project_id} ~{dataset_name}.sample_id_assignment_lock $NAMES_FILE "sample_name:STRING"

    # add sample_name to sample_info_table
    bq --project_id=~{project_id} query --use_legacy_sql=false \
      'INSERT into `~{dataset_name}.~{sample_info_table}` (sample_name) select sample_name from `~{dataset_name}.sample_id_assignment_lock` m where m.sample_name not in (SELECT sample_name FROM `~{dataset_name}.~{sample_info_table}`)'

    # get the current maximum id, or 0 if there are none
    bq --project_id=~{project_id} query --format=csv --use_legacy_sql=false 'SELECT IFNULL(MAX(sample_id),0) FROM `~{dataset_name}.~{sample_info_table}`' > maxid
    offset=$(tail -1 maxid)

    # perform actual id assignment
    bq --project_id=~{project_id} query --format=csv --use_legacy_sql=false --parameter=offset:INTEGER:$offset \
      'UPDATE `~{dataset_name}.~{sample_info_table}` m SET m.sample_id = id_assign.id FROM (SELECT sample_name, @offset + ROW_NUMBER() OVER() as id FROM `~{dataset_name}.~{sample_info_table}` WHERE sample_id IS NULL) id_assign WHERE m.sample_name = id_assign.sample_name;'

    # retrieve the list of assigned ids and samples to update the datamodel
    echo "entity:sample_id,gvs_id" > update.tsv
    bq --project_id=~{project_id} query --format=csv --use_legacy_sql=false -n $num_samples --parameter=offset:INTEGER:$offset \
      'SELECT sample_name, sample_id from `~{dataset_name}.~{sample_info_table}` WHERE sample_id >= @offset' > update.tsv
    cat update.tsv | sed -e 's/sample_id/gvs_id/' -e 's/sample_name/entity:sample_id/' -e 's/,/\t/g' > gvs_ids.tsv

    # get the max id to create tables for
    max_sample_id=$(cat update.tsv | cut -d, -f2 | sort -r -n | head -1)
    python3 -c "from math import ceil; print(ceil($max_sample_id/~{samples_per_table}))" > max_table_id

    # remove the lock table
    bq --project_id=~{project_id} rm -f -t ~{dataset_name}.sample_id_assignment_lock
  >>>
  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:4.1.7.0"
    memory: "3.75 GB"
    disks: "local-disk " + 10 + " HDD"
    cpu: 1
  }
  output {
    File gvs_ids_tsv = "gvs_ids.tsv"
    Int max_table_id = read_int("max_table_id")
  }
}

