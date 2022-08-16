version 1.0

import "GvsCreateTables.wdl" as GvsCreateTables

workflow GvsAssignIds {

  input {
    Boolean go = true
    String dataset_name
    String project_id

    Array[String] external_sample_names
    Boolean samples_are_controls = false

    File? assign_ids_gatk_override
  }

  String sample_info_table = "sample_info"
  String sample_info_schema_json = '[{"name": "sample_name","type": "STRING","mode": "REQUIRED"},{"name": "sample_id","type": "INTEGER","mode": "NULLABLE"},{"name":"is_loaded","type":"BOOLEAN","mode":"NULLABLE"},{"name":"is_control","type":"BOOLEAN","mode":"REQUIRED"},{"name":"withdrawn","type":"TIMESTAMP","mode":"NULLABLE"}]'
  String sample_load_status_json = '[{"name": "sample_id","type": "INTEGER","mode": "REQUIRED"},{"name":"status","type":"STRING","mode":"REQUIRED"}, {"name":"event_timestamp","type":"TIMESTAMP","mode":"REQUIRED"}]'

  call GvsCreateTables.CreateTables as CreateSampleInfoTable {
  	input:
      project_id = project_id,
      dataset_name = dataset_name,
      datatype = "sample_info",
      schema_json = sample_info_schema_json,
      max_table_id = 1,
      superpartitioned = "false",
      partitioned = "false"
  }

  call GvsCreateTables.CreateTables as CreateSampleLoadStatusTable {
    input:
      project_id = project_id,
      dataset_name = dataset_name,
      datatype = "sample_load_status",
      schema_json = sample_load_status_json,
      max_table_id = 1,
      superpartitioned = "false",
      partitioned = "false"
  }

  call CreateCostObservabilityTable {
    input:
      project_id = project_id,
      dataset_name = dataset_name,
  }

  call AssignIds {
    input:
      sample_names = external_sample_names,
      project_id = project_id,
      dataset_name = dataset_name,
      sample_info_table = sample_info_table,
      samples_are_controls = samples_are_controls,
      table_creation_done = CreateSampleInfoTable.done,
      gatk_override = assign_ids_gatk_override
  }

  call GvsCreateTables.CreateBQTables as CreateTablesForMaxId {
    input:
      project_id = project_id,
      dataset_name = dataset_name,
      max_table_id = AssignIds.max_table_id
  }

  output {
    Boolean done = true
    File gvs_ids_tsv = AssignIds.gvs_ids_tsv
  }
}

task AssignIds {
  input {
    String dataset_name
    String project_id

    String sample_info_table
    Array[String] sample_names
    Boolean samples_are_controls
    String table_creation_done

    # runtime
    File? gatk_override
  }
  meta {
    description: "Assigns Ids to samples"
    # Not `volatile: true` since successfully assigned IDs should not be assigned again.
  }

  Int samples_per_table = 4000
  # add labels for DSP Cloud Cost Control Labeling and Reporting
  String bq_labels = "--label service:gvs --label team:variants --label managedby:assign_ids"

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
    bq --project_id=~{project_id} query --use_legacy_sql=false ~{bq_labels} \
      'INSERT into `~{dataset_name}.~{sample_info_table}` (sample_name, is_control) select sample_name, ~{samples_are_controls} from `~{dataset_name}.sample_id_assignment_lock` m where m.sample_name not in (SELECT sample_name FROM `~{dataset_name}.~{sample_info_table}`)'

    # get the current maximum id, or 0 if there are none
    bq --project_id=~{project_id} query --format=csv --use_legacy_sql=false ~{bq_labels} 'SELECT IFNULL(MAX(sample_id),0) FROM `~{dataset_name}.~{sample_info_table}`' > maxid
    offset=$(tail -1 maxid)

    # perform actual id assignment
    bq --project_id=~{project_id} query --format=csv --use_legacy_sql=false ~{bq_labels} --parameter=offset:INTEGER:$offset \
      'UPDATE `~{dataset_name}.~{sample_info_table}` m SET m.sample_id = id_assign.id FROM (SELECT sample_name, @offset + ROW_NUMBER() OVER(order by sample_name) as id FROM `~{dataset_name}.~{sample_info_table}` WHERE sample_id IS NULL) id_assign WHERE m.sample_name = id_assign.sample_name;'

    # retrieve the list of assigned ids and samples to update the datamodel
    echo "entity:sample_id,gvs_id" > update.tsv
    bq --project_id=~{project_id} query --format=csv --use_legacy_sql=false ~{bq_labels} -n $num_samples --parameter=offset:INTEGER:$offset \
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

task CreateCostObservabilityTable {
  input {
    String project_id
    String dataset_name
  }

  String cost_observability_json =
                                 '[ { "name": "call_set_identifier", "type": "STRING", "mode": "REQUIRED", "description": "The name by which we refer to the callset" }, ' +
                                 '  { "name": "step", "type": "STRING", "mode": "REQUIRED", "description": "The name of the core GVS workflow to which this belongs" }, ' +
                                 '  { "name": "call", "type": "STRING", "mode": "NULLABLE", "description": "The WDL call to which this belongs" }, ' +
                                 '  { "name": "shard_identifier", "type": "STRING", "mode": "NULLABLE", "description": "A unique identifier for this shard, may or may not be its index" }, ' +
                                 '  { "name": "call_start_timestamp", "type": "TIMESTAMP", "mode": "REQUIRED", "description": "When the call logging this event was started" }, ' +
                                 '  { "name": "event_timestamp", "type": "TIMESTAMP", "mode": "REQUIRED", "description": "When the observability event was logged" }, ' +
                                 '  { "name": "event_key", "type": "STRING", "mode": "REQUIRED", "description": "The type of observability event being logged" }, ' +
                                 '  { "name": "event_bytes", "type": "INTEGER", "mode": "REQUIRED", "description": "Number of bytes reported for this observability event" } ] '

  meta {
    # not volatile: true, always run this when asked
  }
  command <<<
    set -o xtrace

    echo "project_id = ~{project_id}" > ~/.bigqueryrc
    TABLE="~{dataset_name}.cost_observability"

    # Check that the table has not been created yet
    set +o errexit
    bq show --project_id ~{project_id} $TABLE > /dev/null
    BQ_SHOW_RC=$?
    set -o errexit

    if [ $BQ_SHOW_RC -ne 0 ]; then
      PARTITION_STRING="--time_partitioning_field call_start_timestamp --time_partitioning_type DAY"
      echo "making table $TABLE"
      echo '~{cost_observability_json}' > schema.json
      bq --location=US mk ${PARTITION_STRING} --project_id=~{project_id} $TABLE schema.json
    fi
  >>>
  runtime {
    docker: "us.gcr.io/broad-dsde-methods/variantstore:ah_var_store_2022_08_16"
  }
  output {
    Boolean done = true
  }
}

