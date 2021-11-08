version 1.0

workflow GvsAssignIds {

  input {
    Array[String] external_sample_names
    String project_id
    String dataset_name
    String sample_info_table = "sample_info"
    String sample_info_schema = "sample_name:STRING,sample_id:INTEGER"
    String workspace_namespace
    String workspace_name
    String? service_account_json
    Int? preemptible_tries
    File? gatk_override
    String? docker
  }

  String docker_final = select_first([docker, "us.gcr.io/broad-gatk/gatk:4.1.7.0"])


  call CreateTables as CreateSampleInfoTable {
  	input:
      project_id = project_id,
      dataset_name = dataset_name,
      datatype = "sample_info",
      schema = sample_info_schema,
      service_account_json = service_account_json,
      preemptible_tries = preemptible_tries,
      docker = docker_final
  }

  call AssignIds {
    input:
      sample_names = external_sample_names,
      project_id = project_id,
      dataset_name = dataset_name,
      sample_info_table = sample_info_table,
      table_creation_done = CreateSampleInfoTable.done,
      gatk_override = gatk_override,
      service_account_json = service_account_json,
      docker = docker_final,
  }

  call UpdateDataModel {
    input:
      workspace_namespace = workspace_namespace,
      workspace_name = workspace_name,
      gvs_ids_tsv = AssignIds.gvs_ids_tsv
  }


  output {
    Boolean gvs_ids_created = true
    File gvs_ids_tsv = AssignIds.gvs_ids_tsv
  }
}

task AssignIds {
  input {
    Array[String] sample_names
    String project_id
    String dataset_name
    String sample_info_table
    String? service_account_json
    String table_creation_done
    # runtime
    File? gatk_override
    String docker

  }

  String has_service_account_file = if (defined(service_account_json)) then 'true' else 'false'

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
        gsutil cp ~{service_account_json} local.service_account.json
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
        'INSERT into ~{dataset_name}.~{sample_info_table} (sample_name) select sample_name from ~{dataset_name}.sample_id_assignment_lock m where m.sample_name not in (SELECT sample_name FROM ~{dataset_name}.~{sample_info_table})'

      # get the current maximum id, or 0 if there are none
      bq --project_id=~{project_id} query --format=csv --use_legacy_sql=false "SELECT IFNULL(MAX(sample_id),0) FROM ~{dataset_name}.~{sample_info_table}" > maxid
      offset=$(tail -1 maxid)

      # perform actual id assignment
      bq --project_id=~{project_id} query --format=csv --use_legacy_sql=false \
        "UPDATE ~{dataset_name}.~{sample_info_table} m SET m.sample_id = id_assign.id FROM (SELECT sample_name, $offset + ROW_NUMBER() OVER() as id FROM ~{dataset_name}.~{sample_info_table} WHERE sample_id IS NULL) id_assign WHERE m.sample_name = id_assign.sample_name;"

      # retrieve the list of assigned ids and samples to update the datamodel
      echo "entity:sample_id,gvs_id" > update.tsv
      bq --project_id=~{project_id} query --format=csv --use_legacy_sql=false -n $num_samples \
        "SELECT sample_name, sample_id from ~{dataset_name}.~{sample_info_table} WHERE sample_id >= $offset" > update.tsv
      cat update.tsv | sed -e 's/sample_id/gvs_id/' -e 's/sample_name/entity:sample_id/' -e 's/,/\t/g' > gvs_ids.tsv

      # remove the lock table
      bq --project_id=~{project_id} rm -f -t ~{dataset_name}.sample_id_assignment_lock

  >>>
  runtime {
      docker: docker
      memory: "3.75 GB"
      disks: "local-disk " + 10 + " HDD"
      cpu: 1
  }
  output {
    File gvs_ids_tsv = "gvs_ids.tsv"
  }
}

task UpdateDataModel {
  input {
    String workspace_namespace
    String workspace_name
    File gvs_ids_tsv
  }
  meta {
    description: "Assigns gvs_id attribute in data model"
    volatile: true
  }
  command <<<
      # update the data model
      python3 <<CODE
from firecloud import api as fapi
response = fapi.upload_entities_tsv('~{workspace_namespace}', '~{workspace_name}', '~{gvs_ids_tsv}', 'flexible')
if response.status_code != 200:
  print(response.status_code)
  print(response.text)
  exit 1
CODE
  >>>
  runtime {
      docker: "broadgdac/fiss"
      memory: "3.75 GB"
      disks: "local-disk " + 10 + " HDD"
      cpu: 1
  }
}

# Creates the necessary sample_info table
task CreateTables {
  meta {
    volatile: true
  }

	input {
      String project_id
      String dataset_name
      String datatype
      String? schema
      String? service_account_json

      # runtime
      Int? preemptible_tries
      String docker
    }

    String has_service_account_file = if (defined(service_account_json)) then 'true' else 'false'

  command <<<
    set -e

    if [ ~{has_service_account_file} = 'true' ]; then
      gsutil cp ~{service_account_json} local.service_account.json
      gcloud auth activate-service-account --key-file=local.service_account.json
    fi

    echo "project_id = ~{project_id}" > ~/.bigqueryrc

    TABLE="~{dataset_name}.~{datatype}"

    # Check that the table has not been created yet
    set +e
    bq show --project_id ~{project_id} $TABLE > /dev/null
    BQ_SHOW_RC=$?
    set -e
    if [ $BQ_SHOW_RC -ne 0 ]; then
      echo "making table $TABLE"
      bq --location=US mk --project_id=~{project_id} $TABLE ~{schema}
    fi
  >>>

  output {
    String done = "true"
  }

  runtime {
    docker: docker
    memory: "3 GB"
    disks: "local-disk 10 HDD"
    preemptible: select_first([preemptible_tries, 5])
    cpu: 1
  }
}

