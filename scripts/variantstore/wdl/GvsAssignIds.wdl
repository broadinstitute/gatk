version 1.0

workflow GvsAssignIds {

  input {
    Array[String] external_sample_names
    String project_id
    String dataset_name
    String sample_info_table = "sample_info"
    String sample_info_schema = "sample_name:STRING,sample_id:INTEGER"
    String? service_account_json
    # String? drop_state = "SIXTY"
    Boolean force = false

    Int? preemptible_tries
    File? gatk_override
    String? docker
  }

  String docker_final = select_first([docker, "us.gcr.io/broad-gatk/gatk:4.1.7.0"])


  call SetLock {
    input:
      project_id = project_id,
      dataset_name = dataset_name,
      service_account_json = service_account_json,
      preemptible_tries = preemptible_tries
  }

  call CreateTables as CreateSampleInfoTables {
  	input:
      project_id = project_id,
      dataset_name = dataset_name,
      datatype = "sample_info",
      schema = sample_info_schema,
      service_account_json = service_account_json,
      preemptible_tries = preemptible_tries,
      docker = docker_final
  }

  call AddSamplesToBQ {
    input:
      sample_names = external_sample_names,
      project_id = project_id,
      dataset_name = dataset_name,
      sample_info_table = sample_info_table,
      force = force,
      service_account_json = service_account_json,
      gatk_override = gatk_override,
      table_creation_done = CreateSampleInfoTables.done,
      docker = docker_final,
      preemptible_tries = preemptible_tries
  }

  call AssignIds {
    input:
      service_account_json = service_account_json,
      project_id = project_id,
      dataset_name = dataset_name,
      table_name = sample_info_table,
      gatk_override = gatk_override,
      sample_info_update_done = AddSamplesToBQ.done,
      docker = docker_final,
      preemptible_tries = preemptible_tries

  }
  
  call ReleaseLock {
    input:
      project_id = project_id,
      dataset_name = dataset_name,
      assign_done = AssignIds.done,
      service_account_json = service_account_json,
      preemptible_tries = preemptible_tries
  }

  output {
    Boolean gvs_ids_created = true
  }
}

# we create a bq table as a lock
task SetLock {
  meta {
    volatile: true
  }

  input {
    String project_id
    String dataset_name
    String? service_account_json

    # runtime
    Int? preemptible_tries
  }

  String has_service_account_file = if (defined(service_account_json)) then 'true' else 'false'

  command <<<
    set -e

    if [ ~{has_service_account_file} = 'true' ]; then
      gsutil cp ~{service_account_json} local.service_account.json
      gcloud auth activate-service-account --key-file=local.service_account.json
    fi

    # create the lock table
    bq --project_id=~{project_id} mk ~{dataset_name}.metadata_lock

  >>>

  runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
    memory: "1 GB"
    disks: "local-disk 10 HDD"
    preemptible: select_first([preemptible_tries, 5])
    cpu: 1
  }
}

task ReleaseLock {
  meta {
    volatile: true
  }

  input {
    String project_id
    String dataset_name
    String assign_done
    String? service_account_json

    # runtime
    Int? preemptible_tries
  }

  String has_service_account_file = if (defined(service_account_json)) then 'true' else 'false'

  command <<<
    set -e

    if [ ~{has_service_account_file} = 'true' ]; then
      gsutil cp ~{service_account_json} local.service_account.json
      gcloud auth activate-service-account --key-file=local.service_account.json
    fi

    # remove the lock table
    bq --project_id=~{project_id} rm -f -t ~{dataset_name}.metadata_lock

  >>>

    runtime {
      docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
      memory: "1 GB"
      disks: "local-disk 10 HDD"
      preemptible: select_first([preemptible_tries, 5])
      cpu: 1
    }
}

task AddSamplesToBQ {
  input {
    Array[String] sample_names
    String project_id
    String dataset_name
    String sample_info_table
    String? service_account_json
    String table_creation_done
    Boolean force
    # runtime
    Int? preemptible_tries
    File? gatk_override
    String docker

  }

  String has_service_account_file = if (defined(service_account_json)) then 'true' else 'false'

  meta {
    description: "Adds new sample names to the BQ table"
    volatile: true
  }

  command <<<
      set -ex

      export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

      if [ ~{has_service_account_file} = 'true' ]; then
        gsutil cp ~{service_account_json} local.service_account.json
        gcloud auth activate-service-account --key-file=local.service_account.json
      fi

      NAMES_FILE=~{write_lines(sample_names)}

      echo "project_id = ~{project_id}" > ~/.bigqueryrc

      # first load name into the lock table - will check for dupes before adding to sample_info table
      bq load --project_id=~{project_id} ~{dataset_name}.metadata_lock $NAMES_FILE "sample_name:STRING"

      bq --project_id=~{project_id} query --use_legacy_sql=false --format=csv \
        'SELECT * FROM ~{dataset_name}.~{sample_info_table} s where s.sample_name in (SELECT sample_name FROM ~{dataset_name}.metadata_lock) and s.sample_id is not null' > dupe_sample_names
      ls -l
      wc dupe_sample_names
      sed -i '/^$/d' dupe_sample_names
      wc dupe_sample_names
      if [ -s dupe_sample_names ]; then
        if [ ~{force} = 'true' ]; then
          echo "Ignoring sample names already in sample_info"
          bq --project_id=~{project_id} query --use_legacy_sql=false --format=csv \
            'DELETE ~{dataset_name}.metadata_lock m where m.sample_name in (SELECT sample_name FROM ~{dataset_name}.~{sample_info_table})' 
        else
          echo "Duplicate samples names detected. Aborting run."
          cat dupe_sample_names
          exit 1
        fi
      fi

     bq --project_id=~{project_id} query --use_legacy_sql=false 'INSERT into ~{dataset_name}.~{sample_info_table} (sample_name) select sample_name from ~{dataset_name}.metadata_lock'

      # bq --project_id=~{project_id} ~{dataset_name}.~{sample_info_table} $NAMES_FILE "sample_name:STRING"

  >>>
  runtime {
      docker: docker
      memory: "3.75 GB"
      disks: "local-disk " + 10 + " HDD"
      preemptible: select_first([preemptible_tries, 5])
      cpu: 1
  }
  output {
      Boolean done = true
  }
}

task AssignIds {
  input {
    String project_id
    String dataset_name
    String table_name
    String? service_account_json
    String sample_info_update_done
    # runtime
    Int? preemptible_tries
    File? gatk_override
    String docker

  }

  String has_service_account_file = if (defined(service_account_json)) then 'true' else 'false'

  meta {
    description: "Assigns Ids to samples"
    volatile: true
  }

  command <<<
      set -ex

      export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

      if [ ~{has_service_account_file} = 'true' ]; then
        gsutil cp ~{service_account_json} local.service_account.json
        gcloud auth activate-service-account --key-file=local.service_account.json
      fi

      echo "project_id = ~{project_id}" > ~/.bigqueryrc

      # get the current maximum id, or 0 if there are none
      bq --project_id=~{project_id} query --format=csv --use_legacy_sql=false "SELECT IFNULL(MAX(sample_id),0) FROM ~{dataset_name}.~{table_name}" > maxid
      offset=$(tail -1 maxid)

      # perform actual id assignment
      bq --project_id=~{project_id} query --format=csv --use_legacy_sql=false \
        "UPDATE ~{dataset_name}.~{table_name} m SET m.sample_id = id_assign.id FROM (SELECT sample_name, $offset + ROW_NUMBER() OVER() as id FROM ~{dataset_name}.~{table_name} WHERE sample_id IS NULL) id_assign WHERE m.sample_name = id_assign.sample_name;"
      
  >>>
  runtime {
      docker: docker
      memory: "3.75 GB"
      disks: "local-disk " + 10 + " HDD"
      preemptible: select_first([preemptible_tries, 5])
      cpu: 1
  }
  output {
      Boolean done = true
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
      String? schema
      String? service_account_json

      # runtime
      Int? preemptible_tries
      String docker
    }

    String has_service_account_file = if (defined(service_account_json)) then 'true' else 'false'

  command <<<
    set -x
    set -e

    if [ ~{has_service_account_file} = 'true' ]; then
      gsutil cp ~{service_account_json} local.service_account.json
      gcloud auth activate-service-account --key-file=local.service_account.json
      gcloud config set project ~{project_id}
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

