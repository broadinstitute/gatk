version 1.0

workflow ImportGenomes {

  input {
    Array[File] input_vcfs
    Array[File] input_vcf_indexes
    File interval_list
    String output_directory
    File sample_map
    String project_id
    String dataset_name
    String? pet_schema
    String? vet_schema
    String? metadata_schema
    File? service_account_json
    String? drop_state
    Boolean? drop_state_includes_greater_than = false

    Int? preemptible_tries
    File? gatk_override
    String? docker
  }

  String docker_final = select_first([docker, "us.gcr.io/broad-gatk/gatk:4.1.7.0"])
  String pet_inline_schema = if defined(pet_schema) then "$pet_schema" else "location:INTEGER,sample_id:INTEGER,state:STRING"
  String vet_inline_schema = if defined(vet_schema) then "$vet_schema" else "sample_id:INTEGER,location:INTEGER, \
                                                                             ref:STRING,alt:STRING,AS_RAW_MQ:STRING, \
                                                                             AS_RAW_MQRankSum:STRING,AS_QUALapprox:STRING, \
                                                                             AS_RAW_ReadPosRankSum:STRING,AS_SB_TABLE:STRING, \
                                                                             AS_VarDP:STRING,call_GT:STRING, \
                                                                             call_AD:STRING,call_GQ:INTEGER, \
                                                                             call_PGT:STRING,call_PID:STRING,call_PL:STRING"
  String metadata_inline_schema = if defined(metadata_schema) then "$metadata_schema" else "sample_name:STRING,sample_id:INTEGER, \
                                                                                            interval_list_blob:STRING,inferred_state:STRING"

  call SetLock {
    input:
      output_directory = output_directory,
      service_account_json = service_account_json,
      preemptible_tries = preemptible_tries
  }

  call GetMaxTableId {
    input:
      sample_map = sample_map
  }

  call CreateTables as CreateMetadataTables {
  	input:
      project_id = project_id,
      dataset_name = dataset_name,
      datatype = "metadata",
      max_table_id = GetMaxTableId.max_table_id,
      schema = metadata_inline_schema,
      superpartitioned = "false",
      partitioned = "false",
      uuid = "",
      service_account_json = service_account_json,
      preemptible_tries = preemptible_tries,
      docker = docker_final
  }

  call CreateTables as CreatePetTables {
  	input:
      project_id = project_id,
      dataset_name = dataset_name,
      datatype = "pet",
      max_table_id = GetMaxTableId.max_table_id,
      schema = pet_inline_schema,
      superpartitioned = "true",
      partitioned = "true",
      uuid = "",
      service_account_json = service_account_json,
      preemptible_tries = preemptible_tries,
      docker = docker_final
  }

  call CreateTables as CreateVetTables {
  	input:
      project_id = project_id,
      dataset_name = dataset_name,
      datatype = "vet",
      max_table_id = GetMaxTableId.max_table_id,
      schema = vet_inline_schema,
      superpartitioned = "true",
      partitioned = "true",
      uuid = "",
      service_account_json = service_account_json,
      preemptible_tries = preemptible_tries,
      docker = docker_final
  }

  scatter (i in range(length(input_vcfs))) {
    call CreateImportTsvs {
      input:
        input_vcf = input_vcfs[i],
        input_vcf_index = input_vcf_indexes[i],
        interval_list = interval_list,
        sample_map = sample_map,
        service_account_json = service_account_json,
        drop_state = drop_state,
        drop_state_includes_greater_than = drop_state_includes_greater_than,
        output_directory = output_directory,
        gatk_override = gatk_override,
        docker = docker_final,
        preemptible_tries = preemptible_tries,
        run_uuid = SetLock.run_uuid
    }
  }

  scatter (i in range(GetMaxTableId.max_table_id)) {
    call LoadTable as LoadMetadataTable {
      input:
        project_id = project_id,
        table_id = i + 1,
        dataset_name = dataset_name,
        storage_location = output_directory,
        datatype = "metadata",
        superpartitioned = "false",
        schema = metadata_inline_schema,
        table_creation_done = CreateMetadataTables.done,
        tsv_creation_done = CreateImportTsvs.done,
        service_account_json = service_account_json,
        docker = docker_final,
        run_uuid = SetLock.run_uuid
    }
  }

  scatter (i in range(GetMaxTableId.max_table_id)) {
    call LoadTable as LoadPetTable {
    input:
      project_id = project_id,
      table_id = i + 1,
      dataset_name = dataset_name,
      storage_location = output_directory,
      datatype = "pet",
      superpartitioned = "true",
      schema = pet_inline_schema,
      table_creation_done = CreatePetTables.done,
      tsv_creation_done = CreateImportTsvs.done,
      service_account_json = service_account_json,
      docker = docker_final,
      run_uuid = SetLock.run_uuid
    }
  }

  scatter (i in range(GetMaxTableId.max_table_id)) {
    call LoadTable as LoadVetTable {
    input:
      project_id = project_id,
      table_id = i + 1,
      dataset_name = dataset_name,
      storage_location = output_directory,
      datatype = "vet",
      superpartitioned = "true",
      schema = vet_inline_schema,
      table_creation_done = CreateVetTables.done,
      tsv_creation_done = CreateImportTsvs.done,
      service_account_json = service_account_json,
      docker = docker_final,
      run_uuid = SetLock.run_uuid
    }
  }

  call ReleaseLock {
    input:
      run_uuid = SetLock.run_uuid,
      output_directory = output_directory,
      load_metadata_done = LoadMetadataTable.done,
      load_pet_done = LoadPetTable.done,
      load_vet_done = LoadVetTable.done,
      service_account_json = service_account_json,
      preemptible_tries = preemptible_tries
  }
}

# we create a lock file in the output directory with a uuid for this run of ImportGenomes.
# other tasks (TSV creation, bq load) check that the lock file exists and contains the run_uuid
# specific to this task.
task SetLock {
  meta {
    volatile: true
  }

  input {
    String output_directory
    File? service_account_json

    # runtime
    Int? preemptible_tries
  }

  String has_service_account_file = if (defined(service_account_json)) then 'true' else 'false'

  command <<<
    set -x
    set -e

    if [ ~{has_service_account_file} = 'true' ]; then
      gcloud auth activate-service-account --key-file='~{service_account_json}'
    fi

    # generate uuid for this run
    RUN_UUID=$(dbus-uuidgen)
    echo $RUN_UUID | tee RUN_UUID_STRING

    DIR="~{output_directory}/"

    # check for existing lock file
    LOCKFILE="LOCKFILE"
    HAS_LOCKFILE=$(gsutil ls "${DIR}${LOCKFILE}" | wc -l)
    if [ $HAS_LOCKFILE -gt 0 ]; then
      echo "ERROR: lock file in place. Check whether another run of ImportGenomes with this output directory is in progress or a previous run had an error. If you would like to proceed, run `gsutil rm ${DIR}${LOCKFILE}` and re-run the workflow." 1>&2
      exit 1
    else  # put the lock file in place
      echo "Setting lock file with UUID ${RUN_UUID}"
      echo $RUN_UUID > $LOCKFILE
      gsutil cp $LOCKFILE "${DIR}${LOCKFILE}" || { echo "Error uploading lockfile to ${DIR}${LOCKFILE}" 1>&2 ; exit 1; }
    fi
  >>>

  runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
    memory: "1 GB"
    disks: "local-disk 10 HDD"
    preemptible: select_first([preemptible_tries, 5])
    cpu: 1
  }

  output {
    String run_uuid = read_string("RUN_UUID_STRING")
  }
}

task ReleaseLock {
  meta {
    volatile: true
  }

  input {
    String run_uuid
    String output_directory
    Array[String] load_metadata_done
    Array[String] load_pet_done
    Array[String] load_vet_done
    File? service_account_json

    # runtime
    Int? preemptible_tries
  }

  String has_service_account_file = if (defined(service_account_json)) then 'true' else 'false'

  command <<<
    set -x
    set -e

    if [ ~{has_service_account_file} = 'true' ]; then
      gcloud auth activate-service-account --key-file='~{service_account_json}'
    fi


    LOCKFILE="~{output_directory}/LOCKFILE"
    EXISTING_LOCK_ID=$(gsutil cat ${LOCKFILE})
    CURRENT_RUN_ID="~{run_uuid}"

    if [ ${EXISTING_LOCK_ID} = ${CURRENT_RUN_ID} ]; then
      gsutil rm $LOCKFILE
    else
      echo "ERROR: found mismatched lockfile containing run ${EXISTING_LOCK_ID}, which does not match this run ${CURRENT_RUN_ID}." 1>&2
      exit 1
    fi
  >>>

    runtime {
      docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
      memory: "1 GB"
      disks: "local-disk 10 HDD"
      preemptible: select_first([preemptible_tries, 5])
      cpu: 1
    }
}

task GetMaxTableId {
  input {
    File sample_map
    Int? samples_per_table = 4000

    # runtime
    Int? preemptible_tries
  }

  command <<<
      set -e
      max_sample_id=$(cat ~{sample_map} | cut -d"," -f1 | sort -rn | head -1)
      python -c "from math import ceil; print(ceil($max_sample_id/~{samples_per_table}))"
  >>>
  runtime {
      docker: "python:3.8-slim-buster"
      memory: "1 GB"
      disks: "local-disk 10 HDD"
      preemptible: select_first([preemptible_tries, 5])
      cpu: 1
  }
  output {
      Int max_table_id = read_int(stdout())
  }
}

task CreateImportTsvs {
  input {
    File input_vcf
    File input_vcf_index
    File interval_list
    String output_directory
    File sample_map
    File? service_account_json
    String? drop_state
    Boolean? drop_state_includes_greater_than = false

    # runtime
    Int? preemptible_tries
    File? gatk_override
    String docker

    String? for_testing_only
    String run_uuid
  }

  Int disk_size = if defined(drop_state) then 30 else 75
  String has_service_account_file = if (defined(service_account_json)) then 'true' else 'false'
  # if we are doing a manual localization, we need to set the filename
  String updated_input_vcf = if (defined(service_account_json)) then basename(input_vcf) else input_vcf

  meta {
    description: "Creates a tsv file for import into BigQuery"
    volatile: true
  }

  parameter_meta {
    input_vcf: {
      localization_optional: true
    }
    input_vcf_index: {
      localization_optional: true
    }
  }
  command <<<
      set -e

      # workaround for https://github.com/broadinstitute/cromwell/issues/3647
      export TMPDIR=/tmp

      export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}
      ~{for_testing_only}

      if [ ~{has_service_account_file} = 'true' ]; then
        gcloud auth activate-service-account --key-file='~{service_account_json}'
        gsutil cp ~{input_vcf} .
        gsutil cp ~{input_vcf_index} .
      fi

      # check for existence of the correct lockfile
      LOCKFILE="~{output_directory}/LOCKFILE"
      EXISTING_LOCK_ID=$(gsutil cat ${LOCKFILE}) || { echo "Error retrieving lockfile from ${LOCKFILE}" 1>&2 ; exit 1; }
      CURRENT_RUN_ID="~{run_uuid}"

      if [ ${EXISTING_LOCK_ID} != ${CURRENT_RUN_ID} ]; then
        echo "ERROR: found mismatched lockfile containing run ${EXISTING_LOCK_ID}, which does not match this run ${CURRENT_RUN_ID}." 1>&2
        exit 1
      fi
      
      gatk --java-options "-Xmx7000m" CreateVariantIngestFiles \
        -V ~{updated_input_vcf} \
        -L ~{interval_list} \
        ~{"-IG " + drop_state} \
        --ignore-above-gq-threshold ~{drop_state_includes_greater_than} \
        --mode GENOMES \
        -SNM ~{sample_map} \
        --ref-version 38

      gsutil -m cp metadata_*.tsv ~{output_directory}/metadata_tsvs/
      gsutil -m cp pet_*.tsv ~{output_directory}/pet_tsvs/
      gsutil -m cp vet_*.tsv ~{output_directory}/vet_tsvs/
  >>>
  runtime {
      docker: docker
      memory: "3.75 GB"
      disks: "local-disk " + disk_size + " HDD"
      preemptible: select_first([preemptible_tries, 5])
      cpu: 1
  }
  output {
      String done = "true"
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
      String schema
      String superpartitioned
      String partitioned
      String uuid
      File? service_account_json

      # runtime
      Int? preemptible_tries
      String docker
    }

    String has_service_account_file = if (defined(service_account_json)) then 'true' else 'false'

  command <<<
    set -x
    set -e

    if [ ~{has_service_account_file} = 'true' ]; then
      gcloud auth activate-service-account --key-file='~{service_account_json}'
      gcloud config set project ~{project_id}
    fi

    PREFIX=""
    if [ -n "~{uuid}" ]; then
      PREFIX="~{uuid}_"
    fi

    for TABLE_ID in $(seq 1 ~{max_table_id}); do
      PARTITION_STRING=""
      if [ ~{partitioned} == "true" ]; then
        let "PARTITION_START=(${TABLE_ID}-1)*4000+1"
        let "PARTITION_END=$PARTITION_START+3999"
        let "PARTITION_STEP=1"
        PARTITION_FIELD="sample_id"
        PARTITION_STRING="--range_partitioning=$PARTITION_FIELD,$PARTITION_START,$PARTITION_END,$PARTITION_STEP"
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
        bq --location=US mk ${PARTITION_STRING} --project_id=~{project_id} $TABLE ~{schema}
      fi
    done
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

task LoadTable {
  meta {
    volatile: true
  }

  input {
    String project_id
    String table_id
    String dataset_name
    String storage_location
    String datatype
    String superpartitioned
    String schema
    File? service_account_json
    String table_creation_done
    Array[String] tsv_creation_done
    String run_uuid

    String docker
  }

  String has_service_account_file = if (defined(service_account_json)) then 'true' else 'false'

  command <<<
    set -x
    set -e

    if [ ~{has_service_account_file} = 'true' ]; then
      gcloud auth activate-service-account --key-file='~{service_account_json}'
      gcloud config set project ~{project_id}
    fi

    DIR="~{storage_location}/~{datatype}_tsvs/"
    # check for existence of the correct lockfile
    LOCKFILE="~{storage_location}/LOCKFILE"
    EXISTING_LOCK_ID=$(gsutil cat ${LOCKFILE}) || { echo "Error retrieving lockfile from ${LOCKFILE}" 1>&2 ; exit 1; }
    CURRENT_RUN_ID="~{run_uuid}"

    if [ "${EXISTING_LOCK_ID}" != "${CURRENT_RUN_ID}" ]; then
    echo "ERROR: found mismatched lockfile containing run ${EXISTING_LOCK_ID}, which does not match this run ${CURRENT_RUN_ID}." 1>&2
    exit 1
    fi

    DIR="~{storage_location}/~{datatype}_tsvs/"

    printf -v PADDED_TABLE_ID "%03d" ~{table_id}

    # even for non-superpartitioned tables (e.g. metadata), the TSVs do have the suffix
    FILES="~{datatype}_${PADDED_TABLE_ID}_*"

    NUM_FILES=$(gsutil ls "${DIR}${FILES}" | wc -l)

    if [ ~{superpartitioned} = "true" ]; then
      TABLE="~{dataset_name}.${PREFIX}~{datatype}_${PADDED_TABLE_ID}"
    else
      TABLE="~{dataset_name}.${PREFIX}~{datatype}"
    fi

    if [ $NUM_FILES -gt 0 ]; then
        bq load --location=US --project_id=~{project_id} --skip_leading_rows=1 --source_format=CSV -F "\t" $TABLE $DIR$FILES ~{schema} || exit 1
        echo "ingested ${FILES} file from $DIR into table $TABLE"
        gsutil -m mv $DIR$FILES ${DIR}done/
    else
        echo "no ${FILES} files to process in $DIR"
    fi
  >>>

  runtime {
    docker: docker
    memory: "3 GB"
    disks: "local-disk 10 HDD"
    preemptible: 0
    cpu: 1
  }

  output {
    String done = "true"
  }
}

