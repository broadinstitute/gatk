version 1.0

workflow GvsImportGenomes {

  input {
    Array[File] input_vcfs
    Array[File] input_vcf_indexes
    Array[String] external_sample_names
    File interval_list
    String output_directory
    File sample_map
    String project_id
    String dataset_name
    String? pet_schema = "location:INTEGER,sample_id:INTEGER,state:STRING"
    String? vet_schema = "sample_id:INTEGER,location:INTEGER,ref:STRING,alt:STRING,AS_RAW_MQ:STRING,AS_RAW_MQRankSum:STRING,QUALapprox:STRING,AS_QUALapprox:STRING,AS_RAW_ReadPosRankSum:STRING,AS_SB_TABLE:STRING,AS_VarDP:STRING,call_GT:STRING,call_AD:STRING,call_GQ:INTEGER,call_PGT:STRING,call_PID:STRING,call_PL:STRING"
    String? sample_info_schema = "sample_name:STRING,sample_id:INTEGER,inferred_state:STRING"
    File? service_account_json
    String? drop_state
    Boolean? drop_state_includes_greater_than = false

    Int? preemptible_tries
    File? gatk_override
    String? docker
  }

  String docker_final = select_first([docker, "us.gcr.io/broad-gatk/gatk:4.1.7.0"])

  # how do i return an error if the lengths are not equal?
  Int input_length = length(input_vcfs)
  if (input_length != length(input_vcf_indexes) || input_length != length(external_sample_names)) {
    call TerminateWorkflow {
      input:
        message = 'Input array lengths do not match'
    }
  }

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

  call CreateTables as CreateSampleInfoTables {
  	input:
      project_id = project_id,
      dataset_name = dataset_name,
      datatype = "sample_info",
      max_table_id = GetMaxTableId.max_table_id,
      schema = sample_info_schema,
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
      schema = pet_schema,
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
      schema = vet_schema,
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
        sample_name = external_sample_names[i],
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
    call LoadTable as LoadSampleInfoTable {
      input:
        project_id = project_id,
        table_id = i + 1,
        dataset_name = dataset_name,
        storage_location = output_directory,
        datatype = "sample_info",
        superpartitioned = "false",
        schema = sample_info_schema,
        service_account_json = service_account_json,
        table_creation_done = CreateSampleInfoTables.done,
        tsv_creation_done = CreateImportTsvs.done,
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
      schema = pet_schema,
      service_account_json = service_account_json,
      table_creation_done = CreatePetTables.done,
      tsv_creation_done = CreateImportTsvs.done,
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
      schema = vet_schema,
      service_account_json = service_account_json,
      table_creation_done = CreateVetTables.done,
      tsv_creation_done = CreateImportTsvs.done,
      docker = docker_final,
      run_uuid = SetLock.run_uuid
    }
  }

  call ReleaseLock {
    input:
      run_uuid = SetLock.run_uuid,
      output_directory = output_directory,
      load_sample_info_done = LoadSampleInfoTable.done,
      load_pet_done = LoadPetTable.done,
      load_vet_done = LoadVetTable.done,
      service_account_json = service_account_json,
      preemptible_tries = preemptible_tries
  }

  output {
    Boolean loaded_in_gvs = true
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
      echo "ERROR: lock file in place. Check whether another run of ImportGenomes with this output directory is in progress or a previous run had an error.
            If you would like to proceed, run the following command and re-run the workflow: \
            gsutil rm ${DIR}${LOCKFILE} \
            " && exit 1
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
    Array[String] load_sample_info_done
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
      export GOOGLE_APPLICATION_CREDENTIALS=~{service_account_json}
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
    String sample_name
    File interval_list
    String output_directory
    File sample_map
    File? service_account_json
    String? drop_state
    Boolean? drop_state_includes_greater_than = false

    Boolean call_cache_tsvs = false

    # runtime
    Int? preemptible_tries
    File? gatk_override
    String docker

    String? for_testing_only
    String run_uuid
  }

  Int disk_size = if defined(drop_state) then 30 else 75
  String input_vcf_basename = basename(input_vcf)
  String has_service_account_file = if (defined(service_account_json)) then 'true' else 'false'
  # if we are doing a manual localization, we need to set the filename
  String updated_input_vcf = if (defined(service_account_json)) then input_vcf_basename else input_vcf

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
        export GOOGLE_APPLICATION_CREDENTIALS=~{service_account_json}
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

      # check whether these files have already been generated
      DO_TSV_GENERATION='true'
      if [ ~{call_cache_tsvs} = 'true' ]; then
        echo "Checking for files to call cache"

        declare -a TABLETYPES=("sample_info" "pet" "vet")
        ALL_FILES_EXIST='true'
        for TABLETYPE in ${TABLETYPES[@]}; do
            FILEPATH="~{output_directory}/${TABLETYPE}_tsvs/**${TABLETYPE}_*_~{input_vcf_basename}.tsv"
            # output 1 if no file is found
            result=$(gsutil ls $FILEPATH || echo 1)

            if [ $result == 1 ]; then
              echo "A file matching ${FILEPATH} does not exist"
              ALL_FILES_EXIST='false'
            else
              if [[ $result = *"/done/"* ]]; then
                echo "File ${FILENAME} seems to have been processed already; found at ${result}"
                echo "Something is very wrong!"
                exit 1
              elif [[ $result = "~{output_directory}/${TABLETYPE}_tsvs/set_"* ]]; then
                FILENAME=$(basename $result)
                echo "File ${FILENAME} is in a set directory. Moving out of set directory to ~{output_directory}/${TABLETYPE}_tsvs/${FILENAME}"
                gsutil mv $result "~{output_directory}/${TABLETYPE}_tsvs/"
              fi

            fi
        done

        if [ $ALL_FILES_EXIST = 'true' ]; then
            DO_TSV_GENERATION='false'
            echo "Skipping TSV generation for input file ~{input_vcf_basename} because the output TSV files already exist."
        fi
      fi

      if [ $DO_TSV_GENERATION = 'true' ]; then
          echo "Generating TSVs for input file ~{input_vcf_basename}"

          # TODO in future when we pass the source path or gvs_id as an arg here, use a version of that for the filename too
          gatk --java-options "-Xmx7000m" CreateVariantIngestFiles \
            -V ~{updated_input_vcf} \
            -L ~{interval_list} \
            ~{"-IG " + drop_state} \
            --ignore-above-gq-threshold ~{drop_state_includes_greater_than} \
            --mode GENOMES \
            -SN ~{sample_name} \
            -SNM ~{sample_map} \
            --ref-version 38

          gsutil -m cp sample_info_*.tsv ~{output_directory}/sample_info_tsvs/
          gsutil -m cp pet_*.tsv ~{output_directory}/pet_tsvs/
          gsutil -m cp vet_*.tsv ~{output_directory}/vet_tsvs/
      fi
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
      String? schema
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
      export GOOGLE_APPLICATION_CREDENTIALS=~{service_account_json}
      gcloud auth activate-service-account --key-file='~{service_account_json}'
      gcloud config set project ~{project_id}
    fi

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
        let "PARTITION_END=$PARTITION_START+3999"
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
    String? schema
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

    # even for non-superpartitioned tables (e.g. sample_info), the TSVs do have the suffix
    FILES="~{datatype}_${PADDED_TABLE_ID}_*"

    NUM_FILES=$(gsutil ls "${DIR}${FILES}" | wc -l)

    if [ ~{superpartitioned} = "true" ]; then
      TABLE="~{dataset_name}.${PREFIX}~{datatype}_${PADDED_TABLE_ID}"
    else
      TABLE="~{dataset_name}.${PREFIX}~{datatype}"
    fi

    if [ $NUM_FILES -gt 0 ]; then
        # get list of of pet files and their byte sizes
        echo "Getting file sizes(bytes), paths to each file, and determining sets for chunking."
        echo -e "bytes\tfile_path\tsum_bytes\tset_number" > ~{datatype}_du_sets.txt
        # tr to replace each space -> tab, squeeze (remove) "empty" tabs,
        gsutil du "${DIR}${FILES}" | tr " " "\t" | tr -s "\t" | sed "/~{datatype}_tsvs\/$/d" | awk '{s+=$1}{print $1"\t"$2"\t"s"\t" (1+int(s / 16000000000000))}' >> ~{datatype}_du_sets.txt

        # per set, load table
        for set in $(sed 1d ~{datatype}_du_sets.txt | cut -f4 | sort | uniq)
        do
          # move set data into separate directory
          echo "Moving set $set data into separate directory."
          awk -v set="$set" '$4 == set' ~{datatype}_du_sets.txt | cut -f2 | gsutil -m mv -I "${DIR}set_${set}/" 2> copy.log

          # execute bq load command, get bq load job id, add details per set to tmp file
          echo "Running BigQuery load for set $set."
          bq load --quiet --nosync --location=US --project_id=~{project_id} --skip_leading_rows=1 --source_format=CSV -F "\t" \
            "$TABLE" "${DIR}set_${set}/${FILES}" ~{schema} > status_bq_submission

          cat status_bq_submission | tail -n 1 > status_bq_submission_last_line

          bq_job_id=$(sed 's/.*://' status_bq_submission_last_line)

          echo $bq_job_id
          # add job ID as key and gs path to the data set uploaded as value
          echo -e "${bq_job_id}\t${set}\t${DIR}set_${set}/" >> bq_load_details.tmp
        done

        # for each bq load job submitted, run bq wait, capture success/failure to tmp file
        while IFS="\t" read -r line_bq_load
        do
          bq wait --project_id=~{project_id} $(echo "$line_bq_load" | cut -f1) > bq_wait_status

          # capture SUCCESS or FAILURE, echo to file
          wait_status=$(sed '6q;d' bq_wait_status | tr " " "\t" | tr -s "\t" | cut -f3)
          echo "$wait_status" >> bq_wait_details.tmp
        done < bq_load_details.tmp

        # combine load status and wait status into final report
        paste bq_load_details.tmp bq_wait_details.tmp > bq_final_job_statuses.txt

        # move files from each set into set-level "done" directories
        gsutil -m mv "${DIR}set_${set}/${FILES}" "${DIR}set_${set}/done/" 2> gsutil_mv_done.log

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
    File? manifest_file = "~{datatype}_du_sets.txt"
    File? final_job_statuses = "bq_final_job_statuses.txt"
  }
}

task TerminateWorkflow {
  input {
    String message
  }

  command <<<
      set -e
      echo ~{message}
      exit 1
  >>>
  runtime {
      docker: "python:3.8-slim-buster"
      memory: "1 GB"
      disks: "local-disk 10 HDD"
      preemptible: 3
      cpu: 1
  }
}
