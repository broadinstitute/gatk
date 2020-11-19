version 1.0

workflow ImportGenomes {

  input {
    Array[File] input_vcfs
    Array[File]? input_metrics
    File interval_list
    String output_directory
    File sample_map
    String project_id
    String dataset_name
    File pet_schema
    File vet_schema
    File metadata_schema
    String? drop_state

    Int? preemptible_tries
    File? gatk_override
    String? docker
  }

  String docker_final = select_first([docker, "us.gcr.io/broad-gatk/gatk:4.1.7.0"])

  call GetMaxTableId {
    input:
      sample_map = sample_map
  }

  scatter (i in range(length(input_vcfs))) {
    if (defined(input_metrics)) {
      File input_metric = select_first([input_metrics])[i]
    }

    call CreateImportTsvs {
      input:
        input_vcf = input_vcfs[i],
        interval_list = interval_list,
        input_metrics = input_metric,
        sample_map = sample_map,
        drop_state = drop_state,
        output_directory = output_directory,
        gatk_override = gatk_override,
        docker = docker_final,
        preemptible_tries = preemptible_tries
    }
  }

  scatter (i in range(GetMaxTableId.max_table_id)) {
    call LoadData as LoadMetadataTsvs{
      input:
        done = CreateImportTsvs.done,
        project_id = project_id,
        dataset_name = dataset_name,
        storage_location = output_directory,
        datatype = "metadata",
        numbered = "false",
        partitioned = "false",
        table_id = i + 1,
        schema = metadata_schema,
        preemptible_tries = preemptible_tries,
        docker = docker_final
    }

    call LoadData as LoadPetTsvs{
      input:
        done = CreateImportTsvs.done,
        project_id = project_id,
        dataset_name = dataset_name,
        storage_location = output_directory,
        datatype = "pet",
        table_id = i + 1,
        schema = pet_schema,
        preemptible_tries = preemptible_tries,
        docker = docker_final
    }

    call LoadData as LoadVetTsvs{
      input:
        done = CreateImportTsvs.done,
        project_id = project_id,
        dataset_name = dataset_name,
        storage_location = output_directory,
        datatype = "vet",
        table_id = i + 1,
        schema = vet_schema,
        preemptible_tries = preemptible_tries,
        docker = docker_final
    }
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
    File? input_metrics
    File interval_list
    String output_directory
    File sample_map
    String? drop_state

    # runtime
    Int? preemptible_tries
    File? gatk_override
    String docker

    String? for_testing_only
  }

  Int multiplier = if defined(drop_state) then 4 else 10
  Int disk_size = ceil(size(input_vcf, "GB") * multiplier) + 20

  meta {
    description: "Creates a tsv file for imort into BigQuery"
  }
  parameter_meta {
    input_vcf: {
      localization_optional: true
    }
  }
  command <<<
      set -e

      #workaround for https://github.com/broadinstitute/cromwell/issues/3647
      export TMPDIR=/tmp

      export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}
      ~{for_testing_only}

      gatk --java-options "-Xmx2500m" CreateVariantIngestFiles \
        -V ~{input_vcf} \
        -L ~{interval_list} \
        ~{"-IG " + drop_state} \
        ~{"-QCF " + input_metrics} \
        --mode GENOMES \
        -SNM ~{sample_map} \
        --ref-version 38
        
      gsutil cp metadata_*.tsv ~{output_directory}/metadata_tsvs/
      gsutil cp pet_*.tsv ~{output_directory}/pet_tsvs/
      gsutil cp vet_*.tsv ~{output_directory}/vet_tsvs/
  >>>
  runtime {
      docker: docker
      memory: "10 GB"
      disks: "local-disk " + disk_size + " HDD"
      preemptible: select_first([preemptible_tries, 5])
      cpu: 2
  }
  output {
      File metadata_tsv = glob("metadata_*.tsv")[0]
      File pet_tsv = glob("pet_*.tsv")[0] 
      File vet_tsv = glob("vet_*.tsv")[0]
      String done = "true"
  }
}


task LoadData {
  input {
    String project_id
    String dataset_name
    String storage_location
    String datatype
    Int table_id
    File schema
    String numbered = "true"
    String partitioned = "true"
    String load = "true"
    String uuid = ""

    #input from previous task needed to delay task from running until the other is complete
    Array[String] done

    # runtime
    Int? preemptible_tries
    String docker

    String? for_testing_only
  }

  command <<<
    set -x
    set -e
    ~{for_testing_only}

    DIR="~{storage_location}/~{datatype}_tsvs/"
    PARTITION_STRING=""

    if [ ~{partitioned} == "true" ]; then
      let "PARTITION_START=(~{table_id}-1)*4000+1"
      let "PARTITION_END=$PARTITION_START+3999"
      let "PARTITION_STEP=1"
      PARTITION_FIELD="sample_id"
      PARTITION_STRING="--range_partitioning=$PARTITION_FIELD,$PARTITION_START,$PARTITION_END,$PARTITION_STEP"
    fi

    # we are loading ONLY one table, specified by table_id
    printf -v PADDED_TABLE_ID "_%03d" ~{table_id}
    FILES="~{datatype}${PADDED_TABLE_ID}_*"
    
    if [ ~{numbered} != "true" ]; then
      PADDED_TABLE_ID=""  #override table id to empty string, but it is needed to get the files
    fi

    NUM_FILES=$(gsutil ls $DIR$FILES | wc -l)

    # create the table and load
    PREFIX=""
    if [ -n "~{uuid}" ]; then
        PREFIX="~{uuid}_"
    fi

    TABLE="~{dataset_name}.${PREFIX}~{datatype}${PADDED_TABLE_ID}"
    if [ $NUM_FILES -gt 0 ]; then
      set +e
      bq show --project_id ~{project_id} $TABLE > /dev/null
      BQ_SHOW_RC=$?
      set -e
      if [ $BQ_SHOW_RC -ne 0 ]; then

        # If the dataset does not exist, it will be caught here as bq mk will error out
        echo "making table $TABLE"
        bq --location=US mk ${PARTITION_STRING} --project_id=~{project_id} $TABLE ~{schema}
        #TODO: add a Google Storage Transfer for the table when we make it.
      fi
      #load should be false if using Google Storage Transfer so that the tables will be created by this script, but no data will be uploaded.
      if [ ~{load} = true ]; then
        bq load --location=US --project_id=~{project_id} --skip_leading_rows=1 --source_format=CSV -F "\t" $TABLE $DIR$FILES ~{schema}
        echo "ingested ${FILES} file from $DIR into table $TABLE"
        gsutil mv $DIR$FILES ${DIR}done/
      else
        echo "${FILES} will be ingested from $DIR by Google Storage Transfer"
      fi
    else
      echo "no ${FILES} files to process"
    fi
  >>>
  runtime {
    docker: docker
    memory: "3 GB"
    disks: "local-disk 10 HDD"
    preemptible: select_first([preemptible_tries, 5])
    cpu: 1
  }
}
