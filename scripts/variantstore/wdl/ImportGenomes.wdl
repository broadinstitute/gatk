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
    Boolean? drop_state_includes_greater_than = false

    Int? preemptible_tries
    File? gatk_override
    String? docker
  }

  String docker_final = select_first([docker, "us.gcr.io/broad-gatk/gatk:4.1.7.0"])

  call GetMaxTableId {
    input:
      sample_map = sample_map
  }

  # create tables only requires GetMaxTableId
  call CreateTables as CreateMetadataTables {
  	input:
      project_id = project_id,
      dataset_name = dataset_name,
      storage_location = output_directory,
      datatype = "metadata",
      max_table_id = GetMaxTableId.max_table_id,
      schema = metadata_schema,
      numbered = "false",
      partitioned = "false",
      uuid = "",
      preemptible_tries = preemptible_tries,
      docker = docker_final
  }

  call CreateTables as CreatePetTables {
  	input:
      project_id = project_id,
      dataset_name = dataset_name,
      storage_location = output_directory,
      datatype = "pet",
      max_table_id = GetMaxTableId.max_table_id,
      schema = pet_schema,
      numbered = "false",
      partitioned = "false",
      uuid = "",
      preemptible_tries = preemptible_tries,
      docker = docker_final
  }

  call CreateTables as CreateVetTables {
  	input:
      project_id = project_id,
      dataset_name = dataset_name,
      storage_location = output_directory,
      datatype = "vet",
      max_table_id = GetMaxTableId.max_table_id,
      schema = vet_schema,
      numbered = "false",
      partitioned = "false",
      uuid = "",
      preemptible_tries = preemptible_tries,
      docker = docker_final
  }

  # create the VCFs
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
        drop_state_includes_greater_than = drop_state_includes_greater_than,
        output_directory = output_directory,
        gatk_override = gatk_override,
        docker = docker_final,
        preemptible_tries = preemptible_tries
    }
  }

  # loading tables requires create tables and create import tsvs to be done
  scatter (table_dir_files_str in CreateMetadataTables.table_dir_files_list) {
    call LoadTable as LoadMetadataTable {
      input:
        done = CreateImportTsvs.done,
        table_dir_files_str = table_dir_files_str,
        project_id = project_id,
        schema = metadata_schema,
        load = "true",
        preemptible_tries = preemptible_tries,
        docker = docker
    }
  }

  scatter (table_dir_files_str in CreatePetTables.table_dir_files_list) {
    call LoadTable as LoadPetTable {
      input:
        done = CreateImportTsvs.done,
        table_dir_files_str = table_dir_files_str,
        project_id = project_id,
        schema = pet_schema,
        load = "true",
        preemptible_tries = preemptible_tries,
        docker = docker
    }
  }

  scatter (table_dir_files_str in CreateVetTables.table_dir_files_list) {
    call LoadTable as LoadVetTable {
      input:
        done = CreateImportTsvs.done,
        table_dir_files_str = table_dir_files_str,
        project_id = project_id,
        schema = vet_schema,
        load = "true",
        preemptible_tries = preemptible_tries,
        docker = docker
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
    Boolean? drop_state_includes_greater_than = false

    # runtime
    Int? preemptible_tries
    File? gatk_override
    String docker
  }

  Int multiplier = if defined(drop_state) then 4 else 10
  Int disk_size = ceil(size(input_vcf, "GB") * multiplier) + 20

  meta {
    description: "Creates a tsv file for import into BigQuery"
    volatile: true
  }

  parameter_meta {
    input_vcf: {
      localization_optional: true
    }
  }
  command <<<
      set -e

      # workaround for https://github.com/broadinstitute/cromwell/issues/3647
      export TMPDIR=/tmp

      export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

      gatk --java-options "-Xmx7000m" CreateVariantIngestFiles \
        -V ~{input_vcf} \
        -L ~{interval_list} \
        ~{"-IG " + drop_state} \
        --ignore-above-gq-threshold ~{drop_state_includes_greater_than} \
        ~{"-QCF " + input_metrics} \
        --mode GENOMES \
        -SNM ~{sample_map} \
        --ref-version 38

      gsutil -m cp metadata_*.tsv ~{output_directory}/metadata_tsvs/
      gsutil -m cp pet_*.tsv ~{output_directory}/pet_tsvs/
      gsutil -m cp vet_*.tsv ~{output_directory}/vet_tsvs/
  >>>
  runtime {
      docker: docker
      memory: "10 GB"
      disks: "local-disk " + disk_size + " HDD"
      preemptible: select_first([preemptible_tries, 5])
      cpu: 2
  }
  output {
      String done = "true"
  }
}

# Creates all the tables necessary for the LoadData operation
# As an optimization, I also generate a (table, dir, files) csv file which contains
# most of inputs necessary for the following LoadTable task.
task CreateTables {
	meta {
    	volatile: true
  	}

	input {
      String project_id
      String dataset_name
      String storage_location
      String datatype
      Int max_table_id
      File schema
      String numbered
      String partitioned
      String uuid

      # runtime
      Int? preemptible_tries
      String docker
    }

  command <<<
    set -x
    set -e

    DIR="~{storage_location}/~{datatype}_tsvs/"

    for TABLE_ID in $(seq 1 ~{max_table_id}); do
      PARTITION_STRING=""
      if [ ~{partitioned} == "true" ]; then
        let "PARTITION_START=(${TABLE_ID}-1)*4000+1"
        let "PARTITION_END=$PARTITION_START+3999"
        let "PARTITION_STEP=1"
        PARTITION_FIELD="sample_id"
        PARTITION_STRING="--range_partitioning=$PARTITION_FIELD,$PARTITION_START,$PARTITION_END,$PARTITION_STEP"
      fi

      printf -v PADDED_TABLE_ID "_%03d" ${TABLE_ID}
      FILES="~{datatype}${PADDED_TABLE_ID}_*"

      NUM_FILES=$(gsutil ls $DIR$FILES | wc -l)

      # create the table
      PREFIX=""
      if [ -n "~{uuid}" ]; then
          PREFIX="~{uuid}_"
      fi

      if [ $NUM_FILES -gt 0 ]; then
        if [ ~{numbered} != "true" ]; then
          PADDED_TABLE_ID=""  #override table id to empty string, but it is needed to get the files
        fi

        TABLE="~{dataset_name}.${PREFIX}~{datatype}${PADDED_TABLE_ID}"

        # Check that the table has not been created yet
        set +e
        bq show --project_id ~{project_id} $TABLE > /dev/null
        BQ_SHOW_RC=$?
        set -e
        if [ $BQ_SHOW_RC -ne 0 ]; then
          echo "making table $TABLE"
          bq --location=US mk ${PARTITION_STRING} --project_id=~{project_id} $TABLE ~{schema}
        fi

        echo "$TABLE,$DIR,$FILES" >> table_dir_files.csv
      else
        echo "no ${FILES} files to process"
      fi

    #done
  >>>

  output {
    Array[String] table_dir_files_list = read_lines("table_dir_files.csv")
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
    String table_dir_files_str
    String project_id
    File schema
    String load
    String done

    Int? preemptible_tries
    String docker
  }

  command <<<
    TABLE=$(echo ~{table_dir_files_str} | cut -d, -f1)
    DIR=$(echo ~{table_dir_files_str} | cut -d, -f2)
    FILES=$(echo ~{table_dir_files_str} | cut -d, -f3)

    #load should be false if using Google Storage Transfer so that the tables will be created by this script, but no data will be uploaded.
    if [ ~{load} = true ]; then
      bq load --location=US --project_id=~{project_id} --skip_leading_rows=1 --source_format=CSV -F "\t" $TABLE $DIR$FILES ~{schema} || exit 1
      echo "ingested ${FILES} file from $DIR into table $TABLE"
      gsutil mv $DIR$FILES ${DIR}done/
    else
      echo "${FILES} will be ingested from $DIR by Google Storage Transfer"
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
