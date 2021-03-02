version 1.0

workflow ImportGenomes {

  input {
    Array[File] input_vcfs
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

  call CreateTables as CreateMetadataTables {
  	input:
      project_id = project_id,
      dataset_name = dataset_name,
      storage_location = output_directory,
      datatype = "metadata",
      max_table_id = GetMaxTableId.max_table_id,
      schema = metadata_schema,
      superpartitioned = "false",
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
      superpartitioned = "true",
      partitioned = "true",
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
      superpartitioned = "true",
      partitioned = "true",
      uuid = "",
      preemptible_tries = preemptible_tries,
      docker = docker_final
  }

  scatter (i in range(length(input_vcfs))) {
    call CreateImportTsvs {
      input:
        input_vcf = input_vcfs[i],
        interval_list = interval_list,
        sample_map = sample_map,
        drop_state = drop_state,
        drop_state_includes_greater_than = drop_state_includes_greater_than,
        output_directory = output_directory,
        gatk_override = gatk_override,
        docker = docker_final,
        preemptible_tries = preemptible_tries
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
        schema = metadata_schema,
        table_creation_done = CreateMetadataTables.done,
        tsv_creation_done = CreateImportTsvs.done,
        docker = docker_final
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
      table_creation_done = CreatePetTables.done,
      tsv_creation_done = CreateImportTsvs.done,
      docker = docker_final
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
      table_creation_done = CreateVetTables.done,
      tsv_creation_done = CreateImportTsvs.done,
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
task CreateTables {
	meta {
    	volatile: true
  	}

	input {
      String project_id
      String dataset_name
      String datatype
      Int max_table_id
      File schema
      String superpartitioned
      String partitioned
      String uuid

      # runtime
      Int? preemptible_tries
      String docker
    }

  command <<<
    set -x
    set -e

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
    File schema
    String table_creation_done
    Array[String] tsv_creation_done

    String docker
  }

  command <<<
    set -x
    set -e

    DIR="~{storage_location}/~{datatype}_tsvs/"

    printf -v PADDED_TABLE_ID "%03d" ~{table_id}

    # even for non-superpartitioned tables (e.g. metadata), the TSVs do have the suffix
    FILES="~{datatype}_${PADDED_TABLE_ID}_*"

    if [ ~{superpartitioned} = "true" ]; then
      TABLE="~{dataset_name}.${PREFIX}~{datatype}_${PADDED_TABLE_ID}"
    else
      TABLE="~{dataset_name}.${PREFIX}~{datatype}"
    fi

    bq load --location=US --project_id=~{project_id} --skip_leading_rows=1 --source_format=CSV -F "\t" $TABLE $DIR$FILES ~{schema} || exit 1
    echo "ingested ${FILES} file from $DIR into table $TABLE"
    gsutil mv $DIR$FILES ${DIR}done/
  >>>

  runtime {
    docker: docker
    memory: "3 GB"
    disks: "local-disk 10 HDD"
    preemptible: 0
    cpu: 1
  }
}
