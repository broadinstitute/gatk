version 1.0

workflow ImportArrays {

  input {
    Array[File] input_vcfs
    Array[File]? input_metrics
    String? probe_info_table
    File? probe_info_file
    String output_directory
    File sample_map
    String project_id
    String dataset_name
    File raw_schema
    File sample_list_schema
    #TODO: determine table_id from input sample_map (including looping over multiple table_ids)
    Int table_id

    Int? preemptible_tries
    File? gatk_override
    String? docker
  }

  String docker_final = select_first([docker, "us.gcr.io/broad-gatk/gatk:4.1.7.0"])

  scatter (i in range(length(input_vcfs))) {
    if (defined(input_metrics)) {
      File input_metric = select_first([input_metrics])[i]
    }

    call CreateImportTsvs {
      input:
        input_vcf = input_vcfs[i],
        input_metrics = input_metric,
        probe_info_table = probe_info_table,
        probe_info_file = probe_info_file,
        sample_map = sample_map,
        output_directory = output_directory,
        gatk_override = gatk_override,
        docker = docker_final,
        preemptible_tries = preemptible_tries
    }
  }

  call LoadArrays {
    input:
      sample_tsvs = CreateImportTsvs.sample_tsv,
      project_id = project_id,
      dataset_name = dataset_name,
      storage_location = output_directory,
      table_id = table_id,
      raw_schema = raw_schema,
      sample_list_schema = sample_list_schema,
      preemptible_tries = preemptible_tries,
      docker = docker_final
  }
}


task CreateImportTsvs {
  input {
    File input_vcf
    File? input_metrics
    String? probe_info_table
    File? probe_info_file
    String output_directory
    File sample_map

    # runtime
    Int? preemptible_tries
    File? gatk_override
    String docker

    String? for_testing_only
  }

  Int disk_size = ceil(size(input_vcf, "GB") * 2.5) + 20

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

      gatk --java-options "-Xmx2500m" CreateArrayIngestFiles \
        -V ~{input_vcf} \
        ~{"-QCF " + input_metrics} \
        ~{"--probe-info-file " + probe_info_file} \
        ~{"--probe-info-table " + probe_info_table} \
        -SNM ~{sample_map} \
        --ref-version 37
        
      gsutil cp sample_*.tsv ~{output_directory}/sample_tsvs/
      gsutil cp raw_*.tsv ~{output_directory}/raw_tsvs/
  >>>
  runtime {
      docker: docker
      memory: "4 GB"
      disks: "local-disk " + disk_size + " HDD"
      preemptible: select_first([preemptible_tries, 5])
      cpu: 2
  }
  output {
      File sample_tsv = glob("sample_*.tsv")[0]
      File arraydata_tsv = glob("raw_*.tsv")[0] 
  }
}

task LoadArrays {
  input {
    String project_id
    String dataset_name
    String storage_location
    Int table_id
    File raw_schema
    File sample_list_schema
    String load = "true"
    String uuid = ""

    #input from previous task needed to delay task from running until the other is complete
    Array[String] sample_tsvs

    # runtime
    Int? preemptible_tries
    String docker

    String? for_testing_only
  }

  command <<<
    set -e
    ~{for_testing_only}

    SAMPLE_DIR=~{storage_location}/sample_tsvs/
    RAW_DIR=~{storage_location}/raw_tsvs/

    let "PARTITION_START=(~{table_id}-1)*4000+1"
    let "PARTITION_END=$PARTITION_START+3999"
    let "PARTITION_STEP=1"
    PARTITION_FIELD="sample_id"
    printf -v PADDED_TABLE_ID "%03d" ~{table_id}

    RAW_FILES="raw_${PADDED_TABLE_ID}_*"
    METADATA_FILES="sample_${PADDED_TABLE_ID}_*"

    NUM_RAW_FILES=$(gsutil ls $RAW_DIR${RAW_FILES} | wc -l)
    NUM_METADATA_FILES=$(gsutil ls $SAMPLE_DIR${METADATA_FILES} | wc -l)

    if [ $NUM_RAW_FILES -eq 0 -a $NUM_METADATA_FILES -eq 0 ]; then
      "no files for table ${PADDED_TABLE_ID} to process in ~{storage_location}; exiting"
      exit
    fi

    # create a metadata table and load
    SAMPLE_LIST_TABLE="~{dataset_name}.~{uuid + "_"}sample_list"
    if [ $NUM_METADATA_FILES -gt 0 ]; then
      set +e
      bq ls --project_id ~{project_id} ~{dataset_name} > /dev/null
      set -e
      if [ $? -ne 0 ]; then
        echo "making dataset ~{dataset_name}"
        bq mk --project_id=~{project_id} ~{dataset_name}
      fi
      set +e
      bq show --project_id ~{project_id} $SAMPLE_LIST_TABLE > /dev/null
      BQ_SHOW_RC=$?
      set -e
      if [ $BQ_SHOW_RC -ne 0 ]; then
        echo "making table $SAMPLE_LIST_TABLE"
        bq --location=US mk --project_id=~{project_id} $SAMPLE_LIST_TABLE ~{sample_list_schema}
        #TODO: add a Google Storage Transfer for the table when we make it.
      fi
      #load should be false if using Google Storage Transfer so that the tables will be created by this script, but no data will be uploaded.
      if [ ~{load} = true ]; then
        bq load --location=US --project_id=~{project_id} --skip_leading_rows=1 --null_marker="null" --source_format=CSV -F "\t" $SAMPLE_LIST_TABLE $SAMPLE_DIR$METADATA_FILES ~{sample_list_schema}
        echo "ingested ${METADATA_FILES} file from $SAMPLE_DIR into table $SAMPLE_LIST_TABLE"
      else
        echo "${METADATA_FILES} will be ingested from $SAMPLE_DIR by Google Storage Transfer"
      fi
    else
      echo "no metadata files to process"
    fi

    # create array table
    TABLE="~{dataset_name}.~{uuid + "_"}arrays_${PADDED_TABLE_ID}"
    if [ $NUM_RAW_FILES -gt 0 ]; then
      set +e
      bq show --project_id ~{project_id} $TABLE > /dev/null
      set -e
      if [ $? -ne 0 ]; then
        echo "making table $TABLE"
        bq --location=US mk --range_partitioning=$PARTITION_FIELD,$PARTITION_START,$PARTITION_END,$PARTITION_STEP \
          --project_id=~{project_id} $TABLE ~{raw_schema}
        #TODO: add a Google Storage Transfer for the table when we make it.
      fi
      if [ ~{load} = true ]; then
        bq load --location=US --project_id=~{project_id} --skip_leading_rows=1 --null_marker="null" --source_format=CSV -F "\t" $TABLE $RAW_DIR$RAW_FILES ~{raw_schema}
        echo "ingested ${RAW_FILES} files from $RAW_DIR into table $TABLE"
      else
        echo "${RAW_FILES} will be ingested from $RAW_DIR
         by Google Storage Transfer"
      fi
    else
      echo "no raw data files to process"
    fi
  >>>
  runtime {
    docker: docker
    memory: "4 GB"
    disks: "local-disk 10 HDD"
    preemptible: select_first([preemptible_tries, 5])
    cpu: 2
  }
}
