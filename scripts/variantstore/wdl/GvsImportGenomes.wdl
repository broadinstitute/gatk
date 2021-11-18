version 1.0

workflow GvsImportGenomes {

  input {
    Array[File] input_vcfs
    Array[File] input_vcf_indexes
    Array[String] external_sample_names
    File interval_list
    String project_id
    String dataset_name

    String? service_account_json_path
    String? drop_state = "FORTY"
    Boolean? drop_state_includes_greater_than = false

    Int batch_size = 1

    Int? preemptible_tries
    File? gatk_override = "gs://broad-dsp-spec-ops/scratch/bigquery-jointcalling/jars/ah_var_store_20210914/gatk-package-4.2.0.0-406-ga9206a2-SNAPSHOT-local.jar"
    String? docker
  }

  String docker_final = select_first([docker, "us.gcr.io/broad-gatk/gatk:4.1.7.0"])

  # return an error if the lengths are not equal
  Int input_length = length(input_vcfs)
  if (input_length != length(external_sample_names)) {
    call TerminateWorkflow {
      input:
        message = 'Input array lengths do not match'
    }
  }

  call GetSampleIds {
    input:
      external_sample_names = external_sample_names,
      project_id = project_id,
      dataset_name = dataset_name,
      table_name = "sample_info",
      service_account_json_path = service_account_json_path
  }

  call CheckForDuplicateData {
    input:
      project_id = project_id,
      dataset_name = dataset_name,
      sample_names = external_sample_names,
      service_account_json_path = service_account_json_path
  }

  call CreateFOFNs {
    input:
        input_vcf_list = write_lines(input_vcfs),
        input_vcf_index_list = write_lines(input_vcf_indexes),
        sample_name_list = write_lines(external_sample_names),
        batch_size = batch_size,
  }

  scatter (i in range(length(CreateFOFNs.vcf_batch_vcf_fofns))) {
    call LoadData {
      input:
        project_id = project_id,
        dataset_name = dataset_name,
        input_vcfs = read_lines(CreateFOFNs.vcf_batch_vcf_fofns[i]),
        input_vcf_indexes = read_lines(CreateFOFNs.vcf_batch_vcf_index_fofns[i]),
        sample_names = read_lines(CreateFOFNs.vcf_sample_name_fofns[i]),
        interval_list = interval_list,
        service_account_json_path = service_account_json_path,
        sample_map = GetSampleIds.sample_map,
        drop_state = drop_state,
        drop_state_includes_greater_than = drop_state_includes_greater_than,
        gatk_override = gatk_override,
        docker = docker_final,
        preemptible_tries = preemptible_tries,
        duplicate_check_passed = CheckForDuplicateData.done
    }
  }

  call SetIsLoadedColumn {
    input:
      load_done = LoadData.done,
      dataset_name = dataset_name,
      gvs_ids = GetSampleIds.gvs_ids,
      service_account_json_path = service_account_json_path,
      project_id = project_id,
      preemptible_tries = preemptible_tries
  }


  output {
    Boolean loaded_in_gvs = true
  }
}


task CheckForDuplicateData {
    input {
      String project_id
      String dataset_name
      Array[String] sample_names
      String? service_account_json_path
      # runtime
      Int? preemptible_tries
    }

  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'
  Int num_samples = length(sample_names)

  meta {
    volatile: true
  }

  command <<<
    set -e

    if [ ~{has_service_account_file} = 'true' ]; then
      gsutil cp ~{service_account_json_path} local.service_account.json
      gcloud auth activate-service-account --key-file=local.service_account.json
      gcloud config set project ~{project_id}
    fi

    echo "project_id = ~{project_id}" > ~/.bigqueryrc

    INFO_SCHEMA_TABLE="~{dataset_name}.INFORMATION_SCHEMA.PARTITIONS"
    TEMP_TABLE="~{dataset_name}.sample_dupe_check"
    SAMPLE_INFO_TABLE="~{dataset_name}.sample_info"

    # create a temp table with the sample_names
    bq --project_id=~{project_id} mk ${TEMP_TABLE} "sample_name:STRING"
    NAMES_FILE=~{write_lines(sample_names)}
    bq load --project_id=~{project_id} ${TEMP_TABLE} $NAMES_FILE "sample_name:STRING"

    # check the INFORMATION_SCHEMA.PARTITIONS table to see if any of input sample names/ids have data loaded into their partitions
    # this returns the list of sample names that do already have data loaded
    bq --location=US --project_id=~{project_id} query --format=csv -n ~{num_samples} --use_legacy_sql=false \
      "WITH items as (SELECT s.sample_id, s.sample_name FROM ${TEMP_TABLE} t left outer join ${SAMPLE_INFO_TABLE} s on (s.sample_name = t.sample_name)) " \
      "SELECT i.sample_name FROM ${INFO_SCHEMA_TABLE} p JOIN items i ON (p.partition_id = CAST(i.sample_id AS STRING)) WHERE p.total_logical_bytes > 0 AND table_name like 'pet_%'" | \
      sed -e '/sample_name/d' > duplicates

      # remove the temp table
      bq --project_id=~{project_id} rm -f -t ${TEMP_TABLE}

    # true if there is data in results
    if [ -s duplicates ]; then
      echo "ERROR: Trying to load samples that have already been loaded"
      cat duplicates
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
  output {
      Boolean done = true
      File? duplicates = "duplicates"
  }
}

task CreateFOFNs {
    input {
        File input_vcf_list
        File input_vcf_index_list
        File sample_name_list
        Int batch_size
    }

     command {
         set -e

         split -d -a 5 -l ~{batch_size} ~{input_vcf_list} batched_vcfs.
         split -d -a 5 -l ~{batch_size} ~{input_vcf_index_list} batched_vcf_indexes.
         split -d -a 5 -l ~{batch_size} ~{sample_name_list} batched_sample_names.
     }

     runtime {
         docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
         bootDiskSizeGb: 15
         memory: "3 GB"
         disks: "local-disk 10 HDD"
         preemptible: 3
         cpu: 1
     }

     output {
         Array[File] vcf_batch_vcf_fofns = glob("batched_vcfs.*")
         Array[File] vcf_batch_vcf_index_fofns = glob("batched_vcf_indexes.*")
         Array[File] vcf_sample_name_fofns = glob("batched_sample_names.*")
     }
}

task LoadData {
  input {
    Array[File] input_vcfs
    Array[File] input_vcf_indexes
    Array[String] sample_names
    File interval_list
    File sample_map
    String? service_account_json_path
    String? drop_state
    Boolean? drop_state_includes_greater_than = false

    Boolean duplicate_check_passed

    String project_id
    String dataset_name

    Boolean load_ref_ranges = true
    Boolean load_pet = false
    Boolean load_vet = true

    # runtime
    Int? preemptible_tries
    File? gatk_override
    String docker

  }

  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

  meta {
    description: "Load data into BigQuery using the Write Api"

    # TODO: is this right?  should this call cache?
    volatile: true
  }

  parameter_meta {
    input_vcfs: {
      localization_optional: true
    }
  }

  command <<<
      set -e

      # workaround for https://github.com/broadinstitute/cromwell/issues/3647
      export TMPDIR=/tmp

      export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

      if [ ~{has_service_account_file} = 'true' ]; then
          gsutil cp ~{service_account_json_path} local.service_account.json
          export GOOGLE_APPLICATION_CREDENTIALS=local.service_account.json
          gcloud auth activate-service-account --key-file=local.service_account.json
      fi

      # translate WDL arrays into BASH arrays
      VCFS_ARRAY=(~{sep=" " input_vcfs})
      VCF_INDEXES_ARRAY=(~{sep=" " input_vcf_indexes})
      SAMPLE_NAMES_ARRAY=(~{sep=" " sample_names})

      # loop over the BASH arrays (See https://stackoverflow.com/questions/6723426/looping-over-arrays-printing-both-index-and-value)
      for i in "${!VCFS_ARRAY[@]}"; do
          input_vcf="${VCFS_ARRAY[$i]}"
          input_vcf_basename=$(basename $input_vcf)
          updated_input_vcf=$input_vcf
          input_vcf_index="${VCF_INDEXES_ARRAY[$i]}"
          sample_name="${SAMPLE_NAMES_ARRAY[$i]}"

          # we always do our own localization
          gsutil cp $input_vcf .
          gsutil cp $input_vcf_index .
          updated_input_vcf=$input_vcf_basename

          gatk --java-options "-Xmx2g" CreateVariantIngestFiles \
              -V ${updated_input_vcf} \
              -L ~{interval_list} \
              ~{"-IG " + drop_state} \
              --ignore-above-gq-threshold ~{drop_state_includes_greater_than} \
              --project-id ~{project_id} \
              --dataset-name ~{dataset_name} \
              --output-type BQ \
              --enable-reference-ranges ~{load_ref_ranges} \
              --enable-pet ~{load_pet} \
              --enable-vet ~{load_vet} \
              -SN ${sample_name} \
              -SNM ~{sample_map} \
              --ref-version 38

      done

  >>>
  runtime {
      docker: docker
      memory: "3.75 GB"
      disks: "local-disk 50 HDD"
      preemptible: select_first([preemptible_tries, 5])
      cpu: 1
  }
  output {
      String done = "true"
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

task SetIsLoadedColumn {
  meta {
    volatile: true
  }

  input {
    Array[String] load_done
    String dataset_name
    String project_id
    File gvs_ids
    String? service_account_json_path
    Int? preemptible_tries
  }

  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'
  Array[String] gvs_id_array = read_lines(gvs_ids)

  command <<<
    set -ex

    if [ ~{has_service_account_file} = 'true' ]; then
      gsutil cp ~{service_account_json_path} local.service_account.json
      gcloud auth activate-service-account --key-file=local.service_account.json
      gcloud config set project ~{project_id}
    fi

    echo "project_id = ~{project_id}" > ~/.bigqueryrc

    # set is_loaded to true if there is a corresponding vet table partition with rows for that sample_id
    bq --location=US --project_id=~{project_id} query --format=csv --use_legacy_sql=false \
    "UPDATE ~{dataset_name}.sample_info SET is_loaded = true WHERE sample_id IN (SELECT CAST(partition_id AS INT64) from ~{dataset_name}.INFORMATION_SCHEMA.PARTITIONS WHERE partition_id in ('~{sep="\',\'" gvs_id_array}') AND total_logical_bytes > 0 AND table_name LIKE \"vet_%\")"
  >>>

  runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
    memory: "1 GB"
    disks: "local-disk 10 HDD"
    preemptible: select_first([preemptible_tries, 5])
    cpu: 1
  }

  output {
    String done = "done"
  }
}

task GetSampleIds {
  meta {
    volatile: true
  }

  input {
    Array[String] external_sample_names
    String project_id
    String dataset_name
    String table_name
    String? service_account_json_path
    Int samples_per_table = 4000

    # runtime
    Int? preemptible_tries
  }

  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'
  Int num_samples = length(external_sample_names)

  command <<<
      set -ex
      if [ ~{has_service_account_file} = 'true' ]; then
        gsutil cp ~{service_account_json_path} local.service_account.json
        gcloud auth activate-service-account --key-file=local.service_account.json
      fi

      echo "project_id = ~{project_id}" > ~/.bigqueryrc

      # get the current maximum id, or 0 if there are none
      bq --project_id=~{project_id} query --format=csv --use_legacy_sql=false \
        "SELECT IFNULL(MIN(sample_id),0) as min, IFNULL(MAX(sample_id),0) as max FROM ~{dataset_name}.~{table_name} where sample_name in ('~{sep="\',\'" external_sample_names}')" > results

      # prep for being able to return min table id
      min_sample_id=$(tail -1 results | cut -d, -f1)
      max_sample_id=$(tail -1 results | cut -d, -f2)

      python3 -c "from math import ceil; print(ceil($max_sample_id/~{samples_per_table}))" > max_sample_id
      python3 -c "from math import ceil; print(ceil($min_sample_id/~{samples_per_table}))" > min_sample_id

      bq --project_id=~{project_id} query --format=csv --use_legacy_sql=false -n ~{num_samples} \
        "SELECT sample_id, sample_name FROM ~{dataset_name}.~{table_name} where sample_name in ('~{sep="\',\'" external_sample_names}')" > sample_map

      cut -d, -f1 sample_map > gvs_ids

  >>>
  runtime {
      docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
      memory: "1 GB"
      disks: "local-disk 10 HDD"
      preemptible: select_first([preemptible_tries, 5])
      cpu: 1
  }
  output {
      Int max_table_id = ceil(read_float("max_sample_id"))
      Int min_table_id = ceil(read_float("min_sample_id"))
      File sample_map = "sample_map"
      File gvs_ids = "gvs_ids"
  }
}
