version 1.0 

workflow GvsImportSample {

  input {
    File input_vcf
    String external_sample_name
    File interval_list
    String gvs_sample_id
    String project_id
    String dataset_name
    String? service_account_json_path
    String? drop_state = "SIXTY"
    Boolean? drop_state_includes_greater_than = false

    Int? preemptible_tries
    File? gatk_override = "gs://broad-dsp-spec-ops/scratch/andrea/gatk-ah-writeapi.jar"
    String? docker
  }

  String docker_final = select_first([docker, "us.gcr.io/broad-gatk/gatk:4.1.7.0"])
  
  call CheckForDuplicateData {
    input:
      project_id = project_id,
      dataset_name = dataset_name,
      sample_name = external_sample_name,
      gvs_sample_id = gvs_sample_id,
      service_account_json_path = service_account_json_path,
      preemptible_tries = preemptible_tries
  }

    call ImportSample {
      input:
        input_vcf = input_vcf,
        sample_name = external_sample_name,
        gvs_sample_id = gvs_sample_id,
        project_id = project_id,
        dataset_name = dataset_name,
        interval_list = interval_list,
        service_account_json_path = service_account_json_path,
        drop_state = drop_state,
        drop_state_includes_greater_than = drop_state_includes_greater_than,
        gatk_override = gatk_override,
        docker = docker_final,
        preemptible_tries = preemptible_tries,
        duplicate_check_passed = CheckForDuplicateData.done
    }
  

 
  output {
    Boolean loaded_in_gvs = true
  }
}


task CheckForDuplicateData {
    input{
      String project_id
      String dataset_name
      String sample_name
      String gvs_sample_id
      String? service_account_json_path
      # runtime
      Int? preemptible_tries
  }

  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'


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

    # check the INFORMATION_SCHEMA.PARTITIONS table to see if any of input sample names/ids have data loaded into their partitions
    # this returns the list of sample names that do already have data loaded
    bq --location=US --project_id=~{project_id} query --format=csv --use_legacy_sql=false \
      "SELECT distinct(partition_id) FROM ${INFO_SCHEMA_TABLE} p WHERE (p.partition_id = CAST(~{gvs_sample_id} AS STRING)) AND p.total_logical_bytes > 0 AND table_name like 'pet_%'" | \
      sed -e '/partition_id/d' > duplicates

    # true if there is data in results
    if [ -s duplicates ]; then
      echo "ERROR: Sample ~{sample_name} with gvs_id ~{gvs_sample_id} has already been loaded. Aborting"
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


task ImportSample {
  input {
    File input_vcf
    # File input_vcf_index
    String sample_name
    String project_id
    String dataset_name
    File interval_list
    String gvs_sample_id
    String? service_account_json_path
    String? drop_state
    Boolean? drop_state_includes_greater_than = false

    Boolean duplicate_check_passed

    # runtime
    Int? preemptible_tries
    File? gatk_override
    String docker
  }

  Int disk_size = if defined(drop_state) then 50 else 75
  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

  meta {
    description: "Imports sample data into BQ jointcalling dataset"
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

    if [ ~{has_service_account_file} = 'true' ]; then
        gsutil cp ~{service_account_json_path} local.service_account.json
        export GOOGLE_APPLICATION_CREDENTIALS=local.service_account.json
        gcloud auth activate-service-account --key-file=local.service_account.json
    fi

    # used if we localize the gvcf
    input_vcf_basename=$(basename ~{input_vcf})

    # default if we dont manually localize the gvcf
    updated_input_vcf=~{input_vcf}

    if [ ~{has_service_account_file} = 'true' ]; then
        gsutil cp ~{input_vcf} .
        gsutil cp ~{input_vcf}.tbi .
        updated_input_vcf=$input_vcf_basename
    fi

    # do we need this parameter?
    input_vcf_index="${updated_input_vcf}.tbi"

    gatk --java-options "-Xmx7000m" CreateVariantIngestFiles \
      -V ${updated_input_vcf} \
      -L ~{interval_list} \
      ~{"-IG " + drop_state} \
      --ignore-above-gq-threshold ~{drop_state_includes_greater_than} \
      --project-id ~{project_id} \
      --dataset-name ~{dataset_name} \
      --output-type BQ \
      --mode GENOMES \
      -SN ~{sample_name} \
      --gvs-sample-id ~{gvs_sample_id} \
      --ref-version 38
  >>>
  runtime {
      docker: docker
      memory: "3.75 GB"
      disks: "local-disk " + disk_size + " HDD"
      preemptible: select_first([preemptible_tries, 5])
      cpu: 1
  }
}







