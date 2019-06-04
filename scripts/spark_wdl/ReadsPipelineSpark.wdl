workflow ReadsPipelineSparkWorkflow {
  String ref_fasta
  String input_bam
  String output_vcf
  String known_sites
  String emit_reference_confidence

  String gatk_docker
  String gatk
  String gatk_spark_jar
  File? gcloud_service_account_key_file
  String gcloud_project
  String? gcloud_zone
  String? cluster_name
  String? gatk_gcs_staging
  Int? num_workers
  String? image_version
  String? max_age
  Int executor_cores
  String executor_memory
  String driver_memory
  Boolean? delete_cluster

  call CreateDataprocCluster {
    input:
      gatk_docker = gatk_docker,
      gcloud_service_account_key_file = gcloud_service_account_key_file,
      gcloud_project = gcloud_project,
      gcloud_zone = gcloud_zone,
      cluster_name = cluster_name,
      num_workers = num_workers,
      image_version = image_version,
      max_age = max_age
  }

  call ReadsPipelineSpark {
    input:
      ref_fasta = ref_fasta,
      input_bam = input_bam,
      output_vcf = output_vcf,
      known_sites = known_sites,
      emit_reference_confidence = emit_reference_confidence,
      gatk_docker = gatk_docker,
      gcloud_service_account_key_file = gcloud_service_account_key_file,
      gcloud_project = gcloud_project,
      gatk = gatk,
      gatk_spark_jar = gatk_spark_jar,
      cluster_name = CreateDataprocCluster.resolved_cluster_name,
      gatk_gcs_staging = gatk_gcs_staging,
      executor_cores = executor_cores,
      executor_memory = executor_memory,
      driver_memory = driver_memory
  }

  if (select_first([delete_cluster, true])) {
    call DeleteDataprocCluster {
      input:
        gatk_docker = gatk_docker,
        gcloud_service_account_key_file = gcloud_service_account_key_file,
        gcloud_project = gcloud_project,
        cluster_name = CreateDataprocCluster.resolved_cluster_name,
        job_done = ReadsPipelineSpark.done
    }
  }
}

task ReadsPipelineSpark {
  String ref_fasta
  String input_bam
  String output_vcf
  String known_sites
  String emit_reference_confidence

  String gatk_docker
  File? gcloud_service_account_key_file
  String gcloud_project
  String gatk
  String gatk_spark_jar
  String cluster_name
  String? gatk_gcs_staging
  Int executor_cores
  String executor_memory
  String driver_memory

  command {
    if [ -n "${gcloud_service_account_key_file}" ]; then
      gcloud auth activate-service-account --key-file=${gcloud_service_account_key_file}
      export GOOGLE_APPLICATION_CREDENTIALS=${gcloud_service_account_key_file}
    fi
    gcloud config set project ${gcloud_project}

    export GATK_SPARK_JAR=${gatk_spark_jar}
    if [ -n "${gatk_gcs_staging}" ]; then
      export GATK_GCS_STAGING=${gatk_gcs_staging}
    fi

    ${gatk} \
      ReadsPipelineSpark \
      -R ${ref_fasta} \
      -I ${input_bam} \
      -O ${output_vcf} \
      --known-sites ${known_sites} \
      -pairHMM AVX_LOGLESS_CACHING \
      --max-reads-per-alignment-start 10 \
      -ERC ${emit_reference_confidence} \
      -- \
      --spark-runner GCS \
      --cluster ${cluster_name} \
      --executor-cores ${executor_cores} \
      --executor-memory ${executor_memory} \
      --driver-memory ${driver_memory}
  }

  runtime {
    docker: gatk_docker
  }

  output {
    String raw_vcf = "${output_vcf}"
    Boolean done = true
  }
}

task CreateDataprocCluster {
  String gatk_docker
  File? gcloud_service_account_key_file
  String gcloud_project
  String gcloud_zone = "us-central1-a"
  String? cluster_name
  Int num_workers = 10
  String master_machine_type = "n1-standard-4"
  Int master_boot_disk_size = 1000
  String worker_machine_type = "n1-standard-16"
  Int worker_boot_disk_size = 2000
  String image_version = "1.3"
  String max_age = "3h"

  command {
    if [ -n "${gcloud_service_account_key_file}" ]; then
      gcloud auth activate-service-account --key-file=${gcloud_service_account_key_file}
      export GOOGLE_APPLICATION_CREDENTIALS=${gcloud_service_account_key_file}
    fi
    gcloud config set project ${gcloud_project}

    if [ -z "${cluster_name}" ]; then
        # generate a unique cluster name
        cluster_name="gatk-$(python -c "import uuid; print(uuid.uuid4())")"
    else
        cluster_name=${cluster_name}
    fi
    if ! gcloud dataproc clusters describe "$cluster_name"; then
      # Use beta to make use of scheduled deletion (--max-age)
      # https://cloud.google.com/dataproc/docs/concepts/configuring-clusters/scheduled-deletion
      gcloud beta dataproc clusters create "$cluster_name" \
        --zone ${gcloud_zone} \
        --master-machine-type ${master_machine_type} \
        --master-boot-disk-size ${master_boot_disk_size} \
        --num-workers ${num_workers} \
        --worker-machine-type ${worker_machine_type} \
        --worker-boot-disk-size ${worker_boot_disk_size} \
        --image-version ${image_version} \
        --max-age ${max_age} \
        --project ${gcloud_project}
    fi
    echo "$cluster_name" > cluster_name.txt
  }

  runtime {
    docker: gatk_docker
  }

  output {
    String resolved_cluster_name = read_string("cluster_name.txt")
  }
}

task DeleteDataprocCluster {
  Boolean job_done
  String gatk_docker
  File? gcloud_service_account_key_file
  String gcloud_project
  String cluster_name

  command {
    if [ -n "${gcloud_service_account_key_file}" ]; then
      gcloud auth activate-service-account --key-file=${gcloud_service_account_key_file}
      export GOOGLE_APPLICATION_CREDENTIALS=${gcloud_service_account_key_file}
    fi
    gcloud config set project ${gcloud_project}

    if gcloud dataproc clusters describe "${cluster_name}"; then
      gcloud dataproc clusters delete --quiet "${cluster_name}"
    fi
  }

  runtime {
    docker: gatk_docker
  }

  output {
    Boolean done = true
  }
}