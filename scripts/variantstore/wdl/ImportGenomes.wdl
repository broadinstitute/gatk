version 1.0

import "https://api.firecloud.org/ga4gh/v1/tools/synthetic-microarray-gen:LoadBigQueryData/versions/7/plain-WDL/descriptor" as LoadBigQueryData

workflow ImportGenomes {

  input {
    File input_vcfs_list
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
  Array[String] input_vcfs = read_lines(input_vcfs_list)

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
        drop_state_includes_greater_than = drop_state_includes_greater_than,
        output_directory = output_directory,
        gatk_override = gatk_override,
        docker = docker_final,
        preemptible_tries = preemptible_tries
    }
  }

  call LoadBigQueryData.LoadBigQueryData as LoadMetadataTsvs {
    input:
      done = CreateImportTsvs.done,
      project_id = project_id,
      dataset_name = dataset_name,
      storage_location = output_directory,
      datatype = "metadata",
      numbered = "false",
      partitioned = "false",
      max_table_id = GetMaxTableId.max_table_id,
      schema = metadata_schema,
      preemptible_tries = preemptible_tries,
      docker = docker_final
  }

  call LoadBigQueryData.LoadBigQueryData as LoadPetTsvs {
    input:
      done = CreateImportTsvs.done,
      project_id = project_id,
      dataset_name = dataset_name,
      storage_location = output_directory,
      datatype = "pet",
      max_table_id = GetMaxTableId.max_table_id,
      schema = pet_schema,
      preemptible_tries = preemptible_tries,
      docker = docker_final
  }

  call LoadBigQueryData.LoadBigQueryData as LoadVetTsvs {
    input:
      done = CreateImportTsvs.done,
      project_id = project_id,
      dataset_name = dataset_name,
      storage_location = output_directory,
      datatype = "vet",
      max_table_id = GetMaxTableId.max_table_id,
      schema = vet_schema,
      preemptible_tries = preemptible_tries,
      docker = docker_final
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

    String? for_testing_only
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
      ~{for_testing_only}

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

