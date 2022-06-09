version 1.0

task MergeVCFs {
  input {
    Array[File] input_vcfs
    String gather_type = "BLOCK"
    String output_vcf_name
    String? output_directory
    Int? merge_disk_override
    Int? preemptible_tries
    String? service_account_json_path
  }

  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'
  Int disk_size = if (defined(merge_disk_override)) then merge_disk_override else 100

  parameter_meta {
    input_vcfs: {
      localization_optional: true
    }
  }

  command {
    if [ ~{has_service_account_file} = 'true' ]; then
      gsutil cp ~{service_account_json_path} local.service_account.json
      gcloud auth activate-service-account --key-file=local.service_account.json
      export GOOGLE_APPLICATION_CREDENTIALS=local.service_account.json
    fi

    gatk --java-options -Xmx3g GatherVcfsCloud \
    --ignore-safety-checks --gather-type ~{gather_type} \
    --create-output-variant-index false \
    -I ~{sep=' -I ' input_vcfs} \
    --output ~{output_vcf_name}

    tabix ~{output_vcf_name}

    # Drop trailing slash if one exists
    OUTPUT_GCS_DIR=$(echo ~{output_directory} | sed 's/\/$//')

    if [ -n "$OUTPUT_GCS_DIR" ]; then
      gsutil cp ~{output_vcf_name} $OUTPUT_GCS_DIR/
      gsutil cp ~{output_vcf_name}.tbi $OUTPUT_GCS_DIR/
    fi
  }

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:varstore_d8a72b825eab2d979c8877448c0ca948fd9b34c7_change_to_hwe"
    preemptible: select_first([preemptible_tries, 3])
    memory: "3 GiB"
    disks: "local-disk ~{disk_size} HDD"
  }

  output {
    File output_vcf = "~{output_vcf_name}"
    File output_vcf_index = "~{output_vcf_name}.tbi"
  }
}

task SplitIntervals {
  input {
    File intervals
    File ref_fasta
    File ref_fai
    File ref_dict
    Int scatter_count
    File? interval_weights_bed
    String? intervals_file_extension
    String? split_intervals_extra_args
    Int? split_intervals_disk_size_override
    Int? split_intervals_mem_override
    String? output_gcs_dir
    File? gatk_override
    String? service_account_json_path
  }
  meta {
    # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
  }

  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

  Int disk_size = if (defined(split_intervals_disk_size_override)) then split_intervals_disk_size_override else 10
  Int disk_memory = if (defined(split_intervals_mem_override)) then split_intervals_mem_override else 16
  Int java_memory = disk_memory - 4

  String gatkTool = if (defined(interval_weights_bed)) then 'WeightedSplitIntervals' else 'SplitIntervals'

  parameter_meta {
    intervals: {
      localization_optional: true
    }
    ref_fasta: {
      localization_optional: true
    }
    ref_fai: {
      localization_optional: true
    }
    ref_dict: {
      localization_optional: true
    }
  }

  command {
    set -e
    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

    mkdir interval-files
    gatk --java-options "-Xmx~{java_memory}g" ~{gatkTool} \
      --dont-mix-contigs \
      -R ~{ref_fasta} \
      ~{"-L " + intervals} \
      ~{"--weight-bed-file " + interval_weights_bed} \
      -scatter ~{scatter_count} \
      -O interval-files \
      ~{"--extension " + intervals_file_extension} \
      --interval-file-num-digits 10 \
      ~{split_intervals_extra_args}
    cp interval-files/*.interval_list .

    # Drop trailing slash if one exists
    OUTPUT_GCS_DIR=$(echo ~{output_gcs_dir} | sed 's/\/$//')

    if [ ~{has_service_account_file} = 'true' ]; then
      gsutil cp ~{service_account_json_path} local.service_account.json
      gcloud auth activate-service-account --key-file=local.service_account.json
    fi

    if [ -n "$OUTPUT_GCS_DIR" ]; then
      gsutil -m cp *.interval_list $OUTPUT_GCS_DIR/
    fi
  }

  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:4.2.3.0"
    bootDiskSizeGb: 15
    memory: "~{disk_memory} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: 3
    cpu: 1
  }

  output {
    Array[File] interval_files = glob("*.interval_list")
  }
}

task GetBQTableLastModifiedDatetime {
  input {
    String query_project
    String fq_table
    String? service_account_json_path
  }
  meta {
    # because this is being used to determine if the data has changed, never use call cache
    volatile: true
  }

  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

  # ------------------------------------------------
  # try to get the last modified date for the table in question; fail if something comes back from BigQuery
  # that isn't in the right format (e.g. an error)
  command <<<
    set -o xtrace
    set -o errexit

    if [ ~{has_service_account_file} = 'true' ]; then
      gsutil cp ~{service_account_json_path} local.service_account.json
      gcloud auth activate-service-account --key-file=local.service_account.json
      gcloud config set project ~{query_project}
    fi

    echo "project_id = ~{query_project}" > ~/.bigqueryrc

    # bq needs the project name to be separate by a colon
    DATASET_TABLE_COLON=$(echo ~{fq_table} | sed 's/\./:/')

    LASTMODIFIED=$(bq --location=US --project_id=~{query_project} --format=json show ${DATASET_TABLE_COLON} | python3 -c "import sys, json; print(json.load(sys.stdin)['lastModifiedTime']);")
    if [[ $LASTMODIFIED =~ ^[0-9]+$ ]]; then
      echo $LASTMODIFIED
    else
      exit 1
    fi
  >>>

  output {
    String last_modified_timestamp = read_string(stdout())
  }

  runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
    memory: "3 GB"
    disks: "local-disk 10 HDD"
    preemptible: 3
    cpu: 1
  }
}

task GetBQTablesMaxLastModifiedTimestamp {
  input {
    String query_project
    String data_project
    String data_dataset
    Array[String] table_patterns
    String? service_account_json_path
  }
  meta {
    # because this is being used to determine if the data has changed, never use call cache
    volatile: true
  }

  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

  # ------------------------------------------------
  # try to get the latest last modified timestamp, in epoch microseconds, for all of the tables that match the provided prefixes
  command <<<
    set -e
    if [ ~{has_service_account_file} = 'true' ]; then
    gsutil cp ~{service_account_json_path} local.service_account.json
    gcloud auth activate-service-account --key-file=local.service_account.json
    fi

    echo "project_id = ~{query_project}" > ~/.bigqueryrc

    bq --location=US --project_id=~{query_project} query --format=csv --use_legacy_sql=false \
    "SELECT UNIX_MICROS(MAX(last_modified_time)) last_modified_time FROM \`~{data_project}\`.~{data_dataset}.INFORMATION_SCHEMA.PARTITIONS WHERE table_name like '~{sep="' OR table_name like '" table_patterns}'" > results.txt

    tail -1 results.txt | cut -d, -f1 > max_last_modified_timestamp.txt
  >>>

  output {
    String max_last_modified_timestamp = read_string("max_last_modified_timestamp.txt")
  }

  runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
    memory: "3 GB"
    disks: "local-disk 10 HDD"
    preemptible: 3
    cpu: 1
  }
}

task BuildGATKJarAndCreateDataset {
  input {
    String branch_name
    String dataset_prefix
  }
  meta {
    # Branch may be updated so do not call cache!
    volatile: true
  }

  command <<<
    # Much of this could/should be put into a Docker image.
    set -o errexit -o nounset -o pipefail

    # git and git-lfs
    apt-get -qq update
    apt-get -qq install git git-lfs

    # Java
    apt-get -qq install wget apt-transport-https gnupg
    wget -O - https://packages.adoptium.net/artifactory/api/gpg/key/public | apt-key add -
    echo "deb https://packages.adoptium.net/artifactory/deb $(awk -F= '/^VERSION_CODENAME/{print$2}' /etc/os-release) main" | tee /etc/apt/sources.list.d/adoptium.list
    apt-get -qq update
    apt -qq install -y temurin-11-jdk

    # GATK
    git clone https://github.com/broadinstitute/gatk.git --depth 1 --branch ~{branch_name} --single-branch
    cd gatk
    ./gradlew shadowJar

    branch=$(git symbolic-ref HEAD 2>/dev/null)
    branch=${branch#refs/heads/}

    hash=$(git rev-parse --short HEAD)

    # Rename the GATK jar to embed the branch and hash of the most recent commit on the branch.
    mv build/libs/gatk-package-unspecified-SNAPSHOT-local.jar "build/libs/gatk-${branch}-${hash}-SNAPSHOT-local.jar"

    # Build a dataset name based on the branch name and the git hash of the most recent commit on this branch.
    # Dataset names must be alphanumeric and underscores only. Convert any dashes to underscores, then delete
    # any remaining characters that are not alphanumeric or underscores.
    dataset="$(echo ~{dataset_prefix}_${branch}_${hash} | tr '-' '_' | tr -c -d '[:alnum:]_')"

    bq mk --project_id="spec-ops-aou" "$dataset"

    echo -n "$dataset" > dataset.txt
  >>>

  output {
    Boolean done = true
    File jar = glob("gatk/build/libs/*-SNAPSHOT-local.jar")[0]
    String dataset_name = read_string("gatk/dataset.txt")
  }

  runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:latest"
    disks: "local-disk 500 HDD"
  }
}

task TerminateWorkflow {
  input {
    String message
  }
  meta {
    # Definitely do not call cache this!
    volatile: true
  }

  command <<<
    set -o errexit

    # To avoid issues with special characters within the message, write the message to a file.
    cat > message.txt <<FIN
    ~{message}
    FIN

    # cat the file to stderr as this task is going to fail due to the exit 1.
    cat message.txt >&2
    exit 1
  >>>

  runtime {
    docker: "python:3.8-slim-buster"
    memory: "1 GB"
    disks: "local-disk 10 HDD"
    preemptible: 3
    cpu: 1
  }

  output {
    Boolean done = true
  }
}

task ScaleXYBedValues {
    input {
        Boolean go = true
        File interval_weights_bed
        Float x_bed_weight_scaling
        Float y_bed_weight_scaling
    }
    meta {
        # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
    }
    command <<<
        python3 /app/scale_xy_bed_values.py \
            --input ~{interval_weights_bed} \
            --output "interval_weights_xy_scaled.bed" \
            --xscale ~{x_bed_weight_scaling} \
            --yscale ~{y_bed_weight_scaling} \
    >>>

    output {
        File xy_scaled_bed = "interval_weights_xy_scaled.bed"
        Boolean done = true
    }

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:vs_389_xy_reweighting_2022_05_24"
        maxRetries: 3
        memory: "7 GB"
        preemptible: 3
        cpu: "2"
        disks: "local-disk 500 HDD"
    }
}

task GetNumSamplesLoaded {
  input {
    String fq_sample_table
    String fq_sample_table_lastmodified_timestamp
    String? service_account_json_path
    String project_id
    Boolean control_samples = false
  }
  meta {
    # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
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
    bq query --location=US --project_id=~{project_id} --format=csv --use_legacy_sql=false \
      'SELECT COUNT(*) as num_rows FROM `~{fq_sample_table}` WHERE is_loaded = true and is_control = ~{control_samples}' > num_rows.csv

    NUMROWS=$(python3 -c "csvObj=open('num_rows.csv','r');csvContents=csvObj.read();print(csvContents.split('\n')[1]);")

    [[ $NUMROWS =~ ^[0-9]+$ ]] && echo $NUMROWS || exit 1
  >>>

  output {
    Int num_samples = read_int(stdout())
  }

  runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
    memory: "3 GB"
    disks: "local-disk 10 HDD"
    preemptible: 3
    cpu: 1
  }
}
