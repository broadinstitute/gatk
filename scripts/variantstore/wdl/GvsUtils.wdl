version 1.0

task MergeVCFs {
  input {
    Array[File] input_vcfs
    String gather_type = "BLOCK"
    String output_vcf_name
    String? output_directory
    Int? merge_disk_override
    Int? preemptible_tries
  }

  Int disk_size = if (defined(merge_disk_override)) then merge_disk_override else 100

  parameter_meta {
    input_vcfs: {
      localization_optional: true
    }
  }

  command {
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
    docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:varstore_2022_10_17_2a8c210ac35094997603259fa1cd784486b92e42"
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
  }
  meta {
    # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
  }

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

    if [ -n "$OUTPUT_GCS_DIR" ]; then
      gsutil -m cp *.interval_list $OUTPUT_GCS_DIR/
    fi
  }

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:varstore_2022_10_17_2a8c210ac35094997603259fa1cd784486b92e42"
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
    Boolean go = true
    String query_project
    String fq_table
  }
  meta {
    # because this is being used to determine if the data has changed, never use call cache
    volatile: true
  }

  # ------------------------------------------------
  # try to get the last modified date for the table in question; fail if something comes back from BigQuery
  # that isn't in the right format (e.g. an error)
  command <<<
    set -o xtrace
    set -o errexit

    echo "project_id = ~{query_project}" > ~/.bigqueryrc

    # bq needs the project name to be separate by a colon
    DATASET_TABLE_COLON=$(echo ~{fq_table} | sed 's/\./:/')

    LASTMODIFIED=$(bq --project_id=~{query_project} --format=json show ${DATASET_TABLE_COLON} | python3 -c "import sys, json; print(json.load(sys.stdin)['lastModifiedTime']);")
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
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:404.0.0-alpine"
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
  }
  meta {
    # because this is being used to determine if the data has changed, never use call cache
    volatile: true
  }

  # ------------------------------------------------
  # try to get the latest last modified timestamp, in epoch microseconds, for all of the tables that match the provided prefixes
  command <<<
    set -e

    echo "project_id = ~{query_project}" > ~/.bigqueryrc

    bq --project_id=~{query_project} query --format=csv --use_legacy_sql=false \
    "SELECT UNIX_MICROS(MAX(last_modified_time)) last_modified_time FROM \`~{data_project}\`.~{data_dataset}.INFORMATION_SCHEMA.PARTITIONS WHERE table_name like '~{sep="' OR table_name like '" table_patterns}'" > results.txt

    tail -1 results.txt | cut -d, -f1 > max_last_modified_timestamp.txt
  >>>

  output {
    String max_last_modified_timestamp = read_string("max_last_modified_timestamp.txt")
  }

  runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:404.0.0-alpine"
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
    # Prepend date, time and pwd to xtrace log entries.
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o pipefail -o xtrace

    # git and git-lfs
    apt-get -qq update
    apt-get -qq install git git-lfs

    # The Terra microservices are currently aligned on Temurin as their JDK distribution which is why we use it here.
    # However at least once Temurin unexpectedly became unavailable for download for several hours which broke this
    # task and its dependents. The following block switched the JDK distribution to Amazon Corretto 11 which appeared to
    # work just fine for our purposes during Temurin's brief absence.
    #
    # Corretto Java 11
    # apt-get -qq install wget apt-transport-https gnupg software-properties-common
    # wget -O- https://apt.corretto.aws/corretto.key | apt-key add -
    # add-apt-repository 'deb https://apt.corretto.aws stable main'

    # Temurin Java 11
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

    bq mk --project_id="gvs-internal" "$dataset"

    # add labels for DSP Cloud Cost Control Labeling and Reporting
    bq update --set_label service:gvs gvs-internal:$dataset
    bq update --set_label team:variants gvs-internal:$dataset
    bq update --set_label environment:dev gvs-internal:$dataset
    bq update --set_label managedby:build_gatk_jar_and_create_dataset gvs-internal:$dataset

    echo -n "$dataset" > dataset.txt
  >>>

  output {
    Boolean done = true
    File jar = glob("gatk/build/libs/*-SNAPSHOT-local.jar")[0]
    String dataset_name = read_string("gatk/dataset.txt")
  }

  runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:404.0.0-slim"
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
        docker: "us.gcr.io/broad-dsde-methods/variantstore:2022-11-14-v2-alpine"
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
    String project_id
    String sample_table_timestamp
    Boolean control_samples = false
  }
  meta {
    # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
  }

  command <<<
    set -o errexit -o nounset -o xtrace -o pipefail

    echo "project_id = ~{project_id}" > ~/.bigqueryrc
    bq query --project_id=~{project_id} --format=csv --use_legacy_sql=false '

      SELECT COUNT(*) FROM `~{fq_sample_table}` WHERE
        is_loaded = true AND
        withdrawn IS NULL AND
        is_control = ~{control_samples}

    ' | sed 1d
  >>>

  output {
    Int num_samples = read_int(stdout())
  }

  runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:404.0.0-alpine"
    memory: "3 GB"
    disks: "local-disk 10 HDD"
    preemptible: 3
    cpu: 1
  }
}


task CountSuperpartitions {
    meta {
        description: "Return the number of superpartitions based on the number of vet_% tables in `INFORMATION_SCHEMA.PARTITIONS`."
        # Definitely don't cache this, the values can change while the inputs to this task will not!
        volatile: true
    }
    input {
        String project_id
        String dataset_name
    }
    command <<<
        bq query --location=US --project_id='~{project_id}' --format=csv --use_legacy_sql=false '

            SELECT COUNT(*) FROM `~{project_id}.~{dataset_name}.INFORMATION_SCHEMA.TABLES`
                WHERE table_name LIKE "vet_%"

        ' | sed 1d > num_superpartitions.txt
    >>>
    runtime {
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:404.0.0-alpine"
        disks: "local-disk 500 HDD"
    }
    output {
        Int num_superpartitions = read_int('num_superpartitions.txt')
    }
}

task ValidateFilterSetName {
    input {
        Boolean go = true
        String filter_set_name
        String data_project
        String data_dataset
        String query_project = data_project
        String filter_set_info_timestamp = ""
    }
    meta {
        # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
    }

    # add labels for DSP Cloud Cost Control Labeling and Reporting
    String bq_labels = "--label service:gvs --label team:variants --label managedby:gvs_utils"

    command <<<
        set -o errexit -o nounset -o xtrace -o pipefail

        echo "project_id = ~{query_project}" > ~/.bigqueryrc

        OUTPUT=$(bq --project_id=~{query_project} --format=csv query --use_legacy_sql=false ~{bq_labels} "SELECT filter_set_name as available_filter_set_names FROM \`~{data_project}.~{data_dataset}.filter_set_info\` GROUP BY filter_set_name")
        FILTERSETS=${OUTPUT#"available_filter_set_names"}

        if [[ $FILTERSETS =~ "~{filter_set_name}" ]]; then
            echo "Filter set name '~{filter_set_name}' found."
        else
            echo "ERROR: '~{filter_set_name}' is not an existing filter_set_name. Available in ~{data_project}.~{data_dataset} are"
            echo $FILTERSETS
            exit 1
        fi
    >>>
    output {
        Boolean done = true
    }

    runtime {
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:404.0.0-alpine"
        memory: "3 GB"
        disks: "local-disk 500 HDD"
        preemptible: 3
        cpu: 1
    }
}
