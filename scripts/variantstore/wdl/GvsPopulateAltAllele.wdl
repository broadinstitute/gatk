version 1.0

import "GvsUtils.wdl" as Utils

workflow GvsPopulateAltAllele {
  input {
    # Intentionally unused: this input exists solely to enforce task ordering - the upstream task's `done` output
    # is passed here to prevent this task from running until the upstream task has completed.
    #@ except: UnusedInput
    Boolean go = true
    String dataset_name
    String project_id
    String call_set_identifier
    Int max_alt_allele_shards = 10
    String? variants_docker
    String? cloud_sdk_docker
    String? git_branch_or_tag
    String? git_hash
  }

  String fq_alt_allele_table = "~{project_id}.~{dataset_name}.alt_allele"

  if (!defined(git_hash) || !defined(cloud_sdk_docker) || !defined(variants_docker)) {
    call Utils.GetToolVersions {
      input:
        git_branch_or_tag = git_branch_or_tag,
    }
  }

  String effective_cloud_sdk_docker = select_first([cloud_sdk_docker, GetToolVersions.cloud_sdk_docker])
  String effective_variants_docker = select_first([variants_docker, GetToolVersions.variants_docker])
  String effective_git_hash = select_first([git_hash, GetToolVersions.git_hash])

  # With enable_vrs_ids unset, this reports whether the dataset carries VRS IDs, from its vet tables, and fails if
  # they disagree. alt_allele follows the vet tables: it has vrs_allele_id only when they have vrs_allele_ids.
  call Utils.ResolveEnableVrsIds as ResolveVrsIds {
    input:
      project_id = project_id,
      dataset_name = dataset_name,
      cloud_sdk_docker = effective_cloud_sdk_docker,
  }
  Boolean dataset_has_vrs_ids = ResolveVrsIds.resolved_enable_vrs_ids

  call CreateAltAlleleTable {
    input:
      dataset_name = dataset_name,
      project_id = project_id,
      vrs_ids = dataset_has_vrs_ids,
      cloud_sdk_docker = effective_cloud_sdk_docker,
  }

  call GetMaxSampleId {
    input:
      go = CreateAltAlleleTable.done,
      dataset_name = dataset_name,
      project_id = project_id,
      cloud_sdk_docker = effective_cloud_sdk_docker,
  }

  call GetVetTableNames {
    input:
      dataset_name = dataset_name,
      project_id = project_id,
      max_sample_id = GetMaxSampleId.max_sample_id,
      max_alt_allele_shards = max_alt_allele_shards,
      cloud_sdk_docker = effective_cloud_sdk_docker,
  }

  call Utils.GetBQTableLastModifiedDatetime {
    input:
      go = CreateAltAlleleTable.done,
      project_id = project_id,
      fq_table = fq_alt_allele_table,
      cloud_sdk_docker = effective_cloud_sdk_docker,
  }

  scatter (vet_table_names_file in GetVetTableNames.vet_table_names_files) {
    call PopulateAltAlleleTable {
      input:
        call_set_identifier = call_set_identifier,
        dataset_name = dataset_name,
        project_id = project_id,
        go = CreateAltAlleleTable.done,
        vet_table_names_file = vet_table_names_file,
        last_modified_timestamp = GetBQTableLastModifiedDatetime.last_modified_timestamp,
        max_sample_id = GetMaxSampleId.max_sample_id,
        vrs_ids = dataset_has_vrs_ids,
        variants_docker = effective_variants_docker,
    }
  }

  call VerifySampleCount {
    input:
      go = PopulateAltAlleleTable.done[0],
      dataset_name = dataset_name,
      project_id = project_id,
      cloud_sdk_docker = effective_cloud_sdk_docker,
  }

  output {
    Boolean done = VerifySampleCount.done
    String recorded_git_hash = effective_git_hash
  }
}

task GetMaxSampleId {
  input {
    # Intentionally unused: this input exists solely to enforce task ordering - the upstream task's `done` output
    # is passed here to prevent this task from running until the upstream task has completed.
    #@ except: UnusedInput
    Boolean go = true
    String dataset_name
    String project_id
    String cloud_sdk_docker
  }
  meta {
    # because this is being used to determine the current state of the GVS database, never use call cache
    volatile: true
  }

  command <<<
    # Prepend date, time and pwd to xtrace log entries.
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o pipefail -o xtrace

    echo "project_id = ~{project_id}" > ~/.bigqueryrc
    # bq query --max_rows check: ok one row
    bq --apilog=false query --project_id=~{project_id} --format=csv --use_legacy_sql=false \
    'SELECT IFNULL(MAX(sample_id), 0) AS max_sample_id FROM `~{dataset_name}.alt_allele`' > num_rows.csv

    # remove the header row from the CSV file
    sed -i 1d num_rows.csv
  >>>

  output {
    Int max_sample_id = read_int("num_rows.csv")
  }

  runtime {
    docker: cloud_sdk_docker
    memory: "3 GB"
    disks: "local-disk 10 HDD"
    preemptible: 3
    cpu: 1
  }
}

task GetVetTableNames {
  input {
    String dataset_name
    String project_id
    Int max_sample_id
    Int max_alt_allele_shards
    String cloud_sdk_docker
  }

  meta {
    # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
  }

  # add labels for DSP Cloud Cost Control Labeling and Reporting
  String bq_labels = "--label service:gvs --label team:variants --label managedby:populate_alt_allele"

  command <<<
    # Prepend date, time and pwd to xtrace log entries.
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o pipefail -o xtrace

    echo "project_id = ~{project_id}" > ~/.bigqueryrc

    # if the maximum sample_id value is evenly divisible by 4000, then max_sample_id / 4000 will
    # give us the right vet_* table to start with; otherwise, we need to start with the next table
    if [ $((~{max_sample_id} % 4000)) -eq 0 ]; then
      min_vat_table_num=$((~{max_sample_id} / 4000))
    else
      min_vat_table_num=$(((~{max_sample_id} / 4000) + 1))
    fi

    # use the number calculated from the above math to get the vet_* table names to grab data from
    # set max rows to at least the number of superpartitions (Math.ceil(num_samples / 4000))
    bq --apilog=false query --max_rows 10000000 --project_id=~{project_id} --format=csv --use_legacy_sql=false ~{bq_labels} \
      'SELECT
          table_name
        FROM `~{project_id}.~{dataset_name}.INFORMATION_SCHEMA.TABLES`
        WHERE
          REGEXP_CONTAINS(table_name, "^vet_[0-9]+$") AND
          CAST(SUBSTRING(table_name, length("vet_") + 1) AS INT64) >= '"${min_vat_table_num}" > vet_tables.csv

    # remove the header row from the CSV file, count the number of tables and divide them up into
    # no more than max_alt_allele_shards files (which is the estimated number of BQ queries we can run at once without
    # reserving flex slots)
    sed -i 1d vet_tables.csv
    num_tables=$(cat vet_tables.csv | wc -l)
    if [ $((num_tables % ~{max_alt_allele_shards})) -eq 0 ]; then
      num_tables_per_file=$((num_tables / ~{max_alt_allele_shards}))
    else
      num_tables_per_file=$(((num_tables / ~{max_alt_allele_shards}) + 1))
    fi
    split -l $num_tables_per_file vet_tables.csv vet_tables_
  >>>
  runtime {
    docker: cloud_sdk_docker
    memory: "3 GB"
    disks: "local-disk 10 HDD"
    preemptible: 3
    cpu: 1
  }

  output {
    Array[File] vet_table_names_files = glob("vet_tables_*")
  }
}

task CreateAltAlleleTable {
  input {
    # Intentionally unused: this input exists solely to enforce task ordering - the upstream task's `done` output
    # is passed here to prevent this task from running until the upstream task has completed.
    #@ except: UnusedInput
    Boolean go = true
    String dataset_name
    String project_id
    # Whether the dataset carries VRS IDs, as ResolveEnableVrsIds reads it from the vet tables.
    Boolean vrs_ids
    String cloud_sdk_docker
  }
  meta {
    # should always be run; if the table already exists, no harm no foul
    volatile: true
  }

  # add labels for DSP Cloud Cost Control Labeling and Reporting
  String bq_labels = "--label service:gvs --label team:variants --label managedby:create_alt_allele"

  command <<<
    # Prepend date, time and pwd to xtrace log entries.
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o pipefail -o xtrace

    echo "project_id = ~{project_id}" > ~/.bigqueryrc

    # alt_allele has the vrs_allele_id column exactly when the vet tables have vrs_allele_ids. An existing alt_allele
    # that disagrees cannot be fixed by adding or ignoring the column, so fail. alt_allele is derived entirely from the
    # vet tables, so dropping it and rerunning rebuilds it to match.
    bq --apilog=false query --project_id=~{project_id} --format=csv --use_legacy_sql=false --max_rows=1 ~{bq_labels} '
      SELECT
        COUNT(*) > 0 AS table_exists,
        COUNTIF(column_name = "vrs_allele_id") > 0 AS has_vrs_allele_id
      FROM `~{project_id}.~{dataset_name}.INFORMATION_SCHEMA.COLUMNS`
      WHERE table_name = "alt_allele"' | sed 1d > alt_allele_columns.csv

    table_exists=$(cut -d, -f1 alt_allele_columns.csv)
    has_vrs_allele_id=$(cut -d, -f2 alt_allele_columns.csv)
    vrs_ids="~{vrs_ids}"
    where="~{project_id}:~{dataset_name}"

    if [ "${table_exists}" = "true" ] && [ "${has_vrs_allele_id}" != "${vrs_ids}" ]; then
      # Turn off xtrace so the message is easy to read in the task log.
      set +o xtrace
      if [ "${vrs_ids}" = "true" ]; then
        echo "ERROR: the vet tables in ${where} have the vrs_allele_ids column, so the dataset carries VRS IDs," >&2
        echo "but its alt_allele table has no vrs_allele_id column. It was built before the dataset had VRS IDs." >&2
      else
        echo "ERROR: the vet tables in ${where} have no vrs_allele_ids column, so the dataset does not carry" >&2
        echo "VRS IDs, but its alt_allele table has a vrs_allele_id column." >&2
      fi
      echo "alt_allele is built entirely from the vet tables. Drop it with" >&2
      echo "    bq rm -f -t ${where}.alt_allele" >&2
      echo "and rerun GvsPopulateAltAllele to rebuild it to match them." >&2
      exit 1
    fi

    # bq query --max_rows check: ok create table
    bq --apilog=false query --project_id=~{project_id} --format=csv --use_legacy_sql=false ~{bq_labels} \
    'CREATE TABLE IF NOT EXISTS `~{project_id}.~{dataset_name}.alt_allele` (
      location INT64,
      sample_id INT64,
      ref STRING,
      allele STRING,
      allele_pos INT64,
      call_GT STRING,
      call_GQ INT64,
      as_raw_mq STRING,
      raw_mq INT64,
      as_raw_mqranksum STRING,
      raw_mqranksum_x_10 INT64,
      as_qualapprox STRING,
      qualapprox STRING,
      qual INT64,
      as_raw_readposranksum STRING,
      raw_readposranksum_x_10 INT64,
      as_sb_table STRING,
      sb_ref_plus INT64,
      sb_ref_minus INT64,
      sb_alt_plus INT64,
      sb_alt_minus INT64,
      call_AD STRING,
      ref_ad INT64,
      ad INT64~{if vrs_ids then ", vrs_allele_id STRING" else ""}
    ) PARTITION BY RANGE_BUCKET(location, GENERATE_ARRAY(0, 25000000000000, 1000000000000))
    CLUSTER BY location, sample_id;'

  >>>
  runtime {
    docker: cloud_sdk_docker
    memory: "3 GB"
    disks: "local-disk 10 HDD"
    cpu: 1
  }

  output {
    Boolean done = true
  }
}

task PopulateAltAlleleTable {
  input {
    # Intentionally unused: this input exists solely to enforce task ordering - the upstream task's `done` output
    # is passed here to prevent this task from running until the upstream task has completed.
    #@ except: UnusedInput
    Boolean go

    String dataset_name
    String project_id

    File vet_table_names_file
    String call_set_identifier
    Int max_sample_id
    # Intentionally unused: passed solely to bust WDL call-caching when the referenced BigQuery table has been modified.
    #@ except: UnusedInput
    String last_modified_timestamp
    # Whether the dataset carries VRS IDs; if so, alt_allele has vrs_allele_id and the script fills it.
    Boolean vrs_ids

    String variants_docker
  }
  meta {
    # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
  }

  Array[String] vet_table_names = read_lines(vet_table_names_file)

  command <<<
    # Prepend date, time and pwd to xtrace log entries.
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o pipefail -o xtrace

    VET_TABLES_ARRAY=(~{sep=" " vet_table_names})

    for vet_table in "${VET_TABLES_ARRAY[@]}"; do
      python3 /app/populate_alt_allele_table.py \
        --call_set_identifier ~{call_set_identifier} \
        --query_project ~{project_id} \
        --vet_table_name ${vet_table} \
        --fq_dataset ~{project_id}.~{dataset_name} \
        --max_sample_id ~{max_sample_id} ~{if vrs_ids then "--include_vrs_allele_id" else ""}
    done
  >>>
  runtime {
    docker: variants_docker
    memory: "3 GB"
    disks: "local-disk 10 HDD"
    cpu: 1
    noAddress: true
  }

  output {
    Boolean done = true
  }
}

task VerifySampleCount {
  input {
    # Intentionally unused: this input exists solely to enforce task ordering - the upstream task's `done` output
    # is passed here to prevent this task from running until the upstream task has completed.
    #@ except: UnusedInput
    Boolean go = true
    String dataset_name
    String project_id
    String cloud_sdk_docker
  }
  meta {
    # because this is being used to determine the current state of the GVS database, never use call cache
    volatile: true
  }

  command <<<
    # Prepend date, time and pwd to xtrace log entries.
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o pipefail -o xtrace

    # Determine the number of distinct sample_ids in the alt_allele table
    bq --apilog=false query --project_id=~{project_id} --format=csv --use_legacy_sql=false \
      'SELECT COUNT(DISTINCT sample_id) AS num_samples FROM `~{dataset_name}.alt_allele`' | sed 1d > num_samples_in_aa.csv

    # Determine the number of distinct sample_ids in the sample_info table
    bq --apilog=false query --project_id=~{project_id} --format=csv --use_legacy_sql=false \
    'SELECT COUNT(DISTINCT sample_id) AS num_samples FROM `~{dataset_name}.sample_info` WHERE withdrawn IS NULL' | sed 1d > num_samples_in_si.csv

    # Compare the two counts
    set +o errexit
    diff num_samples_in_aa.csv num_samples_in_si.csv
    rc=$?
    set -o errexit
    if [[ $rc -ne 0 ]]; then
      # Just log this, don't fail the task. Because we may withdraw samples at any time, it's reasonable that the counts may not match.
      echo "Note: The number of samples (`cat num_samples_in_aa.csv`) in ~{dataset_name}.alt_allele differs from that (`cat num_samples_in_si.csv`) in ~{dataset_name}.sample_info!"
    fi
  >>>

  output {
    Boolean done = true
  }

  runtime {
    docker: cloud_sdk_docker
    memory: "3 GB"
    disks: "local-disk 10 HDD"
    preemptible: 3
    cpu: 1
  }
}
