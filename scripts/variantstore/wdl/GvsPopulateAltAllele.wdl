version 1.0

import "GvsUtils.wdl" as Utils

workflow GvsPopulateAltAllele {
  input {
    Boolean go = true
    String dataset_name
    String project_id
    String call_set_identifier
  }

  String fq_alt_allele_table = "~{project_id}.~{dataset_name}.alt_allele"

  call CreateAltAlleleTable {
    input:
      dataset_name = dataset_name,
      project_id = project_id
  }

  call GetMaxSampleId {
    input:
      go = CreateAltAlleleTable.done,
      dataset_name = dataset_name,
      project_id = project_id
  }

  call GetVetTableNames {
    input:
      dataset_name = dataset_name,
      project_id = project_id,
      max_sample_id = GetMaxSampleId.max_sample_id
  }

  call Utils.GetBQTableLastModifiedDatetime {
    input:
      go = CreateAltAlleleTable.done,
      query_project = project_id,
      fq_table = fq_alt_allele_table
  }

  scatter (idx in range(length(GetVetTableNames.vet_tables))) {
    call PopulateAltAlleleTable {
      input:
        call_set_identifier = call_set_identifier,
        dataset_name = dataset_name,
        project_id = project_id,
        create_table_done = CreateAltAlleleTable.done,
        vet_table_name = GetVetTableNames.vet_tables[idx],
        last_modified_timestamp = GetBQTableLastModifiedDatetime.last_modified_timestamp,
        max_sample_id = GetMaxSampleId.max_sample_id
    }
  }

  output {
    Array[String] vet_tables_loaded = PopulateAltAlleleTable.done
    Boolean done = true
  }
}

task GetMaxSampleId {
  input {
    Boolean go = true
    String dataset_name
    String project_id
  }
  meta {
    # because this is being used to determine the current state of the GVS database, never use call cache
    volatile: true
  }

  command <<<
    set -o errexit -o nounset -o xtrace -o pipefail

    echo "project_id = ~{project_id}" > ~/.bigqueryrc
    bq query --location=US --project_id=~{project_id} --format=csv --use_legacy_sql=false \
    'SELECT IFNULL(MAX(sample_id), 0) AS max_sample_id FROM `~{dataset_name}.alt_allele`' > num_rows.csv

    # remove the header row from the CSV file
    sed -i 1d num_rows.csv
  >>>

  output {
    Int max_sample_id = read_int("num_rows.csv")
  }

  runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:398.0.0"
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
  }

  meta {
    # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
  }

  # add labels for DSP Cloud Cost Control Labeling and Reporting
  String bq_labels = "--label service:gvs --label team:variants --label managedby:create_alt_allele"

  command <<<
    set -e
    echo "project_id = ~{project_id}" > ~/.bigqueryrc

    # if the maximum sample_id value is evenly divisible by 4000, then max_sample_id / 4000 will
    # give us the right vet_* table to start with; otherwise, we need to start with the next table
    if [ $((~{max_sample_id} % 4000)) -eq 0 ]; then
      min_vat_table_num=$((~{max_sample_id} / 4000))
    else
      min_vat_table_num=$(((~{max_sample_id} / 4000) + 1))
    fi

    # use the number calculated from the above math to get the vet_* table names to grab data from
    bq query --location=US --project_id=~{project_id} --format=csv --use_legacy_sql=false ~{bq_labels} \
    "SELECT table_name FROM \`~{project_id}.~{dataset_name}.INFORMATION_SCHEMA.TABLES\` WHERE table_name LIKE 'vet_%' AND CAST(SUBSTRING(table_name, length('vet_') + 1) AS INT64) >= ${min_vat_table_num}" > vet_tables.csv

    # remove the header row from the CSV file
    sed -i 1d vet_tables.csv
  >>>
  runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:398.0.0"
    memory: "3 GB"
    disks: "local-disk 10 HDD"
    preemptible: 3
    cpu: 1
  }

  output {
    Array[String] vet_tables = read_lines("vet_tables.csv")
  }
}

task CreateAltAlleleTable {
  input {
    Boolean go = true
    String dataset_name
    String project_id
  }
  meta {
    # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
  }

  # add labels for DSP Cloud Cost Control Labeling and Reporting
  String bq_labels = "--label service:gvs --label team:variants --label managedby:create_alt_allele"

  command <<<
    set -e

    echo "project_id = ~{project_id}" > ~/.bigqueryrc
    bq query --location=US --project_id=~{project_id} --format=csv --use_legacy_sql=false ~{bq_labels} \
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
      ad INT64
    ) PARTITION BY RANGE_BUCKET(location, GENERATE_ARRAY(0, 25000000000000, 1000000000000))
    CLUSTER BY location, sample_id;'

  >>>
  runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:398.0.0"
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
    String dataset_name
    String project_id

    String create_table_done
    String vet_table_name
    String call_set_identifier
    Int max_sample_id

    String last_modified_timestamp
  }
  meta {
    # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
  }

  command <<<
    set -e

    python3 /app/populate_alt_allele_table.py \
      --call_set_identifier ~{call_set_identifier} \
      --query_project ~{project_id} \
      --vet_table_name ~{vet_table_name} \
      --fq_dataset ~{project_id}.~{dataset_name} \
      --max_sample_id ~{max_sample_id} \
  >>>
  runtime {
    docker: "us.gcr.io/broad-dsde-methods/variantstore:ah_var_store_2022_08_22"
    memory: "3 GB"
    disks: "local-disk 10 HDD"
    cpu: 1
  }

  output {
    String done = "~{vet_table_name}"
  }
}
