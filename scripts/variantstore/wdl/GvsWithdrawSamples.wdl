version 1.0

workflow GvsWithdrawSamples {

  input {
    String dataset_name
    String project_id

    String sample_name_column_name_in_file
    File sample_names_to_include_file
    # should be in the format "2022-01-01 00:00:00 UTC"
    String withdrawn_timestamp
  }

  call WithdrawSamples {
    input:
      project_id = project_id,
      dataset_name = dataset_name,
      sample_name_column_name_in_file = sample_name_column_name_in_file,
      sample_names_to_include_file = sample_names_to_include_file,
      withdrawn_timestamp = withdrawn_timestamp
  }

  output {
    Int num_rows_updated = WithdrawSamples.num_rows_updated
  }
}

task WithdrawSamples {
  input {
    String project_id
    String dataset_name

    String sample_name_column_name_in_file
    File sample_names_to_include_file
    # should be in the format "2022-01-01 00:00:00 UTC"
    String withdrawn_timestamp
  }

  meta {
    description: "Given a list of samples for a callset, withdraw samples from GVS that are not included by marking them as 'withdrawn' in the sample_info table with a passed in timestamp."
    # Might not be strictly necessary to make this volatile, but just in case:
    volatile: true
  }

  command <<<
    set -o errexit -o nounset -o xtrace -o pipefail
    echo "project_id = ~{project_id}" > ~/.bigqueryrc

    # get just the sample_name values from sample_names_to_include_file based on the
    # sample_name_column_name_in_file into sample_names.tsv
    col_num=$(head -n1 ~{sample_names_to_include_file} | tr '\t' '\n' | grep -Fxn '~{sample_name_column_name_in_file}' | cut -f1 -d:)
    awk "{print \$$col_num}" ~{sample_names_to_include_file} | sed 1d > sample_names.tsv

    # make sure that we end up with some samples to make the temp table, warn and exit if empty
    num_samples=$(cat sample_names.tsv | wc -l)
    if [ $num_samples -eq 0 ]; then
      echo "No sample names for callset produced. Exiting"
      exit 0
    fi

    # create the temp table (expires in 1 day)
    bq --apilog=false --project_id=~{project_id} mk --expiration=86400 ~{dataset_name}.current_callset_samples "sample_name:STRING"
    # populate the temp table
    bq --apilog=false load --project_id=~{project_id} -F "tab" ~{dataset_name}.current_callset_samples sample_names.tsv

    # join on the temp table to figure out which samples should be marked as withdrawn
    bq --apilog=false --project_id=~{project_id} query --format=csv --use_legacy_sql=false \
      "UPDATE \`~{dataset_name}.sample_info\` AS samples SET withdrawn = '~{withdrawn_timestamp}' \
        WHERE NOT EXISTS \
        (SELECT * \
        FROM \`~{project_id}.~{dataset_name}.current_callset_samples\` AS callset \
        WHERE \
        samples.sample_name = callset.sample_name) \
        AND NOT samples.is_control \
        AND withdrawn IS NULL" > log_message.txt

    cat log_message.txt | sed -e 's/Number of affected rows: //' > rows_updated.txt
  >>>
  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:4.2.5.0"
    memory: "3.75 GB"
    disks: "local-disk 10 HDD"
    cpu: 1
  }
  output {
    Int num_rows_updated = read_int("rows_updated.txt")
  }
}

