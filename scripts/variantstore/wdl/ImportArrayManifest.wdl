version 1.0

workflow ImportArrayManifest {

  input {
    File extended_manifest_csv
    File manifest_schema_json
    String project_id
    String dataset_name
    String? table_name
 
    Int? preemptible_tries
    String? docker
  }

  String docker_final = select_first([docker, "us.gcr.io/broad-gatk/gatk:4.1.7.0"])
 
  call CreateManifestCsv {
    input:
      extended_manifest_csv = extended_manifest_csv,
      preemptible_tries = preemptible_tries,
      docker = docker_final
  }

  call LoadManifest {
    input:
      project_id = project_id,
      dataset_name = dataset_name,
      table_name = table_name,
      manifest_schema_json = manifest_schema_json,
      manifest_csv = CreateManifestCsv.manifest_csv,
      preemptible_tries = preemptible_tries,
      docker = docker_final
  }
  output {
    File manifest_csv = CreateManifestCsv.manifest_csv
    File manifest_ingest_csv = CreateManifestCsv.manifest_ingest_csv
    File manifest_sub_csv = CreateManifestCsv.manifest_sub_csv
    File manifest_proc_csv = CreateManifestCsv.manifest_proc_csv
  }
}

task LoadManifest {
  input {
    String project_id
    String dataset_name
    String? table_name
    File manifest_csv
    File manifest_schema_json
    # runtime
    Int? preemptible_tries
    String docker
    # String to add command for testing only. Can be ignored otherwise.
    String? for_testing_only
  }

  String ingest_table = dataset_name + "." + select_first([table_name, "probe_info"])

  parameter_meta {
    manifest_schema_json: {
      localization_optional: false
    }
  }
   
    command <<<
      set +e
      ~{for_testing_only}
      bq ls --project_id ~{project_id} ~{dataset_name} > /dev/null
      if [ $? -ne 0 ]; then
        echo "making dataset ~{project_id}.~{dataset_name}"
        bq mk --project_id=~{project_id} ~{dataset_name}
      fi
      bq show --project_id ~{project_id} ~{ingest_table} > /dev/null
      if [ $? -ne 0 ]; then
        echo "making table ~{ingest_table}"
        # create a site info table and load - schema and TSV header need to be the same order
        bq --location=US mk --project_id=~{project_id} ~{ingest_table} ~{manifest_schema_json}
      fi
      set -e

      bq load --location=US --project_id=~{project_id} --null_marker "null" --source_format=CSV ~{ingest_table} ~{manifest_csv} ~{manifest_schema_json}
    >>>
    runtime {
      docker: docker
      memory: "4 GB"
      disks: "local-disk " + 20 + " HDD"
      preemptible: select_first([preemptible_tries, 5])
      cpu: 2
  }

}

task CreateManifestCsv {
  input {
    File extended_manifest_csv

    # runtime
    Int? preemptible_tries
    String docker
  }

  Int disk_size = ceil(size(extended_manifest_csv, "GB") * 2.5) + 20

  meta {
    description: "Creates a tsv file for imort into BigQuery"
  }
  parameter_meta {
    extended_manifest_csv: {
      localization_optional: false
    }
  }
  command <<<
    set -e

    TMP_SORTED="manifest_ingest_sorted.csv"
    TMP_SUB="manifest_ingest_sub.csv"
    TMP_PROC="manifest_ingest_processed.csv"
    TMP="manifest_ingest.csv"

    # put integers in front of the chromosomes that are not numbered so that they end up ordered by X, Y and MT
    sed 's/,X,/,23X,/g; s/,Y,/,24Y,/g; s/,MT,/,25MT,/g' ~{extended_manifest_csv} > $TMP_SUB

    # sort the probes by chrom, position and then name so there is a specific ordering when we assign integers
    sort -t , -k23n,23 -k24n,24 -k2,2 $TMP_SUB > $TMP_SORTED

    # checking for != "build37Flag" skips the header row (we don't want that numbered)
    # only process rows with 29 fields - this skips some header info fields
    # also skip entries that are flagged, not matched or have index conflict
    awk -F ',' 'NF==29 && ($29!="ILLUMINA_FLAGGED" && $29!="INDEL_NOT_MATCHED" && $29!="INDEL_CONFLICT" && $29!="build37Flag") { flag=$29; if ($29=="PASS") flag="null"; print ++id","$2","$9","$23","$24","$25","$26","$27","flag }' $TMP_SORTED > $TMP_PROC

    # remove the integer prefixes for chromosomes X, Y and MT
    sed 's/,23X,/,X,/g; s/,24Y,/,Y,/g; s/,25MT,/,MT,/g' $TMP_PROC > $TMP

    echo "created file for ingest $TMP"
  >>>
  runtime {
      docker: docker
      memory: "4 GB"
      disks: "local-disk " + disk_size + " HDD"
      preemptible: select_first([preemptible_tries, 5])
      cpu: 2
  }
  output {
      File manifest_csv = "manifest_ingest.csv"
      File manifest_ingest_csv = "manifest_ingest_sorted.csv"
      File manifest_sub_csv = "manifest_ingest_sub.csv"
      File manifest_proc_csv = "manifest_ingest_processed.csv"
  
  }
}
