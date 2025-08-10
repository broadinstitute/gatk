version 1.0

import "../../wdl/GvsUtils.wdl" as GvsUtils

workflow InspectDroppedDuplicateVCFs {
    input {
        String project_id
        String dataset_name
        Int query_mem_gib = 14
        Boolean retry = false
    }

    meta {
        description: "Reads relevant lines of GVCFs for samples with duplicates dropped from the VAT and loads into BigQuery."
    }

    call GvsUtils.GetToolVersions {}

    call QueryGVCFPaths {
        input:
            project_id = project_id,
            dataset_name = dataset_name,
            mem_gib = query_mem_gib,
            retry = retry,
            variants_docker = GetToolVersions.variants_docker,
    }

    scatter (path_shard in QueryGVCFPaths.path_shards) {
        call ReadGVCFs {
            input:
                paths_json = path_shard,
                variants_docker = GetToolVersions.variants_docker,
        }
    }

    call GvsUtils.MergeJSONs {
        input:
            input_files = ReadGVCFs.gvcf_content_json,
            variants_docker = GetToolVersions.variants_docker,
    }

    if (!retry) {
        call UploadGVCFContent {
            input:
                project_id = project_id,
                dataset_name = dataset_name,
                merged_tsv = MergeJSONs.merged_tsv,
                variants_docker = GetToolVersions.variants_docker,
        }
    }

    output {
        File merged_json = MergeJSONs.merged_json
        File? merged_tsv = MergeJSONs.merged_tsv
    }
}

task QueryGVCFPaths {
    input {
        String project_id
        String dataset_name
        String variants_docker
        Boolean retry
        Int mem_gib
        Int split_batch_size = 300
        Int split_suffix_size = 3
    }

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        if [[ "~{retry}" = "true" ]]
        then
            bq --apilog=false query --use_legacy_sql=false --max_rows=100000000 --project_id=~{project_id} \
                --format=json '

            SELECT
                *
            FROM
                `~{dataset_name}.dropped_duplicates_gvcf_content`
            ORDER BY sample_id, location, ref, alt

            ' > raw_paths.json
        else
            bq --apilog=false query --use_legacy_sql=false --max_rows=100000000 --project_id=~{project_id} \
                --format=json '

            SELECT
                si.sample_id,
                si.sample_name,
                ddaa.location,
                ddaa.ref,
                ddaa.allele,
                dt.reblocked_gvcf,
                dt.gvcf_path
            FROM
                `~{dataset_name}.dropped_duplicates_alt_allele` ddaa
            JOIN
                `~{dataset_name}.sample_info` si
            ON
                ddaa.sample_id = si.sample_id
            JOIN
                `~{dataset_name}.sample_data_table` dt
            ON
                si.sample_name = dt.research_id
            ORDER BY si.sample_id, ddaa.location, ddaa.ref, ddaa.allele

            ' > raw_paths.json
        fi

        jq '.[]' raw_paths.json | jq --compact-output . > paths.json

        # The above queries return a JSON array which we reformat into a file with one JSON object per line so it can be
        # split and scattered over downstream.

        # ~{split_batch_size} lines per file, ~{split_suffix_size} letter suffix
        split -l ~{split_batch_size} -a ~{split_suffix_size} paths.json path_shard_.
    >>>
    runtime {
        docker: variants_docker
        memory: mem_gib + " GiB"
        preemptible: 2
    }
    output {
        File raw_paths_json = "raw_paths.json"
        File paths_json = "paths.json"
        Array[File] path_shards = glob("path_shard_.*")
    }
}

task ReadGVCFs {
    input {
        String variants_docker
        File paths_json
    }

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        # Specify a GCS auth token for htslib/bcftools. Note this token only works for a limited time, but that's fine
        # for this use case. If a persistent solution is required, see
        # https://broadinstitute.slack.com/archives/C0544AAC70D/p1696360070640409

        # Let's not log the token.
        set +o xtrace
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        set -o xtrace

        # Wrap the JSON objects in an array to make a valid JSON file.
        jq --slurp '.' ~{paths_json} > array.json

        python3 /app/process_gvcfs_for_dropped_duplicates.py array.json > gvcf_content.json
    >>>
    runtime {
        docker: variants_docker
        preemptible: 2
    }
    output {
        File gvcf_content_json = "gvcf_content.json"
    }
}

task UploadGVCFContent {
    input {
        String project_id
        String dataset_name
        File merged_tsv
        String variants_docker
    }

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        python /app/reorder_dropped_duplicate_columns.py ~{merged_tsv} > reordered.tsv

        bq --apilog=false load --source_format=CSV --field_delimiter="\t" --skip_leading_rows=1 \
            --project_id=~{project_id} \
            --schema "sample_id:INTEGER,sample_name:STRING,location:INTEGER,ref:STRING,alt:STRING,chr:STRING,position:INTEGER,gvcf_path:STRING,reblocked_gvcf:STRING,gvcf_line:STRING,reblocked_gvcf_line:STRING" \
            ~{dataset_name}.dropped_duplicates_gvcf_content \
            reordered.tsv
    >>>
    runtime {
        docker: variants_docker
    }
    output {
        Boolean done = true
    }
}
