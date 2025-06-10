version 1.0

import "../../wdl/GvsUtils.wdl" as GvsUtils

workflow SearchGVCFsForUnmappedVIDs {
    input {
        String project_id
        String dataset_name
    }

    meta {
        description: "Reads relevant lines of GVCFs for samples with unmapped VIDs in the VAT, creates a TSV suitable for loading into BigQuery."
    }

    call GvsUtils.GetToolVersions {}

    call QueryGVCFPaths {
        input:
            project_id = project_id,
            dataset_name = dataset_name,
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

    call UploadGVCFContent {
        input:
            project_id = project_id,
            dataset_name = dataset_name,
            merged_tsv = MergeJSONs.merged_tsv,
            variants_docker = GetToolVersions.variants_docker,
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
    }

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        bq --apilog=false query --use_legacy_sql=false --max_rows=100000000 --project_id=~{project_id} \
            --format=json '

        SELECT
            si.sample_id,
            si.sample_name,
            map.chr,
            map.input_position,
            map.input_ref,
            map.input_alt,
            dt.reblocked_gvcf,
            dt.gvcf_path
        FROM
            `~{dataset_name}.pseudo_vid_sample_id` psi
        JOIN
            `~{dataset_name}.sample_info` si
        ON
            psi.sample_id = si.sample_id
        JOIN
            `~{dataset_name}.sample_data_table` dt
        ON
            si.sample_name = dt.research_id
        JOIN
            `~{dataset_name}.pseudo_vid_mapping` map
        ON
            psi.location = map.input_location
            AND psi.ref = map.input_ref
            AND psi.alt = map.input_alt

        ' | jq '.[]' | jq --compact-output . > paths.json

        # The above query returns a JSON array which we reformat into a file with one JSON object per line so it can be
        # split and scattered over downstream.

        split -l 300 paths.json path_shard_.
    >>>
    runtime {
        docker: variants_docker
    }
    output {
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
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        # Wrap the JSON objects in an array to make a valid JSON file.
        jq --slurp '.' ~{paths_json} > array.json

        python3 /app/process_gvcf_variants.py array.json > gvcf_content.json
    >>>
    runtime {
        docker: variants_docker
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

        python /app/reorder_gvcf_content_cols.py ~{merged_tsv} > reordered.tsv

        bq --apilog=false load --source_format=CSV --field_delimiter="\t" --skip_leading_rows=1 \
            --project_id=~{project_id} \
            --schema "sample_name:STRING,sample_id:INTEGER,chr:STRING,input_position:INTEGER,input_ref:STRING,input_alt:STRING,gvcf_path:STRING,reblocked_gvcf:STRING,gvcf_line:STRING,reblocked_gvcf_line:STRING" \
            ~{dataset_name}.pseudo_vid_gvcf_content \
            reordered.tsv
    >>>
    runtime {
        docker: variants_docker
    }
    output {
        Boolean done = true
    }
}