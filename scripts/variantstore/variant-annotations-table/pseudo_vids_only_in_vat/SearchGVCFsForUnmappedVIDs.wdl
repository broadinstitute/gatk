version 1.0

import "../../wdl/GvsUtils.wdl" as GvsUtils

workflow SearchGVCFsForUnmappedVIDs {
    input {
        String project_id
        String dataset_name
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

        -- DEBUG
        LIMIT 200

        ' | jq '.[]' | jq --compact-output . > paths.json

        split -l 100 paths.json path_shard_
    >>>
    runtime {
        docker: variants_docker
    }
    output {
        File paths_json = "paths.json"
        Array[File] path_shards = glob("path_shard_*")
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

        # For htslib in bcftools to access GCS. Note this token will only work for a limited time. If a persistent
        # solution is required, see https://broadinstitute.slack.com/archives/C0544AAC70D/p1696360070640409
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
