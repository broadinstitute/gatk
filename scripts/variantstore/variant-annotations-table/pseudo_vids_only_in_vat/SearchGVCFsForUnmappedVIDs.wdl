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

    scatter (json_shard in QueryGVCFPaths.unmapped_vids_gvcf_shards) {
        call ReadGVCFs {
            input:
                unmapped_vid_gcf_paths_json = json_shard,
                variants_docker = GetToolVersions.variants_docker,
        }
    }

    call GvsUtils.MergeJSONs {
        input:
            input_files = ReadGVCFs.unmapped_vids_gvcf_json,
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

        ' | jq '.[]' | jq --compact-output . > unmapped_vid_gvcf_paths.json

        split -l 100 unmapped_vid_gvcf_paths.json unmapped_vid_gvcf_shard_
    >>>
    runtime {
        docker: variants_docker
    }
    output {
        File unmapped_vids_gvcf_paths_json = "unmapped_vid_gvcf_paths.json"
        Array[File] unmapped_vids_gvcf_shards = glob("unmapped_vid_gvcf_shard_*")
    }
}

task ReadGVCFs {
    input {
        String variants_docker
        File unmapped_vid_gcf_paths_json
    }

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        # Required for htslib in bcftools to access GCS.
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        python3 /app/process_gvcf_variants.py ~{unmapped_vid_gcf_paths_json} > unmapped_vid_gvcfs.json
    >>>
    runtime {
        docker: variants_docker
    }
    output {
        File unmapped_vids_gvcf_json = "unmapped_vid_gvcfs.json"
    }
}
