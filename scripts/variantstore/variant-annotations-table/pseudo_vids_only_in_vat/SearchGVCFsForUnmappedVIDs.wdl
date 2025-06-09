version 1.0

import "../../wdl/GvsUtils.wdl" as GvsUtils

workflow SearchGVCFsForUnmappedVIDs {
    input {
        String project_id
        String dataset_name
    }

    call GvsUtils.GetToolVersions {}


    call SearchGVCFsForUnmappedVIDsTask {
        input:
            project_id = project_id,
            dataset_name = dataset_name,
            variants_docker = GetToolVersions.variants_docker,
    }

    output {
        File unmapped_vids_gvcf_json = SearchGVCFsForUnmappedVIDsTask.unmapped_vids_gvcf_json
    }
}

task SearchGVCFsForUnmappedVIDsTask {
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
            --format=prettyjson '

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
            `echo.pseudo_vid_sample_id` psi
        JOIN
            `echo.sample_info` si
        ON
            psi.sample_id = si.sample_id
        JOIN
            `echo.sample_data_table` dt
        ON
            si.sample_name = dt.research_id
        JOIN
            `echo.pseudo_vid_mapping` map
        ON
            psi.location = map.input_location
            AND psi.ref = map.input_ref
            AND psi.alt = map.input_alt

        ' > unmapped_vid_gvcfs.json

        python3 /app/process_gvcf_variants.py unmapped_vid_gvcfs.json
    >>>
    runtime {
        docker: variants_docker
    }
    output {
        File unmapped_vids_gvcf_json = "unmapped_vid_gvcfs.json"
    }
}