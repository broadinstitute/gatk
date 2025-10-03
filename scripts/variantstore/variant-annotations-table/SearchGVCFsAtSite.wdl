version 1.0

import "../wdl/GvsUtils.wdl" as GvsUtils

workflow SearchGVCFsAtSite {
    input {
        String project_id
        String dataset_name
        String locus
        Array[String] echo_sample_names
        Array[String] foxtrot_sample_names
    }

    meta {
        description: "Reads relevant lines of GVCFs for samples with unmapped VIDs in the VAT and loads into BigQuery."
    }

    call GvsUtils.GetToolVersions {}

    call QueryGVCFPaths as QueryEcho {
        input:
            project_id = project_id,
            dataset_name = dataset_name,
            sample_names = echo_sample_names,
            variants_docker = GetToolVersions.variants_docker,
    }

    call QueryGVCFPaths as QueryFoxtrot {
        input:
            project_id = project_id,
            dataset_name = dataset_name,
            sample_names = foxtrot_sample_names,
            variants_docker = GetToolVersions.variants_docker,
    }

    call ReadGVCFs as ReadEchoGVCFs {
        input:
            paths_json = QueryEcho.paths_json,
            locus = locus,
            variants_docker = GetToolVersions.variants_docker,
    }

    call ReadGVCFs as ReadFoxtrotGVCFs {
        input:
            paths_json = QueryFoxtrot.paths_json,
            locus = locus,
            variants_docker = GetToolVersions.variants_docker,
    }

    output {
        File echo_content = ReadEchoGVCFs.gvcf_content_json
        File foxtrot_content = ReadFoxtrotGVCFs.gvcf_content_json
    }
}

task QueryGVCFPaths {
    input {
        String project_id
        String dataset_name
        Array[String] sample_names
        String variants_docker
    }

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        formatted_sample_names=$(echo ~{sep=' ' sample_names} | tr ' ' "\n" | sed 's/.*/"&"/' | paste -sd ',' -)

        bq --apilog=false query --use_legacy_sql=false --max_rows=100000000 --project_id=~{project_id} \
            --format=json '

        SELECT
            si.sample_id,
            si.sample_name,
            dt.reblocked_gvcf,
            dt.gvcf_path
        FROM
            `~{dataset_name}.sample_info` si
        JOIN
            `~{dataset_name}.reblocking_data_table` dt
        ON
            si.sample_name = dt.research_id
        WHERE si.sample_id in ('"${formatted_sample_names}"')

        ' > paths.json
    >>>
    runtime {
        docker: variants_docker
    }
    output {
        File paths_json = "paths.json"
    }
}

task ReadGVCFs {
    input {
        File paths_json
        String locus
        String variants_docker
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

        python3 /app/dst_2716_vat_investigation_read_site.py ~{paths_json} --locus ~{locus} > gvcf_content.json
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