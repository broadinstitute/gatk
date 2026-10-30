version 1.0

import "../wdl/GvsUtils.wdl" as GvsUtils

workflow SearchGVCFsAtSite {
    input {
        String project_id
        String dataset_name
        String locus
        Array[String] sample_names
        String sample_group_name
        String gvcf_lines_table_name
        String gvcf_info_fields_table_name
    }

    meta {
        description: "Reads relevant lines of GVCFs for the specified sample names at a specified locus"
    }

    call GvsUtils.GetToolVersions {}

    call QueryGVCFPaths {
        input:
            project_id = project_id,
            dataset_name = dataset_name,
            sample_names = sample_names,
            variants_docker = GetToolVersions.variants_docker,
    }

    call ReadGVCFs {
        input:
            paths_json = QueryGVCFPaths.paths_json,
            locus = locus,
            sample_group_name = sample_group_name,
            variants_docker = GetToolVersions.variants_docker,
    }

    # Intentionally unused: this workflow exists solely to run this task as a side effect; there are no outputs to consume.
    #@ except: UnusedCall
    call UploadGVCFContent {
        input:
            project_id = project_id,
            dataset_name = dataset_name,
            gvcf_lines_tsv = ReadGVCFs.gvcf_lines_tsv,
            gvcf_info_fields_tsv = ReadGVCFs.gvcf_info_fields_tsv,
            gvcf_lines_table_name = gvcf_lines_table_name,
            gvcf_info_fields_table_name = gvcf_info_fields_table_name,
            variants_docker = GetToolVersions.variants_docker,
    }

    output {
        File gvcf_content = ReadGVCFs.gvcf_content_json
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
        WHERE si.sample_name in ('"${formatted_sample_names}"')

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
        String sample_group_name
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

        # Generate a TSV file including embedded lines from the gVCFs (both raw and reblocked) in a format suitable for
        # loading into BigQuery.
        jq --raw-output '(.[0] | keys_unsorted) as $keys | $keys, map([.[ $keys[] ]|tostring])[] | @tsv' gvcf_content.json > gvcf_content_unordered.tsv

        cat gvcf_content_unordered.tsv | awk 'BEGIN {OFS="\t";} {print $3, $4, $1, $2, $5, $6, "~{sample_group_name}"}' > gvcf_lines.tsv

        # Generate a TSV file that splits out the info fields in a format suitable for loading in to BigQuery, e.g.
        # sample_id	sample_name	sample_group   vcf	format	info
        # 444444	5555555	foxtrot_gq40	raw	GT	0/0
        # 444444	5555555	foxtrot_gq40	raw	AD	41,0
        # 444444	5555555	foxtrot_gq40	raw	DP	41
        # 444444	5555555	foxtrot_gq40	raw	GQ	83
        # 444444	5555555	foxtrot_gq40	raw	MIN_DP	35
        # 444444	5555555	foxtrot_gq40	raw	PL	0,83,1485
        # 444444	5555555	foxtrot_gq40	raw	SPL	0,83,255
        # 444444	5555555	foxtrot_gq40	raw	ICNT	40,0
        # 444444	5555555	foxtrot_gq40	reblocked	GT	0/0
        # 444444	5555555	foxtrot_gq40	reblocked	DP	33
        # 444444	5555555	foxtrot_gq40	reblocked	GQ	40
        python3 /app/dst_2716_split_vcf_info_fields.py gvcf_lines.tsv > gvcf_info_fields.tsv
    >>>
    runtime {
        docker: variants_docker
    }
    output {
        File gvcf_content_json = "gvcf_content.json"
        File gvcf_lines_tsv = "gvcf_lines.tsv"
        File gvcf_info_fields_tsv = "gvcf_info_fields.tsv"
    }
}

task UploadGVCFContent {
    input {
        String project_id
        String dataset_name
        File gvcf_lines_tsv
        File gvcf_info_fields_tsv
        String gvcf_lines_table_name
        String gvcf_info_fields_table_name
        String variants_docker
    }

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        bq --apilog=false load --source_format=CSV --field_delimiter="\t" --skip_leading_rows=1 \
            --project_id=~{project_id} \
            --schema "sample_id:INTEGER,sample_name:STRING,gvcf_path:STRING,reblocked_gvcf:STRING,gvcf_line:STRING,reblocked_gvcf_line:STRING,sample_group_name:STRING" \
            ~{dataset_name}.~{gvcf_lines_table_name}\
            ~{gvcf_lines_tsv}

        bq --apilog=false load --source_format=CSV --field_delimiter="\t" --skip_leading_rows=1 \
            --project_id=~{project_id} \
            --schema "sample_id:INTEGER,sample_name:STRING,sample_group_name:STRING,vcf:STRING,format:STRING,info:STRING" \
            ~{dataset_name}.~{gvcf_info_fields_table_name}\
            ~{gvcf_info_fields_tsv}
    >>>
    runtime {
        docker: variants_docker
    }
    output {
        Boolean done = true
    }
}