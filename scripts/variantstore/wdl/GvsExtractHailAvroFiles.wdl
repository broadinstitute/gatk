version 1.0

workflow GvsExtractHailAvroFiles {
    input {
        String project_id
        String dataset
        String write_prefix
        String filter_set_name
    }

    call OutputPaths { input: go = true }

    call ExtractAvroFiles {
        input:
            project_id = project_id,
            dataset = dataset,
            filter_set_name = filter_set_name,
            write_peer = OutputPaths.out
    }
}


task OutputPaths {
    input {
        Boolean go = true
    }
    command <<<
    >>>
    output {
        File out = stdout()
    }
    runtime {
        docker: "ubuntu:latest"
    }
}


task ExtractAvroFiles {
    input {
        String project_id
        String dataset
        String write_peer
        String filter_set_name
    }
    command <<<
        set -o errexit -o nounset -o xtrace -o pipefail
        echo "project_id = ~{project_id}" > ~/.bigqueryrc

        echo "write_peer is ~{write_peer}"

        write_prefix=$(dirname ~{write_peer})

        bq query --nouse_legacy_sql --project_id=~{project_id} "
            EXPORT DATA OPTIONS(
            uri='${write_prefix}/avro/sample_mapping/sample_mapping_*.avro', format='AVRO', compression='SNAPPY') AS
            SELECT sample_id, sample_name, "40",
            "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list" as intervals_file
            FROM \`~{project_id}.~{dataset}.sample_info\`
            ORDER BY sample_id
        "

        bq query --nouse_legacy_sql --project_id=~{project_id} "
            EXPORT DATA OPTIONS(
            uri='${write_prefix}/avro/vqsr_filtering_data_*.avro', format='AVRO', compression='SNAPPY') AS
            SELECT location, type as model, ref, alt, vqslod, yng_status
            FROM \`~{project_id}.~{dataset}.filter_set_info`
            WHERE filter_set_name = '~{filter_set_name}'
            ORDER BY location
        "

        bq query --nouse_legacy_sql --project_id=~{project_id} "
            EXPORT DATA OPTIONS(
            uri='${write_prefix}/avro/site_filtering_data/site_filtering_data_*.avro', format='AVRO', compression='SNAPPY') AS
            SELECT location, filters
            FROM \`~{project_id}.~{dataset}.filter_set_sites\`
            WHERE filter_set_name = '~{filter_set_name}'
            ORDER BY location
        "

        bq query --nouse_legacy_sql --project_id=~{project_id} "
            -- TODO handle superpartitioning, i.e. > 1 vet / ref_ranges table
            EXPORT DATA OPTIONS(
            uri='${write_prefix}/avro/vet/vet_001_*.avro', format='AVRO', compression='SNAPPY') AS
            SELECT location, sample_id, ref, REPLACE(alt,',<NON_REF>','') alt, call_GT as GT, call_AD as AD, call_GQ as GQ, cast(SPLIT(call_pl,',')[OFFSET(0)] as int64) as RGQ
            FROM \`~{project_id}.~{dataset}.vet_001\`
            ORDER BY location
        "

        bq query --nouse_legacy_sql --project_id=~{project_id} "
            -- TODO handle superpartitioning, i.e. > 1 vet / ref_ranges table
            EXPORT DATA OPTIONS(
            uri='${write_prefix}/avro/ref_ranges/ref_ranges_001_*.avro', format='AVRO', compression='SNAPPY') AS
            SELECT location, sample_id, length, state
            FROM \`~{project_id}.~{dataset}.ref_ranges_001\`
            ORDER BY location
        "
    >>>

    output {
        Boolean done = true
    }

    runtime {
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:latest"
        disk: "local-disk 200 HDD"
    }
}