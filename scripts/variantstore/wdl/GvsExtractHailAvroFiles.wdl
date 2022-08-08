version 1.0

workflow GvsExtractHailAvroFiles {
    input {
        String project_id
        String dataset
        String write_prefix
    }

    Array[String] array = ["first", "second", "third"]
    File f = write_lines(array)

    call ExtractAvroFiles {
        input:
            project_id = project_id,
            dataset = dataset,
            write_peer = f
    }
}


task ExtractAvroFiles {
    input {
        String project_id
        String dataset
        String write_peer
    }
    command <<<
        set -o errexit -o nounset -o xtrace -o pipefail
        echo "project_id = ~{project_id}" > ~/.bigqueryrc

        echo "write_peer is ~{write_peer}"

        write_prefix=$(dirname ~{write_peer})

        bq query --nouse_legacy_sql --project_id=~{project_id} '
            -- TODO handle superpartitioning, i.e. > 1 vet / ref_ranges table
            EXPORT DATA OPTIONS(uri="${write_prefix}/vet_001_*.avro", format="AVRO", compression="SNAPPY") AS
            SELECT location, sample_id, ref, REPLACE(alt,",<NON_REF>","") alt, call_GT as GT, call_AD as AD, call_GQ as GQ, cast(SPLIT(call_pl,",")[OFFSET(0)] as int64) as RGQ
            FROM `~{project_id}.~{dataset}.vet_001`
            ORDER BY location
        '
    >>>

    output {
        Boolean done = true
    }

    runtime {
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:latest"
        disk: "local-disk 200 HDD"
    }
}