version 1.0

import "../wdl/GvsUtils.wdl" as Utils

workflow GvsCreateParticipantMappingTable {
    input {
        String project_id
        String dataset
        String vat_table_name
        String mapping_table_name
        Boolean chr20_only = false
    }

    call Utils.GetToolVersions

    call CreateParticipantMappingTable {
        input:
        variants_docker = GetToolVersions.variants_docker,
        project_id = project_id,
        dataset = dataset,
        vat_table_name = vat_table_name,
        mapping_table_name = mapping_table_name,
        chr20_only = chr20_only,
    }
}


task CreateParticipantMappingTable {
    input {
        String variants_docker
        String project_id
        String dataset
        String vat_table_name
        String mapping_table_name
        Boolean chr20_only
    }
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        # bq query --max_rows check: ok, results going into new participant mapping table

        bq --apilog=false query --nouse_legacy_sql --project_id=~{project_id} '

   CREATE TEMP FUNCTION vidToLocation(vid string)
    RETURNS int64
    AS (
        (CASE SPLIT(vid, "-")[OFFSET(0)]
                            WHEN "X" THEN 23
                            WHEN "Y" THEN 24
                            ELSE CAST(SPLIT(vid, "-")[OFFSET(0)] AS int64) END) * 1000000000000 +
                    CAST(SPLIT(vid, "-")[OFFSET(1)] AS int64)
    );

    CREATE TABLE `~{project_id}.~{dataset}.~{mapping_table_name}` AS
    SELECT vat.vid as vid, ARRAY_AGG(SAFE_CAST(si.sample_name as INT64) IGNORE NULLS) AS person_ids
        FROM `~{project_id}.~{dataset}.alt_allele` AS aa
                JOIN `~{project_id}.~{dataset}.sample_info` AS si
                    ON aa.sample_id = si.sample_id
                JOIN
            (SELECT vid,
                vidToLocation(vid) AS location,
                SPLIT(vid, "-")[OFFSET(2)] AS ref_allele,
                SPLIT(vid, "-")[OFFSET(3)] AS alt_allele
            FROM `~{project_id}.~{dataset}.~{vat_table_name}`
            GROUP BY vid, location) AS vat
        ON
            vat.ref_allele = aa.ref AND
            vat.alt_allele = aa.allele AND
            vat.location = aa.location
    ~{if (chr20_only) then "WHERE aa.location >= 20 * (1000 * 1000 * 1000 * 1000) AND aa.location < 21 * (1000 * 1000 * 1000 * 1000)" else ""}
    GROUP BY vat.vid
    ORDER BY
        vidToLocation(vat.vid),
        SPLIT(vat.vid, "-")[OFFSET(2)],
        SPLIT(vat.vid, "-")[OFFSET(3)]
        '

    bq update --project_id=~{project_id} --clustering_fields=vid ~{dataset}.~{mapping_table_name}

    >>>
    runtime {
        docker: variants_docker
    }
    output {
        Boolean done = true
    }
}
