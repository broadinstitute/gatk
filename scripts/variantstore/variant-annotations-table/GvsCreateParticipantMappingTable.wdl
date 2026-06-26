version 1.0

import "../wdl/GvsUtils.wdl" as Utils
import "../structs/Range.wdl" as Range

workflow GvsCreateParticipantMappingTable {
    input {
        String project_id
        String dataset
        String vat_table_name
        String participant_mapping_table_name
        Range? range_filter
    }

    call Utils.GetToolVersions

    # Intentionally unused: this workflow exists solely to run this task as a side effect; there are no outputs to consume.
    #@ except: UnusedCall
    call CreateParticipantMappingTable {
        input:
        variants_docker = GetToolVersions.variants_docker,
        project_id = project_id,
        dataset = dataset,
        vat_table_name = vat_table_name,
        participant_mapping_table_name = participant_mapping_table_name,
        range_filter = range_filter,
    }
}


task CreateParticipantMappingTable {
    input {
        String variants_docker
        String project_id
        String dataset
        String vat_table_name
        String participant_mapping_table_name
        Range? range_filter
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

    CREATE TABLE `~{project_id}.~{dataset}.~{participant_mapping_table_name}` AS
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
    ~{if (defined(range_filter)) then "WHERE aa.location >= ~{select_first([range_filter]).startLocation} AND aa.location < ~{select_first([range_filter]).endLocation}" else ""}
    GROUP BY vat.vid
    ORDER BY
        vidToLocation(vat.vid),
        SPLIT(vat.vid, "-")[OFFSET(2)],
        SPLIT(vat.vid, "-")[OFFSET(3)]
        '

    bq update --project_id=~{project_id} --clustering_fields=vid ~{dataset}.~{participant_mapping_table_name}

    >>>
    runtime {
        docker: variants_docker
    }
    output {
        Boolean done = true
    }
}
