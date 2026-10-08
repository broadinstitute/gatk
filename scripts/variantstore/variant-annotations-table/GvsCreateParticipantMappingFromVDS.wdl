version 1.0

# Participant mapping built from the VDS (VS-2013). hail_create_participant_mapping_inputs.py writes the carrier pairs
# and the VDS keys that are not left-aligned; this workflow gives those keys their VAT VIDs and builds, in the given
# dataset, `<participant_mapping_table_name>_base` (het_ids, hom_ids, hemi_ids per VID), its provenance table
# `<participant_mapping_table_name>_provenance`, and the view `<participant_mapping_table_name>` with the delivered
# (vid, person_ids) shape.
import "../wdl/GvsUtils.wdl" as Utils

workflow GvsCreateParticipantMappingFromVDS {
    input {
        String project_id
        String dataset_name
        String fq_vat_table
        String participant_mapping_table_name
        # The Hail script's --output prefix.
        String mapping_inputs_path
        String reference_name = "hg38"

        String? git_branch_or_tag
        String? basic_docker
        String? variants_docker
    }

    parameter_meta {
        fq_vat_table: {
            help: "The VAT, as project.dataset.table. Every VID the mapping files a carrier under must appear in it."
        }
        mapping_inputs_path: {
            help: "GCS prefix passed as --output to hail_create_participant_mapping_inputs.py, holding keys_to_normalize.vcf.bgz and pairs.parquet."
        }
        participant_mapping_table_name: {
            help: "Name of the delivered view. The base and provenance tables, and the staging tables, take it as a prefix. None may already exist, apart from the staging tables, which are replaced."
        }
    }

    String mapping_inputs = sub(mapping_inputs_path, "/$", "")

    if (!defined(basic_docker) || !defined(variants_docker)) {
        call Utils.GetToolVersions {
            input:
                git_branch_or_tag = git_branch_or_tag,
        }
    }

    String effective_basic_docker = select_first([basic_docker, GetToolVersions.basic_docker])
    String effective_variants_docker = select_first([variants_docker, GetToolVersions.variants_docker])

    call Utils.GetReference {
        input:
            reference_name = reference_name,
            basic_docker = effective_basic_docker,
    }

    call NormalizeKeys {
        input:
            keys_to_normalize_vcf = "~{mapping_inputs}/keys_to_normalize.vcf.bgz",
            ref = GetReference.reference.reference_fasta,
            variants_docker = effective_variants_docker,
    }

    call LoadMappingInputs {
        input:
            project_id = project_id,
            dataset_name = dataset_name,
            participant_mapping_table_name = participant_mapping_table_name,
            pairs_parquet_path = "~{mapping_inputs}/pairs.parquet",
            key_to_vid = NormalizeKeys.key_to_vid,
            variants_docker = effective_variants_docker,
    }

    call CreateMappingTables {
        input:
            project_id = project_id,
            dataset_name = dataset_name,
            fq_vat_table = fq_vat_table,
            participant_mapping_table_name = participant_mapping_table_name,
            pairs_table = LoadMappingInputs.pairs_table,
            key_to_vid_table = LoadMappingInputs.key_to_vid_table,
            variants_docker = effective_variants_docker,
    }

    call CheckMappingTables {
        input:
            project_id = project_id,
            fq_vat_table = fq_vat_table,
            base_table = CreateMappingTables.base_table,
            provenance_table = CreateMappingTables.provenance_table,
            variants_docker = effective_variants_docker,
    }

    output {
        File key_to_vid = NormalizeKeys.key_to_vid
        File acceptance_checks = CheckMappingTables.acceptance_checks
        String base_table = CreateMappingTables.base_table
        String provenance_table = CreateMappingTables.provenance_table
        String view = CreateMappingTables.view
    }
}

task NormalizeKeys {
    input {
        File keys_to_normalize_vcf
        File ref
        String variants_docker
    }

    Int disk_size = ceil(size(ref, "GB") * 2 + size(keys_to_normalize_vcf, "GB") * 5) + 20

    # Writes key_to_vid.tsv, one headerless `vds_key<TAB>vid` row per input key.
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        # The same normalization RemoveDuplicatesFromSitesOnlyVCF in GvsCreateVATfromVDS.wdl applies, against the same
        # reference, so each key gets the VID the VAT gave its allele. ID passes through untouched.
        bcftools norm -m- --check-ref w -f ~{ref} ~{keys_to_normalize_vcf} -O v -o normalized.vcf

        # VID format: contig without the chr prefix.
        bcftools query -f '%ID\t%CHROM-%POS-%REF-%ALT\n' normalized.vcf \
            | awk -F'\t' -v OFS='\t' '{ sub(/^chr/, "", $2); print }' > key_to_vid.tsv

        # Every key is biallelic, so norm -m- neither splits nor drops; anything else means the input is not what the
        # Hail script writes.
        input_keys=$(bcftools view -H ~{keys_to_normalize_vcf} | wc -l)
        output_rows=$(wc -l < key_to_vid.tsv)
        distinct_keys=$(cut -f1 key_to_vid.tsv | sort -u | wc -l)
        if [[ ${input_keys} -ne ${output_rows} || ${input_keys} -ne ${distinct_keys} ]]
        then
            echo "Expected one row per key: ${input_keys} keys in, ${output_rows} rows and ${distinct_keys} distinct keys out." >&2
            exit 1
        fi

        # The Hail script flags a key only when its indel can move left, so every one must come out at a new VID.
        unmoved=$(awk -F'\t' '$1 == $2' key_to_vid.tsv | wc -l)
        if [[ ${unmoved} -ne 0 ]]
        then
            echo "${unmoved} flagged keys normalized to themselves, e.g.:" >&2
            awk -F'\t' '$1 == $2' key_to_vid.tsv | head -5 >&2
            exit 1
        fi
    >>>

    runtime {
        docker: variants_docker
        maxRetries: 3
        memory: "4 GB"
        preemptible: 3
        cpu: 2
        disks: "local-disk " + disk_size + " HDD"
    }

    output {
        File key_to_vid = "key_to_vid.tsv"
    }
}

task LoadMappingInputs {
    input {
        String project_id
        String dataset_name
        String participant_mapping_table_name
        String pairs_parquet_path
        File key_to_vid
        String variants_docker
    }

    String pairs_table_name = "~{participant_mapping_table_name}_pairs"
    String key_to_vid_table_name = "~{participant_mapping_table_name}_key_to_vid"

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        # Staging tables, replaced on a rerun. Load jobs are free; the pairs are on the order of 2.5e12 rows genome-wide.
        # The glob skips Spark's _SUCCESS marker.
        bq --apilog=false load --project_id=~{project_id} --replace --source_format=PARQUET \
            ~{dataset_name}.~{pairs_table_name} '~{pairs_parquet_path}/*.parquet'

        bq --apilog=false load --project_id=~{project_id} --replace --source_format=CSV --field_delimiter='\t' \
            ~{dataset_name}.~{key_to_vid_table_name} ~{key_to_vid} vds_key:STRING,vid:STRING
    >>>

    runtime {
        docker: variants_docker
    }

    output {
        String pairs_table = pairs_table_name
        String key_to_vid_table = key_to_vid_table_name
    }
}

task CreateMappingTables {
    input {
        String project_id
        String dataset_name
        String fq_vat_table
        String participant_mapping_table_name
        String pairs_table
        String key_to_vid_table
        String variants_docker
    }

    String base_table_name = "~{participant_mapping_table_name}_base"
    String provenance_table_name = "~{participant_mapping_table_name}_provenance"
    String fq_dataset = "~{project_id}.~{dataset_name}"

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        cat > create_mapping_tables.sql <<'SQL'
        DECLARE non_vat_vids INT64;
        DECLARE non_vat_examples ARRAY<STRING>;

        CREATE TEMP FUNCTION vidToLocation(vid STRING)
        RETURNS INT64
        AS (
            (CASE SPLIT(vid, '-')[OFFSET(0)]
                WHEN 'X' THEN 23
                WHEN 'Y' THEN 24
                ELSE CAST(SPLIT(vid, '-')[OFFSET(0)] AS INT64) END) * 1000000000000 +
            CAST(SPLIT(vid, '-')[OFFSET(1)] AS INT64)
        );

        -- Carrier counts per VDS key, from one scan of the pairs. Everything except the base table reads this rather
        -- than the pairs. A key is the VID unless bcftools norm gave it another.
        CREATE TEMP TABLE key_counts AS
        SELECT p.vds_key,
               COALESCE(k.vid, p.vds_key) AS vid,
               p.n_het, p.n_hom, p.n_hemi
        FROM (
            SELECT vds_key,
                   COUNTIF(zygosity = 'het') AS n_het,
                   COUNTIF(zygosity = 'hom') AS n_hom,
                   COUNTIF(zygosity = 'hemi') AS n_hemi
            FROM `~{fq_dataset}.~{pairs_table}`
            GROUP BY vds_key
        ) AS p
        LEFT JOIN `~{fq_dataset}.~{key_to_vid_table}` AS k USING (vds_key);

        -- Every VID a carrier is filed under must be a VAT VID. One that is not means the VDS and the VAT disagree about
        -- which variants exist, and the mapping would hand out an ID nobody can look up.
        SET (non_vat_vids, non_vat_examples) = (
            SELECT AS STRUCT COUNT(*), ARRAY_AGG(c.vid ORDER BY c.vid LIMIT 10)
            FROM (SELECT DISTINCT vid FROM key_counts) AS c
            LEFT JOIN (SELECT DISTINCT vid FROM `~{fq_vat_table}`) AS vat USING (vid)
            WHERE vat.vid IS NULL
        );
        IF non_vat_vids > 0 THEN
            RAISE USING MESSAGE = FORMAT('%d mapped VIDs are not VAT VIDs, e.g. %s',
                                         non_vat_vids, ARRAY_TO_STRING(non_vat_examples, ', '));
        END IF;

        -- Clustered at creation: a bq update after the fact does not take. No ORDER BY; see VS-2013.
        CREATE TABLE `~{fq_dataset}.~{base_table_name}`
        CLUSTER BY vid
        AS
        WITH person_zygosity AS (
            -- One row per (vid, person). A person appears once per key, so only VIDs that more than one key normalizes
            -- to can repeat a person here. Such a person keeps the strongest zygosity any of their keys gives them,
            -- so no person lands in two of the three arrays.
            SELECT COALESCE(k.vid, p.vds_key) AS vid,
                   p.person_id,
                   MAX(CASE p.zygosity WHEN 'het' THEN 1 WHEN 'hemi' THEN 2 WHEN 'hom' THEN 3 END) AS zygosity_rank
            FROM `~{fq_dataset}.~{pairs_table}` AS p
            LEFT JOIN `~{fq_dataset}.~{key_to_vid_table}` AS k USING (vds_key)
            -- By position: an unqualified vid here could resolve to k.vid, which is NULL for nearly every key.
            GROUP BY 1, 2
        )
        SELECT vid,
               -- ARRAY_AGG over nothing but NULLs is NULL, and ARRAY_CONCAT with a NULL is NULL.
               IFNULL(ARRAY_AGG(IF(zygosity_rank = 1, person_id, NULL) IGNORE NULLS), []) AS het_ids,
               IFNULL(ARRAY_AGG(IF(zygosity_rank = 3, person_id, NULL) IGNORE NULLS), []) AS hom_ids,
               IFNULL(ARRAY_AGG(IF(zygosity_rank = 2, person_id, NULL) IGNORE NULLS), []) AS hemi_ids
        FROM person_zygosity
        GROUP BY vid;

        -- One row per normalized key that has carriers. input_location, input_ref and input_alt join to alt_allele on
        -- (location, ref, allele), which explains a participant filed under a VID whose own location has no alt_allele
        -- row for them.
        CREATE TABLE `~{fq_dataset}.~{provenance_table_name}`
        CLUSTER BY vid
        AS
        SELECT vid,
               vidToLocation(vds_key) AS input_location,
               SPLIT(vds_key, '-')[OFFSET(2)] AS input_ref,
               SPLIT(vds_key, '-')[OFFSET(3)] AS input_alt,
               n_het, n_hom, n_hemi
        FROM key_counts
        WHERE vid != vds_key;

        CREATE VIEW `~{fq_dataset}.~{participant_mapping_table_name}` AS
        SELECT vid, ARRAY_CONCAT(het_ids, hom_ids, hemi_ids) AS person_ids
        FROM `~{fq_dataset}.~{base_table_name}`;
        SQL

        # bq query --max_rows check: ok, results go to the new tables
        bq --apilog=false query --nouse_legacy_sql --project_id=~{project_id} "$(cat create_mapping_tables.sql)"
    >>>

    runtime {
        docker: variants_docker
    }

    output {
        String base_table = "~{fq_dataset}.~{base_table_name}"
        String provenance_table = "~{fq_dataset}.~{provenance_table_name}"
        String view = "~{fq_dataset}.~{participant_mapping_table_name}"
    }
}

# The VS-2013 acceptance criteria, against the tables as built. Fails the workflow, leaving the tables in place for
# inspection, if any is not met. One scan of the base table's arrays; everything else reads vid columns.
task CheckMappingTables {
    input {
        String project_id
        String fq_vat_table
        String base_table
        String provenance_table
        String variants_docker
    }

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        cat > acceptance_checks.sql <<'SQL'
        -- Only the VAT VIDs on contigs the mapping covers, so a run restricted to whole contigs is checked as a whole.
        WITH contigs AS (
            SELECT DISTINCT SPLIT(vid, '-')[OFFSET(0)] AS contig FROM `~{base_table}`
        ),
        vat AS (
            SELECT vid, ANY_VALUE(gvs_all_sc) AS sc
            FROM `~{fq_vat_table}`
            WHERE SPLIT(vid, '-')[OFFSET(0)] IN (SELECT contig FROM contigs)
            GROUP BY vid
        ),
        mapping AS (
            SELECT vid,
                   ARRAY_LENGTH(het_ids) + ARRAY_LENGTH(hom_ids) + ARRAY_LENGTH(hemi_ids) AS n,
                   (SELECT COUNT(DISTINCT person_id)
                    FROM UNNEST(ARRAY_CONCAT(het_ids, hom_ids, hemi_ids)) AS person_id) AS n_distinct
            FROM `~{base_table}`
        ),
        provenance AS (
            SELECT DISTINCT vid FROM `~{provenance_table}`
        ),
        j AS (
            SELECT m.vid IS NOT NULL AS in_mapping,
                   v.vid IS NOT NULL AS in_vat,
                   p.vid IS NOT NULL AS in_provenance,
                   m.n, m.n_distinct, v.sc
            FROM mapping AS m
            FULL JOIN vat AS v ON v.vid = m.vid
            LEFT JOIN provenance AS p ON p.vid = COALESCE(m.vid, v.vid)
        )
        SELECT
            COUNTIF(in_mapping) AS mapped_vids,
            COUNTIF(in_provenance) AS provenance_vids,
            -- Criterion 1. Outside the provenance VIDs the mapping count is exactly gvs_all_sc; inside, where synonym
            -- carriers are unioned in, it may exceed it but never fall below.
            COUNTIF(in_mapping AND in_vat AND NOT in_provenance AND n = sc) AS exact_outside_provenance,
            COUNTIF(in_mapping AND in_vat AND NOT in_provenance AND n != sc) AS mismatch_outside_provenance,
            COUNTIF(in_provenance AND n > sc) AS above_inside_provenance,
            COUNTIF(in_provenance AND n < sc) AS below_inside_provenance,
            -- Criterion 2. Every mapped VID is a VAT VID, and every VAT VID with carriers has a mapping row.
            COUNTIF(in_mapping AND NOT in_vat) AS mapped_not_in_vat,
            COUNTIF(in_vat AND NOT in_mapping AND sc > 0) AS vat_carriers_without_row,
            -- Criterion 3. No person appears twice for one VID, within an array or across them.
            COUNTIF(n != n_distinct) AS vids_with_duplicate_person
        FROM j
        SQL

        # bq query --max_rows check: ok, a single row of counts
        bq --apilog=false query --nouse_legacy_sql --project_id=~{project_id} --format=json \
            "$(cat acceptance_checks.sql)" | jq '.[0] | map_values(tonumber)' > acceptance_checks.json
        cat acceptance_checks.json

        failures=$(jq -r 'to_entries[]
            | select(.key | IN("mismatch_outside_provenance", "below_inside_provenance", "mapped_not_in_vat",
                               "vat_carriers_without_row", "vids_with_duplicate_person"))
            | select(.value != 0) | "\(.key)=\(.value)"' acceptance_checks.json)
        if [[ -n "${failures}" ]]
        then
            echo "Acceptance checks failed: ${failures//$'\n'/, }" >&2
            exit 1
        fi
    >>>

    runtime {
        docker: variants_docker
    }

    output {
        File acceptance_checks = "acceptance_checks.json"
    }
}
