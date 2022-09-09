version 1.0

workflow GvsCallsetStatistics {
    input {
        String project_id
        String dataset_name
        String filter_set_name
        String extract_prefix
        String output_statistics_table = "~{extract_prefix}_sample_metrics"
        String output_aggregate_statistics_table = "~{extract_prefix}_sample_metrics_agg"
    }

    call CreateTables {
        input:
            project_id = project_id,
            dataset_name = dataset_name,
            output_statistics_table = output_statistics_table,
            output_aggregate_statistics_table = output_aggregate_statistics_table
    }

    # Only collect statistics for the autosomal chromosomes, the first 22 in our location scheme.
    scatter(chrom in range(22)) {
        call CollectStatisticsForChromosome {
            input:
                go = CreateTables.done,
                project_id = project_id,
                dataset_name = dataset_name,
                filter_set_name = filter_set_name,
                extract_prefix = extract_prefix,
                output_statistics_table = output_statistics_table,
                chromosome = chrom + 1 # 0-based ==> 1-based
        }
    }

    call AggregateStatisticsAcrossChromosomes {
        input:
            go = CollectStatisticsForChromosome.done[0],
            project_id = project_id,
            dataset_name = dataset_name,
            filter_set_name = filter_set_name,
            extract_prefix = extract_prefix,
            output_statistics_table = output_statistics_table,
            output_aggregate_statistics_table = output_aggregate_statistics_table
    }
}

task CreateTables {
    input {
        String project_id
        String dataset_name
        String output_statistics_table
        String output_aggregate_statistics_table
    }
    command <<<
        set -o errexit -o nounset -o xtrace -o pipefail

        set +o errexit
        bq --project_id=~{project_id} show ~{dataset_name}.~{output_statistics_table}
        BQ_SHOW_RC=$?
        set -o errexit

        if [ $BQ_SHOW_RC -eq 0 ]; then
            echo "Output statistics table already exists, exiting."
            exit 1
        fi

        set +o errexit
        bq --project_id=~{project_id} show ~{dataset_name}.~{output_aggregate_statistics_table}
        BQ_SHOW_RC=$?
        set -o errexit

        if [ $BQ_SHOW_RC -eq 0 ]; then
            echo "Output aggregate statistics table already exists, exiting."
            exit 1
        fi

        cat > schema.json <<FIN
        [
            {
                "name": "filter_set_name",
                "type": "STRING",
                "mode": "REQUIRED"
            },
            {
                "name": "sample_id",
                "type": "INTEGER",
                "mode": "REQUIRED"
            },
            {
                "name": "variant_entries",
                "type": "INTEGER",
                "mode": "REQUIRED"
            },
            {
                "name": "del_count",
                "type": "INTEGER",
                "mode": "REQUIRED"
            },
            {
                "name": "ins_count",
                "type": "INTEGER",
                "mode": "REQUIRED"
            },
            {
                "name": "snp_count",
                "type": "INTEGER",
                "mode": "REQUIRED"
            },
            {
                "name": "ti_count",
                "type": "INTEGER",
                "mode": "REQUIRED"
            },
            {
                "name": "tv_count",
                "type": "INTEGER",
                "mode": "REQUIRED"
            },
            {
                "name": "snp_het_count",
                "type": "INTEGER",
                "mode": "REQUIRED"
            },
            {
                "name": "snp_homvar_count",
                "type": "INTEGER",
                "mode": "REQUIRED"
            },
            {
                "name": "indel_het_count",
                "type": "INTEGER",
                "mode": "REQUIRED"
            },
            {
                "name": "indel_homvar_count",
                "type": "INTEGER",
                "mode": "REQUIRED"
            },
            {
                "name": "singleton",
                "type": "INTEGER",
                "mode": "REQUIRED"
            },
            {
                "name": "pass_qc",
                "type": "STRING",
                "mode": "OPTIONAL"
            },

        ]
        FIN

        bq mk --table ~{project_id}:~{dataset_name}.~{output_statistics_table} schema.json
        bq mk --table ~{project_id}:~{dataset_name}.~{output_aggregate_statistics_table} schema.json
    >>>
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:ah_var_store_2022_08_22"
    }
    output {
        Boolean done = true
    }
}

task CollectStatisticsForChromosome {
    input {
        Boolean go = true
        String project_id
        String dataset_name
        String filter_set_name
        String extract_prefix
        String output_statistics_table
        Int chromosome
    }
    command <<<
        set -o errexit -o nounset -o xtrace -o pipefail

        echo "project_id = ~{project_id}" > ~/.bigqueryrc

        bq query --location=US --project_id=~{project_id} --format=csv --use_legacy_sql=false '
        CREATE TEMPORARY FUNCTION titv(ref STRING, allele STRING)
        RETURNS STRING
            LANGUAGE js AS """
                if ( ref.length > 1 || allele.length > 1) {
                    return "other";
                } else if ( (ref == "A" && allele == "G") ||
                            (ref == "G" && allele == "A") ||
                            (ref == "C" && allele == "T") ||
                            (ref == "T" && allele == "C") ) {
                    return "ti";
                } else {
                    return "tv";
                }
        """;

        CREATE TEMPORARY FUNCTION type(ref STRING, allele STRING, gt_str STRING)
        RETURNS STRING
            LANGUAGE js AS """

        alts = allele.split(",")

        // get the the non-reference allele indexes
        ai = gt_str.replace("|","/").split("/").filter(i => i != "0");

        // the the distinct set of lengths of the alternates
        alt_lengths = new Set(ai.map(i => alts[parseInt(i)-1].length))

        if (alt_lengths.size > 1) {
            return "complex"
        } else {
            // get first (only) element
            al = alt_lengths.keys().next().value

            if ( ref.length == al && al == 1) {
                return "snp"
            } else if (ref.length > al) {
                return "del"
            } else if (ref.length < al) {
                return "ins"
            } else {
                return "other"
            }
        }
        """;

        INSERT `~{project_id}.~{dataset_name}.~{output_statistics_table}` (
            filter_set_name,
            sample_id,
            variant_entries,
            del_count,
            ins_count,
            snp_count,
            ti_count,
            tv_count,
            snp_het_count,
            snp_homvar_count,
            indel_het_count,
            indel_homvar_count,
            singleton,
            pass_qc
        )
        SELECT "~{filter_set_name}" filter_set_name,
               sample_id,
               count(1) variant_entries,
               SUM(CASE WHEN type = "del" THEN 1 ELSE 0 END) del_count,
               SUM(CASE WHEN type = "ins" THEN 1 ELSE 0 END) ins_count,
               SUM(CASE WHEN type = "snp" THEN 1 ELSE 0 END) snp_count,
               SUM(CASE WHEN type = "snp" AND titv = "ti" THEN 1 ELSE 0 END) ti_count, # TODO: minimize alleles
               SUM(CASE WHEN type = "snp" AND titv = "tv" THEN 1 ELSE 0 END) tv_count, # TODO: minimize alleles
               SUM(CASE WHEN type = "snp" AND gt_type = "het" THEN 1 ELSE 0 END) snp_het_count,
               SUM(CASE WHEN type = "snp" AND gt_type = "homvar" THEN 1 ELSE 0 END) snp_homvar_count,
               SUM(CASE WHEN type IN ("ins","del") AND gt_type = "het" THEN 1 ELSE 0 END) indel_het_count,
               SUM(CASE WHEN type IN ("ins","del") AND gt_type = "homvar" THEN 1 ELSE 0 END) indel_homvar_count,
               COUNTIF(not in_gnomad) singleton,
               null AS pass_qc
        FROM (
            SELECT sample_id,
                   type(ref, alt, call_GT) as type,
                   CASE WHEN INSTR(call_GT, "0") > 0 THEN "het" ELSE "homvar" END as gt_type,
                   titv(ref, alt) as titv,
                   CASE WHEN gnomad.location IS NULL THEN false ELSE true END in_gnomad
            FROM `~{project_id}.~{dataset_name}.~{extract_prefix}__VET_DATA` v
            LEFT JOIN `spec-ops-aou.gvs_public_reference_data.gnomad_v3_sites` gnomad ON (v.location = gnomad.location)
            WHERE call_GT != "./."
            AND v.location >= ~{chromosome}0000000000000
            AND v.location < ~{chromosome + 1}000000000000) GROUP BY 1,2

        '

    >>>
    output {
        Boolean done = true
    }
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:ah_var_store_2022_08_22"
    }
}

task AggregateStatisticsAcrossChromosomes {
    input {
        Boolean go
        String project_id
        String dataset_name
        String filter_set_name
        String extract_prefix
        String output_statistics_table
        String output_aggregate_statistics_table
    }
    command <<<
        bq query --location=US --project_id=~{project_id} --format=csv --use_legacy_sql=false '

        INSERT `~{project_id}.~{dataset_name}.~{output_aggregate_statistics_table}` (
            filter_set_name,
            sample_id,
            variant_entries,
            del_count,
            ins_count,
            snp_count,
            ti_count,
            tv_count,
            snp_het_count,
            snp_homvar_count,
            indel_het_count,
            indel_homvar_count,
            singleton,
            pass_qc
        )

        SELECT "~{filter_set_name}" filter_set_name,
                sample_id,
                SUM(variant_entries) variant_entries,
                SUM(CASE WHEN type = "del" THEN 1 ELSE 0 END) del_count,
                SUM(CASE WHEN type = "ins" THEN 1 ELSE 0 END) ins_count,
                SUM(CASE WHEN type = "snp" THEN 1 ELSE 0 END) snp_count,
                SUM(CASE WHEN type = "snp" AND titv = "ti" THEN 1 ELSE 0 END) ti_count, # TODO: minimize alleles
                SUM(CASE WHEN type = "snp" AND titv = "tv" THEN 1 ELSE 0 END) tv_count, # TODO: minimize alleles
                SUM(CASE WHEN type = "snp" AND gt_type = "het" THEN 1 ELSE 0 END) snp_het_count,
                SUM(CASE WHEN type = "snp" AND gt_type = "homvar" THEN 1 ELSE 0 END) snp_homvar_count,
                SUM(CASE WHEN type IN ("ins","del") AND gt_type = "het" THEN 1 ELSE 0 END) indel_het_count,
                SUM(CASE WHEN type IN ("ins","del") AND gt_type = "homvar" THEN 1 ELSE 0 END) indel_homvar_count,
                SUM(singleton) singleton,
                null AS pass_qc
        FROM `~{project_id}.~{dataset_name}.~{output_statistics_table}` GROUP BY 1,2

        '
    >>>
    output {
        Boolean done = true
    }
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:ah_var_store_2022_08_22"
    }
}
