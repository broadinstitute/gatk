version 1.0

import "GvsUtils.wdl" as Utils

workflow GvsCallsetStatistics {
    input {
        String? git_branch_or_tag
        String? git_hash
        String project_id
        String dataset_name
        String filter_set_name
        Int num_chunks = 2
        String extract_prefix
        String metrics_table = "~{extract_prefix}_sample_metrics"
        String aggregate_metrics_table = "~{extract_prefix}_sample_metrics_aggregate"
        String statistics_table = "~{extract_prefix}_statistics"
        String? basic_docker
        String? cloud_sdk_docker
    }

    # Always call `GetToolVersions` to get the git hash for this run as this is a top-level-only WDL (i.e. there are
    # no calling WDLs that might supply `git_hash`).
    if (!defined(git_hash) || !defined(basic_docker) || !defined(cloud_sdk_docker)) {
      call Utils.GetToolVersions {
          input:
              git_branch_or_tag = git_branch_or_tag,
      }
    }

    String effective_basic_docker = select_first([basic_docker, GetToolVersions.basic_docker])
    String effective_cloud_sdk_docker = select_first([cloud_sdk_docker, GetToolVersions.cloud_sdk_docker])
    String effective_git_hash = select_first([git_hash, GetToolVersions.git_hash])

    if (num_chunks <= 0 || num_chunks > 50) {
        call Utils.TerminateWorkflow as NumChunksError {
            input:
                message = "The input parameter 'num_chunks' must be >= 1 and < 50!",
                basic_docker = effective_basic_docker,
        }
    }

    call Utils.ValidateFilterSetName {
        input:
            project_id = project_id,
            fq_filter_set_info_table = "~{project_id}.~{dataset_name}.filter_set_info",
            filter_set_name = filter_set_name,
            cloud_sdk_docker = effective_cloud_sdk_docker
    }

    call CreateTables {
        input:
            go = ValidateFilterSetName.done,
            project_id = project_id,
            dataset_name = dataset_name,
            metrics_table = metrics_table,
            aggregate_metrics_table = aggregate_metrics_table,
            statistics_table = statistics_table,
            cloud_sdk_docker = effective_cloud_sdk_docker,
    }

    # Only collect statistics for the autosomal chromosomes, the first 22 in our location scheme.
    scatter(chrom in range(22)) {
        call CollectMetricsForChromosome {
            input:
                go = CreateTables.done,
                project_id = project_id,
                dataset_name = dataset_name,
                filter_set_name = filter_set_name,
                num_chunks = num_chunks,
                extract_prefix = extract_prefix,
                metrics_table = metrics_table,
                cloud_sdk_docker = effective_cloud_sdk_docker,
                chromosome = chrom + 1 # 0-based ==> 1-based,
        }
    }

    call AggregateMetricsAcrossChromosomes {
        input:
            go = CollectMetricsForChromosome.done[0],
            project_id = project_id,
            dataset_name = dataset_name,
            filter_set_name = filter_set_name,
            extract_prefix = extract_prefix,
            metrics_table = metrics_table,
            aggregate_metrics_table = aggregate_metrics_table,
            cloud_sdk_docker = effective_cloud_sdk_docker,
    }

    call CollectStatistics {
        input:
            go = AggregateMetricsAcrossChromosomes.done,
            project_id = project_id,
            dataset_name = dataset_name,
            filter_set_name = filter_set_name,
            extract_prefix = extract_prefix,
            metrics_table = metrics_table,
            aggregate_metrics_table = aggregate_metrics_table,
            statistics_table = statistics_table,
            cloud_sdk_docker = effective_cloud_sdk_docker,
    }

    call ExportToCSV {
        input:
          project_id = project_id,
          dataset_name = dataset_name,
          statistics_table = statistics_table,
          go = CollectStatistics.done,
          cloud_sdk_docker = effective_cloud_sdk_docker,
    }

    output {
        File callset_statistics = ExportToCSV.callset_statistics
        String recorded_git_hash = effective_git_hash
    }
}

task CreateTables {
    input {
        # Intentionally unused: this input exists solely to enforce task ordering — the upstream task's `done` output
        # is passed here to prevent this task from running until the upstream task has completed.
        #@ except: UnusedInput
        Boolean go = true
        String project_id
        String dataset_name
        String metrics_table
        String aggregate_metrics_table
        String statistics_table
        String cloud_sdk_docker
    }
    meta {
        # Always check that these tables exist
        volatile: true
    }
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        apk add jq

        set +o errexit
        bq --apilog=false --project_id=~{project_id} show ~{dataset_name}.~{metrics_table}
        BQ_SHOW_METRICS=$?
        set -o errexit

        set +o errexit
        bq --apilog=false --project_id=~{project_id} show ~{dataset_name}.~{aggregate_metrics_table}
        BQ_SHOW_METRICS_AGG=$?
        set -o errexit

        set +o errexit
        bq --apilog=false --project_id=~{project_id} show ~{dataset_name}.~{statistics_table}
        BQ_SHOW_STATISTICS=$?
        set -o errexit

        # Schemas extracted programatically: https://stackoverflow.com/a/66987934
        #
        # After cleaning up header and quotes:
        #
        # cat raw.json | jq -M '[.[]|{name,type,mode}' > clean.json

        cat > metrics_schema.json <<FIN
        [
          {
            "name": "filter_set_name",
            "type": "STRING",
            "mode": "NULLABLE"
          },
          {
            "name": "sample_id",
            "type": "INT64",
            "mode": "NULLABLE"
          },
          {
            "name": "chromosome",
            "type": "INT64",
            "mode": "NULLABLE"
          },
          {
            "name": "variant_entries",
            "type": "INT64",
            "mode": "NULLABLE"
          },
          {
            "name": "del_count",
            "type": "INT64",
            "mode": "NULLABLE"
          },
          {
            "name": "ins_count",
            "type": "INT64",
            "mode": "NULLABLE"
          },
          {
            "name": "snp_count",
            "type": "INT64",
            "mode": "NULLABLE"
          },
          {
            "name": "ti_count",
            "type": "INT64",
            "mode": "NULLABLE"
          },
          {
            "name": "tv_count",
            "type": "INT64",
            "mode": "NULLABLE"
          },
          {
            "name": "snp_het_count",
            "type": "INT64",
            "mode": "NULLABLE"
          },
          {
            "name": "snp_homvar_count",
            "type": "INT64",
            "mode": "NULLABLE"
          },
          {
            "name": "indel_het_count",
            "type": "INT64",
            "mode": "NULLABLE"
          },
          {
            "name": "indel_homvar_count",
            "type": "INT64",
            "mode": "NULLABLE"
          },
          {
            "name": "singleton",
            "type": "INT64",
            "mode": "NULLABLE"
          }
        ]
        FIN

        # The aggregate metrics schema is the same as the non-aggregate metrics schema delta the `chromosome` field.
        jq '[ .[] | select(.name != "chromosome") ]' metrics_schema.json > metrics_aggregate_schema.json

        cat > statistics_schema.json <<FIN
        [
          {
            "name": "sample_id",
            "type": "INT64",
            "mode": "NULLABLE"
          },
          {
            "name": "sample_name",
            "type": "STRING",
            "mode": "NULLABLE"
          },
          {
            "name": "del_count",
            "type": "INT64",
            "mode": "NULLABLE"
          },
          {
            "name": "ins_count",
            "type": "INT64",
            "mode": "NULLABLE"
          },
          {
            "name": "snp_count",
            "type": "INT64",
            "mode": "NULLABLE"
          },
          {
            "name": "singleton",
            "type": "INT64",
            "mode": "NULLABLE"
          },
          {
            "name": "ins_del_ratio",
            "type": "FLOAT64",
            "mode": "NULLABLE"
          },
          {
            "name": "ti_tv_ratio",
            "type": "FLOAT64",
            "mode": "NULLABLE"
          },
          {
            "name": "snp_het_homvar_ratio",
            "type": "FLOAT64",
            "mode": "NULLABLE"
          },
          {
            "name": "indel_het_homvar_ratio",
            "type": "FLOAT64",
            "mode": "NULLABLE"
          }
        ]
        FIN

        # Make any tables that need making
        if [ $BQ_SHOW_METRICS -ne 0 ]; then
            bq --apilog=false mk --table ~{project_id}:~{dataset_name}.~{metrics_table} metrics_schema.json
        fi

        if [ $BQ_SHOW_METRICS_AGG -ne 0 ]; then
            bq --apilog=false mk --table ~{project_id}:~{dataset_name}.~{aggregate_metrics_table} metrics_aggregate_schema.json
        fi

        if [ $BQ_SHOW_STATISTICS -ne 0 ]; then
            bq --apilog=false mk --table ~{project_id}:~{dataset_name}.~{statistics_table} statistics_schema.json
        fi
    >>>
    runtime {
        docker: cloud_sdk_docker
        disks: "local-disk 500 HDD"
    }
    output {
        Boolean done = true
    }
}

task CollectMetricsForChromosome {
    input {
        # Intentionally unused: this input exists solely to enforce task ordering — the upstream task's `done` output
        # is passed here to prevent this task from running until the upstream task has completed.
        #@ except: UnusedInput
        Boolean go = true
        String project_id
        String dataset_name
        String filter_set_name
        Int num_chunks = 2
        String extract_prefix
        String metrics_table
        Int chromosome
        String cloud_sdk_docker
    }
    meta {
        # This table is expected to exist and be empty. Always run to confirm it wasn't externally deleted, and don't
        # trust the leftover contents of any previous executions.
        volatile: true
    }
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        echo "project_id = ~{project_id}" > ~/.bigqueryrc

        # bq query --max_rows check: ok one row
        bq --apilog=false query --project_id=~{project_id} --format=csv --use_legacy_sql=false '
            SELECT COUNT(*) from `~{project_id}.~{dataset_name}.~{metrics_table}` WHERE chromosome = ~{chromosome}' | sed 1d > existing_row_count.txt

        existing_row_count=$(cat existing_row_count.txt)

        if [ $existing_row_count -gt 0 ]; then
            echo "Found $existing_row_count rows in '~{project_id}.~{dataset_name}.~{metrics_table}' for chromosome ~{chromosome}, exiting."
            exit 1
        fi

        chunk_endpoints=()
        if [ ~{num_chunks} -gt 1 ]; then
            # bq query Get the mid-point and min of the (found) locations for this chromosome
            bq --apilog=false query --project_id=~{project_id} --format=csv --use_legacy_sql=false '
                SELECT CAST((max(location) - min(location)) / ~{num_chunks} AS INT64), min(location)
                    FROM `~{project_id}.~{dataset_name}.~{extract_prefix}__VET_DATA`
                    WHERE location >= ~{chromosome}000000000000 and location < ~{chromosome + 1}000000000000' | sed 1d > locations.txt

            chunk_width=$(cut -f 1 -d ',' locations.txt)
            echo "chunk_width = $chunk_width"
            min_location=$(cut -f 2 -d ',' locations.txt)
            echo "min_location = $min_location"

            for i in $(seq 1 $((~{num_chunks} - 1))); do
                chunk_endpoints+=($((chunk_width * i + min_location)))
                echo "Added chunk endpoint ${chunk_endpoints[$((${#chunk_endpoints[@]} - 1))]}"
            done
        fi
        # Set the final endpoint to be beyond the start of the next chromosome
        chunk_endpoints+=(~{chromosome + 1}000000000000)
        echo "Added final chunk endpoint (the beginning of the next chromosome) = ~{chromosome + 1}000000000000"

        echo "Number of chunk endpoints: ${#chunk_endpoints[@]}"

        # Iterate over all chunks
        start_location=~{chromosome}000000000000
        for ((i=0; i<${#chunk_endpoints[@]}; i++)); do
            end_location=${chunk_endpoints[i]}
            echo "Running query for >= $start_location to < $end_location"

            # bq query --max_rows check: ok insert (elaborate one)
            bq --apilog=false query --project_id=~{project_id} --use_legacy_sql=false '
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

            INSERT `~{project_id}.~{dataset_name}.~{metrics_table}` (
                filter_set_name,
                sample_id,
                chromosome,
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
                singleton
            )
            SELECT "~{filter_set_name}" filter_set_name,
                   sample_id,
                   ~{chromosome},
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
                   COUNTIF(not in_gnomad) singleton
            FROM (
                SELECT sample_id,
                       type(ref, alt, call_GT) as type,
                       CASE WHEN INSTR(call_GT, "0") > 0 THEN "het" ELSE "homvar" END as gt_type,
                       titv(ref, alt) as titv,
                       CASE WHEN gnomad.location IS NULL THEN false ELSE true END in_gnomad
                FROM `~{project_id}.~{dataset_name}.~{extract_prefix}__VET_DATA` v
                LEFT JOIN `aou-genomics-curation-prod.gvs_public_reference_data.gnomad_v3_sites` gnomad ON (v.location = gnomad.location)
                WHERE call_GT != "./."
                AND v.location >= '$start_location'
                AND v.location < '$end_location') GROUP BY 1,2
            '
            start_location=$end_location
        done

    >>>
    output {
        Boolean done = true
    }
    runtime {
        docker: cloud_sdk_docker
        disks: "local-disk 500 HDD"
    }
}

task AggregateMetricsAcrossChromosomes {
    input {
        # Intentionally unused: this input exists solely to enforce task ordering — the upstream task's `done` output
        # is passed here to prevent this task from running until the upstream task has completed.
        #@ except: UnusedInput
        Boolean go
        String project_id
        String dataset_name
        String filter_set_name
        String extract_prefix
        String metrics_table
        String aggregate_metrics_table
        String cloud_sdk_docker
    }
    meta {
        # This table is expected to exist and be empty. Always run to confirm it wasn't externally deleted, and don't
        # trust the leftover contents of any previous executions.
        volatile: true
    }
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        # bq query --max_rows check: ok one row
        bq --apilog=false query --project_id=~{project_id} --format=csv --use_legacy_sql=false '
            SELECT COUNT(*) from `~{project_id}.~{dataset_name}.~{aggregate_metrics_table}`
        ' | sed 1d > existing_row_count.txt

        existing_row_count=$(cat existing_row_count.txt)

        if [ $existing_row_count -gt 0 ]; then
            echo "Found $existing_row_count rows in '~{project_id}.~{dataset_name}.~{aggregate_metrics_table}', exiting."
            exit 1
        fi

        # bq query --max_rows check: ok insert
        bq --apilog=false query --project_id=~{project_id} --use_legacy_sql=false '
        INSERT `~{project_id}.~{dataset_name}.~{aggregate_metrics_table}` (
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
            singleton
        )
        SELECT "~{filter_set_name}" filter_set_name,
            sample_id,
            SUM(variant_entries) variant_entries,
            SUM(del_count) del_count,
            SUM(ins_count) ins_count,
            SUM(snp_count) snp_count,
            SUM(ti_count) ti_count,
            SUM(tv_count) tv_count,
            SUM(snp_het_count) snp_het_count,
            SUM(snp_homvar_count) snp_homvar_count,
            SUM(indel_het_count) indel_het_count,
            SUM(indel_homvar_count) indel_homvar_count,
            SUM(singleton) singleton
        FROM `~{project_id}.~{dataset_name}.~{metrics_table}` GROUP BY 1,2

        '
    >>>
    output {
        Boolean done = true
    }
    runtime {
        docker: cloud_sdk_docker
        disks: "local-disk 500 HDD"
    }
}

task CollectStatistics {
    input {
        # Intentionally unused: this input exists solely to enforce task ordering — the upstream task's `done` output
        # is passed here to prevent this task from running until the upstream task has completed.
        #@ except: UnusedInput
        Boolean go
        String project_id
        String dataset_name
        String filter_set_name
        String extract_prefix
        String metrics_table
        String aggregate_metrics_table
        String statistics_table
        String cloud_sdk_docker
    }
    meta {
        # This table is expected to exist and be empty. Always run to confirm it wasn't externally deleted, and don't
        # trust the leftover contents of any previous executions.
        volatile: true
    }
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        # bq query --max_rows check: ok one row
        bq --apilog=false query --project_id=~{project_id} --format=csv --use_legacy_sql=false '
            SELECT COUNT(*) from `~{project_id}.~{dataset_name}.~{statistics_table}`
        ' | sed 1d > existing_row_count.txt

        existing_row_count=$(cat existing_row_count.txt)

        if [ $existing_row_count -gt 0 ]; then
            echo "Found $existing_row_count rows in '~{project_id}.~{dataset_name}.~{statistics_table}', exiting."
            exit 1
        fi

        # bq query --max_rows check: ok insert
        bq --apilog=false query --project_id=~{project_id} --format=csv --use_legacy_sql=false '
        INSERT `~{project_id}.~{dataset_name}.~{statistics_table}` (
            sample_id,
            sample_name,
            del_count,
            ins_count,
            snp_count,
            singleton,
            ins_del_ratio,
            ti_tv_ratio,
            snp_het_homvar_ratio,
            indel_het_homvar_ratio
        )

        SELECT
          amt.sample_id,
          si.sample_name,
          del_count,
          ins_count,
          snp_count,
          singleton,
          (ins_count / del_count) as ins_del_ratio,
          (ti_count / tv_count) as ti_tv_ratio,
          (snp_het_count / snp_homvar_count) snp_het_homvar_ratio,
          (indel_het_count / indel_homvar_count) as indel_het_homvar_ratio
        FROM `~{project_id}.~{dataset_name}.~{aggregate_metrics_table}` amt
        JOIN `~{project_id}.~{dataset_name}.sample_info` si ON (amt.sample_id = si.sample_id)
        WHERE amt.filter_set_name = "~{filter_set_name}"
        AND si.withdrawn IS NULL
        ORDER BY 1

        '
    >>>
    output {
        Boolean done = true
    }
    runtime {
        docker: cloud_sdk_docker
        disks: "local-disk 500 HDD"
    }
}

task ExportToCSV {
    input {
        String project_id
        String dataset_name
        String statistics_table
        # Intentionally unused: this input exists solely to enforce task ordering — the upstream task's `done` output
        # is passed here to prevent this task from running until the upstream task has completed.
        #@ except: UnusedInput
        Boolean go
        String cloud_sdk_docker
    }
    meta {
        # Many upstream dependencies and this is inexpensive anyway
        volatile: true
    }
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        # bq query --max_rows check: max rows set to at least the number of samples
        bq --apilog=false query --use_legacy_sql=false --project_id=~{project_id} --format=csv --max_rows 1000000000 '

          SELECT * FROM `~{project_id}.~{dataset_name}.~{statistics_table}` ORDER BY SAMPLE_NAME

        ' > '~{statistics_table}.csv'
    >>>
    output {
        File callset_statistics = "~{statistics_table}.csv"
    }
    runtime {
        docker: cloud_sdk_docker
        disks: "local-disk 500 HDD"
    }
}
