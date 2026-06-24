version 1.0

import "../wdl/GvsUtils.wdl" as Utils

workflow GvsCreateVATFilesFromBigQuery {
    input {
        String project_id
        String? git_branch_or_tag
        String? git_hash
        String dataset_name
        String vat_table_name

        String output_path
        Int? merge_vcfs_disk_size_override
        String? cloud_sdk_docker
        String? cloud_sdk_slim_docker
    }

    Array[String] contig_array = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"]

    if (!defined(git_hash) || !defined(cloud_sdk_docker) || !defined(cloud_sdk_slim_docker)) {
        call Utils.GetToolVersions {
            input:
                git_branch_or_tag = git_branch_or_tag,
        }
    }

    String effective_cloud_sdk_docker = select_first([cloud_sdk_docker, GetToolVersions.cloud_sdk_docker])
    String effective_cloud_sdk_slim_docker = select_first([cloud_sdk_slim_docker, GetToolVersions.cloud_sdk_slim_docker])
    String effective_git_hash = select_first([git_hash, GetToolVersions.git_hash])

    scatter(i in range(length(contig_array)) ) {
        call BigQueryExportVat {
            input:
                contig = contig_array[i],
                project_id = project_id,
                dataset_name = dataset_name,
                output_path = output_path,
                vat_table = vat_table_name,
                cloud_sdk_docker = effective_cloud_sdk_docker,
        }
    }

    call MergeVatTSVs {
        input:
            export_done = BigQueryExportVat.done,
            contig_array = contig_array,
            output_path = output_path,
            merge_vcfs_disk_size_override = merge_vcfs_disk_size_override,
            project_id = project_id,
            dataset_name = dataset_name,
            vat_table = vat_table_name,
            cloud_sdk_slim_docker = effective_cloud_sdk_slim_docker,
    }

    output {
        File final_tsv_file = MergeVatTSVs.tsv_file
        File final_tsv_file_index = MergeVatTSVs.tsv_file_index
        String recorded_git_hash = effective_git_hash
    }
}

################################################################################


task BigQueryExportVat {
    input {
        String contig
        String project_id
        String dataset_name
        String vat_table
        String output_path
        String cloud_sdk_docker
    }

    String export_path = output_path + "export/" + contig + "/*.tsv.gz"

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        echo "project_id = ~{project_id}" > ~/.bigqueryrc

        # Write out a query to a file using an uninterpolated here doc. The single quotes around the delimiter tell bash
        # not to interpolate within the here doc. Without the single quotes bash would evaluate anything that looked
        # like a bash expression in the query body, which would be problematic here for backticks.
        cat > query.sql <<'FIN'

        DECLARE dynamic_vat_query STRING;
        DECLARE vat_query STRING;
        DECLARE export_query STRING;

        -- `dynamic_vat_query` is a concatenation of three expressions: two literal strings sandwiching an inner query
        -- of `INFORMATION_SCHEMA`. This inner query determines the names and types of columns in the VAT table.
        -- Depending on whether a column type is primitive or `ARRAY`, the `CASE` statement returns a query fragment
        -- that either returns the value of the field (for a primitive) or an expression that `STRING_AGG`s the array
        -- values into a single comma-separated string. Note this logic is building a query string for determining the
        -- shape of the VAT table based on the contents of `INFORMATION_SCHEMA`, but that query is not actually executed
        -- until the *second* `EXECUTE IMMEDIATE` statement below.
        SET dynamic_vat_query = """

        SELECT 'SELECT ' ||
        (
            SELECT
            STRING_AGG(
                CASE
                WHEN data_type LIKE 'ARRAY%' THEN FORMAT('(SELECT STRING_AGG(CAST(x AS STRING), ",") FROM UNNEST(%s) x) AS %s', column_name, column_name)
                ELSE column_name
                END,
            ', ')
            FROM
                `~{project_id}.~{dataset_name}.INFORMATION_SCHEMA.COLUMNS`
            WHERE
                table_name = '~{vat_table}'
        ) ||
        ' FROM `~{project_id}.~{dataset_name}.~{vat_table}` WHERE contig = "~{contig}" ORDER BY position'

        """;

        -- Run the query constructed above to materialize the string for querying the VAT table and store the result in
        -- the variable `vat_query`.
        EXECUTE IMMEDIATE dynamic_vat_query INTO vat_query;

        -- Build a final export query by concatenating the VAT query built above with some export syntax.
        -- The combination of tab delimiter and GZIP compression creates tsv.gz files.
        SET export_query = """

        EXPORT DATA OPTIONS(
            uri="~{export_path}",
            format="CSV",
            compression="GZIP",
            overwrite=true,
            header=false,
            field_delimiter="\t") AS

        """ || vat_query;

        -- Print out the export query for diagnostic purposes.
        SELECT export_query;

        -- Execute the export query.
        EXECUTE IMMEDIATE export_query;
        FIN

        # Feed the query file to `bq` as input.
        # bq query --max_rows check: ok export
        bq --apilog=false query --nouse_legacy_sql --project_id=~{project_id} < query.sql > export_query.sql

    >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: cloud_sdk_docker
        memory: "2 GB"
        preemptible: 3
        cpu: "1"
        disks: "local-disk 100 HDD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        Boolean done = true
        File query = "query.sql"
        File export_query = "export_query.sql"
    }
}

task MergeVatTSVs {
    input {
        Array[Boolean] export_done
        Array[String] contig_array
        String project_id
        String dataset_name
        String vat_table
        String output_path

        Int? merge_vcfs_disk_size_override
        String cloud_sdk_slim_docker
    }

    File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"

    # going large with the default to make gsutil -m cp really zippy
    Int disk_size = if (defined(merge_vcfs_disk_size_override)) then select_first([merge_vcfs_disk_size_override]) else 500

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        # Kick off the monitoring script
        bash ~{monitoring_script} > monitoring.log &

        apt-get update
        apt-get install --assume-yes tabix

        # custom function to prepend the current datetime to an echo statement "borrowed" from ExtractAnAcAfFromVCF
        echo_date () { echo "`date "+%Y/%m/%d %H:%M:%S"` $1"; }

        mkdir TSVs
        contigs=( ~{sep=' ' contig_array} )
        files="header.gz"

        cat > query.sql <<'FIN'

        SELECT STRING_AGG(column_name, "\t") FROM
            `~{project_id}.~{dataset_name}.INFORMATION_SCHEMA.COLUMNS`
        WHERE
            table_name = '~{vat_table}'

        FIN

        # Feed the query file to `bq` as input.
        # bq query --max_rows check: ok 1 row
        bq --apilog=false query --nouse_legacy_sql --format csv --project_id=~{project_id} < query.sql | sed 1d > header.tsv

        echo_date "looping over contigs: $contigs"
        for i in "${contigs[@]}"
        do
            echo_date "copying files from ~{output_path}export/$i/*.tsv.gz"
            gcloud storage cp ~{output_path}export/$i/*.tsv.gz TSVs/
            echo_date "concatenating local tsv.gz files"

            # the internet says that * is deterministic, see https://serverfault.com/questions/122737/in-bash-are-wildcard-expansions-guaranteed-to-be-in-order
            cat TSVs/*.tsv.gz > vat_$i.tsv.gz

            echo_date "removing now concatenated files"
            rm TSVs/*.tsv.gz
            files="$files vat_$i.tsv.gz"
        done

        echo_date "making header.gz"
        # NOTE: Contents of tsvs exported from BigQuery are tab-separated, the header must also be tab-separated!
        cat header.tsv | gzip > header.gz

        echo_date "bgzipping concatenated file"
        cat $(echo $files) | gunzip | bgzip -@ $(nproc) > vat_complete.bgz.tsv.gz
        rm -f $(echo $files)

        echo_date "indexing with tabix"
        contig_col=$(head -1 header.tsv | tr '\t' '\n' | grep -n "^contig$" | cut -d: -f1)
        position_col=$(head -1 header.tsv | tr '\t' '\n' | grep -n "^position$" | cut -d: -f1)
        echo_date "contig column: $contig_col, position column: $position_col"
        tabix -s "$contig_col" -b "$position_col" -e "$position_col" -S 1 vat_complete.bgz.tsv.gz

        echo_date "copying bgzipped file and index to ~{output_path}"
        gcloud storage cp vat_complete.bgz.tsv.gz ~{output_path}
        gcloud storage cp vat_complete.bgz.tsv.gz.tbi ~{output_path}
    >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: cloud_sdk_slim_docker
        memory: "8 GB"
        preemptible: 3
        cpu: 8
        disks: "local-disk ~{disk_size} SSD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        File tsv_file = "vat_complete.bgz.tsv.gz"
        File tsv_file_index = "vat_complete.bgz.tsv.gz.tbi"
        File monitoring_log = "monitoring.log"
    }
}
