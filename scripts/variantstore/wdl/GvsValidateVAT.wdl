version 1.0

import "GvsWarpTasks.wdl" as Tasks

workflow GvsValidateVatTable {
    input {
        String query_project_id
        String default_dataset
        String vat_table_name
        String? service_account_json_path
    }

    String fq_vat_table = "~{query_project_id}.~{default_dataset}.~{vat_table_name}"

    call GetBQTableLastModifiedDatetime {
        input:
            query_project = query_project_id,
            fq_table = fq_vat_table,
            service_account_json_path = service_account_json_path
    }

    call EnsureVatTableHasVariants {
        input:
            query_project_id = query_project_id,
            fq_vat_table = fq_vat_table,
            service_account_json_path = service_account_json_path,
            last_modified_timestamp = GetBQTableLastModifiedDatetime.last_modified_timestamp
    }

    call SpotCheckForExpectedTranscripts {
        input:
            query_project_id = query_project_id,
            fq_vat_table = fq_vat_table,
            service_account_json_path = service_account_json_path,
            last_modified_timestamp = GetBQTableLastModifiedDatetime.last_modified_timestamp
    }

    call SchemaOnlyOneRowPerNullTranscript {
        input:
            query_project_id = query_project_id,
            fq_vat_table = fq_vat_table,
            service_account_json_path = service_account_json_path,
            last_modified_timestamp = GetBQTableLastModifiedDatetime.last_modified_timestamp
    }

    call EnsemblTranscripts {
        input:
            query_project_id = query_project_id,
            fq_vat_table = fq_vat_table,
            service_account_json_path = service_account_json_path,
            last_modified_timestamp = GetBQTableLastModifiedDatetime.last_modified_timestamp
    }


    call NonzeroAcAn {
        input:
            query_project_id = query_project_id,
            fq_vat_table = fq_vat_table,
            service_account_json_path = service_account_json_path,
            last_modified_timestamp = GetBQTableLastModifiedDatetime.last_modified_timestamp
    }

    # once there is more than one check, they will be gathered into this workflow output, in the format
    # [{ValidationRule1: "PASS/FAIL Extra info from this test"},
    #  {ValidationRule2: "PASS/FAIL Extra from this test"}]

    output {
      Array[Map[String, String]] validation_results = [
        EnsureVatTableHasVariants.result,
        SpotCheckForExpectedTranscripts.result,
        SchemaOnlyOneRowPerNullTranscript.result,
        EnsemblTranscripts.result,
        NonzeroAcAn.result
      ]
    }
}

task GetBQTableLastModifiedDatetime {
    # because this is being used to determine if the data has changed, never use call cache
    meta {
        volatile: true
    }

    input {
        String query_project
        String fq_table
        String? service_account_json_path
    }

    String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

    # ------------------------------------------------
    # try to get the last modified date for the table in question; fail if something comes back from BigQuery
    # that isn't in the right format (e.g. an error)
    command <<<
        set -e

        if [ ~{has_service_account_file} = 'true' ]; then
            gsutil cp ~{service_account_json_path} local.service_account.json
            gcloud auth activate-service-account --key-file=local.service_account.json
            gcloud config set project ~{query_project}
        fi

        echo "project_id = ~{query_project}" > ~/.bigqueryrc

        # bq needs the project name to be separate by a colon
        DATASET_TABLE_COLON=$(echo ~{fq_table} | sed 's/\./:/')

        LASTMODIFIED=$(bq --location=US --project_id=~{query_project} --format=json show ${DATASET_TABLE_COLON} | python3 -c "import sys, json; print(json.load(sys.stdin)['lastModifiedTime']);")
        if [[ $LASTMODIFIED =~ ^[0-9]+$ ]]; then
            echo $LASTMODIFIED
        else
            exit 1
        fi
    >>>

    output {
        String last_modified_timestamp = read_string(stdout())
    }

    runtime {
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
        memory: "3 GB"
        disks: "local-disk 10 HDD"
        preemptible: 3
        cpu: 1
    }

}

task EnsureVatTableHasVariants {
    input {
        String query_project_id
        String fq_vat_table
        String? service_account_json_path
        String last_modified_timestamp
    }

    String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

    command <<<
        set -e
        if [ ~{has_service_account_file} = 'true' ]; then
            gsutil cp ~{service_account_json_path} local.service_account.json
            gcloud auth activate-service-account --key-file=local.service_account.json
            gcloud config set project ~{query_project_id}
        fi
        echo "project_id = ~{query_project_id}" > ~/.bigqueryrc

        bq query --nouse_legacy_sql --project_id=~{query_project_id} --format=csv 'SELECT COUNT (DISTINCT vid) AS count FROM ~{fq_vat_table}' > bq_variant_count.csv

        NUMVARS=$(python3 -c "csvObj=open('bq_variant_count.csv','r');csvContents=csvObj.read();print(csvContents.split('\n')[1]);")

        # if the result of the bq call and the csv parsing is a series of digits, then check that it isn't 0
        if [[ $NUMVARS =~ ^[0-9]+$ ]]; then
            if [[ $NUMVARS = "0" ]]; then
                echo "FAIL: The VAT table ~{fq_vat_table} has no variants in it." > validation_results.txt
            else
                echo "PASS: The VAT table ~{fq_vat_table} has $NUMVARS variants in it." > validation_results.txt
            fi
        # otherwise, something is off, so return the output from the bq query call
        else
            echo "Something went wrong. The attempt to count the variants returned: " $(cat bq_variant_count.csv) > validation_results.txt
        fi
    >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
        memory: "1 GB"
        preemptible: 3
        cpu: "1"
        disks: "local-disk 100 HDD"
    }
    # ------------------------------------------------
    # Output: {"Name of validation rule": "PASS/FAIL plus additional validation results"}
    output {
        Map[String, String] result = {"EnsureVatTableHasVariants": read_string('validation_results.txt')}
    }
}


task SpotCheckForExpectedTranscripts {
    input {
        String query_project_id
        String fq_vat_table
        String? service_account_json_path
        String last_modified_timestamp
    }

    String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

    command <<<
        set -e

        if [ ~{has_service_account_file} = 'true' ]; then
            gsutil cp ~{service_account_json_path} local.service_account.json
            gcloud auth activate-service-account --key-file=local.service_account.json
            gcloud config set project ~{query_project_id}
        fi
        echo "project_id = ~{query_project_id}" > ~/.bigqueryrc

        bq query --nouse_legacy_sql --project_id=~{query_project_id} --format=csv 'SELECT
            contig,
            position,
            vid,
            gene_symbol,
            variant_consequence
        FROM
            ~{fq_vat_table},
            UNNEST(consequence) AS variant_consequence
        WHERE
            contig = "chr19" AND
            position >= 35740407 AND
            position <= 35740469 AND
            variant_consequence NOT IN ("downstream_gene_variant","upstream_gene_variant") AND
            gene_symbol NOT IN ("IGFLR1","AD000671.2")' > bq_query_output.csv

        # get number of lines in bq query output
        NUMRESULTS=$(awk 'END{print NR}' bq_query_output.csv)

        # if the result of the query has any rows, that means there were unexpected transcripts at the
        # specified location, so report those back in the output
        if [[ $NUMRESULTS = "0" ]]; then
            echo "PASS: The VAT table ~{fq_vat_table} only has the expected transcripts at the tested location ('IGFLR1' and 'AD000671.2' in chromosome 19, between positions 35,740,407 - 35,740,469)." > validation_results.txt
        else
           echo "FAIL: The VAT table ~{fq_vat_table} had unexpected transcripts at the tested location: [csv output follows] " > validation_results.txt
            cat bq_query_output.csv >> validation_results.txt
        fi
    >>>

    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
        memory: "1 GB"
        preemptible: 3
        cpu: "1"
        disks: "local-disk 100 HDD"
    }

    output {
        Map[String, String] result = {"SpotCheckForExpectedTranscripts": read_string('validation_results.txt')}

    }
}

task EnsemblTranscripts {
    input {
        String query_project_id
        String fq_vat_table
        String? service_account_json_path
        String last_modified_timestamp
    }
    # Every transcript_source is Ensembl or null

    String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

    command <<<
        if [ ~{has_service_account_file} = 'true' ]; then
            gsutil cp ~{service_account_json_path} local.service_account.json
            gcloud auth activate-service-account --key-file=local.service_account.json
            gcloud config set project ~{query_project_id}
        fi
        echo "project_id = ~{query_project_id}" > ~/.bigqueryrc

        bq query --nouse_legacy_sql --project_id=~{query_project_id} --format=csv
        'SELECT
            contig,
            position,
            vid,
            transcript,
            transcript_source
        FROM
            ~{fq_vat_table},
        WHERE
            transcript IS NOT NULL AND
            transcript_source != "Ensembl"' > bq_transcript_output.csv


            # get number of lines in bq query output
            NUMRESULTS=$(awk 'END{print NR}' bq_transcript_output.csv)

            # if the result of the query has any rows, that means there were unexpected transcripts (not from Ensembl), so report those back in the output
            if [[ $NUMRESULTS = "0" ]]; then
                echo "PASS: The VAT table ~{fq_vat_table} only has the expected Ensembl transcripts"  > validation_results.txt
            else
                echo "FAIL: The VAT table ~{fq_vat_table} had unexpected transcripts (not from Ensembl): [csv output follows] " > validation_results.txt
                cat bq_transcript_output.csv >> validation_results.txt
            fi
    >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
        memory: "1 GB"
        preemptible: 3
        cpu: "1"
        disks: "local-disk 100 HDD"
    }
    # ------------------------------------------------
    # Output: {"Name of validation rule": "PASS/FAIL plus additional validation results"}
    output {
        Map[String, String] result = {"EnsemblTranscripts": read_string('validation_results.txt')}
    }
}

task NonzeroAcAn {
    input {
        String query_project_id
        String fq_vat_table
        String? service_account_json_path
        String last_modified_timestamp
    }
    # No row has AC of zero or AN of zero.

    String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

    command <<<
        if [ ~{has_service_account_file} = 'true' ]; then
            gsutil cp ~{service_account_json_path} local.service_account.json
            gcloud auth activate-service-account --key-file=local.service_account.json
            gcloud config set project ~{query_project_id}
        fi
        echo "project_id = ~{query_project_id}" > ~/.bigqueryrc

        bq query --nouse_legacy_sql --project_id=~{query_project_id} --format=csv
        'SELECT
            contig,
            position,
            vid,
            gvs_all_ac,
            gvs_all_an
        FROM
            ~{fq_vat_table},
        WHERE
            gvs_all_ac IS NULL OR
            gvs_all_ac == 0 OR
            gvs_all_an IS NULL OR
            gvs_all_an == 0 ' > bq_ac_an_output.csv


            # get number of lines in bq query output
            NUMRESULTS=$(awk 'END{print NR}' bq_ac_an_output.csv)

            # if the result of the query has any rows, that means there were unexpected rows with either an AC of zero or AN of zero, so report those back in the output
            if [[ $NUMRESULTS = "0" ]]; then
                echo "PASS: The VAT table ~{fq_vat_table} only has no rows with AC of zero or AN of zero"  > validation_results.txt
            else
                echo "FAIL: The VAT table ~{fq_vat_table} had unexpected rows with AC of zero or AN of zero: [csv output follows] " > validation_results.txt
                cat bq_ac_an_output.csv >> validation_results.txt
            fi
    >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
        memory: "1 GB"
        preemptible: 3
        cpu: "1"
        disks: "local-disk 100 HDD"
    }
    # ------------------------------------------------
    # Output: {"Name of validation rule": "PASS/FAIL plus additional validation results"}
    output {
        Map[String, String] result = {"NonzeroAcAn": read_string('validation_results.txt')}
    }
}

task GetBQTableLastModifiedDatetime {
    # because this is being used to determine if the data has changed, never use call cache
    meta {
        volatile: true
    }

task SchemaOnlyOneRowPerNullTranscript {
    input {
        String query_project_id
        String fq_vat_table
        String? service_account_json_path
        String last_modified_timestamp
    }

    String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

    command <<<
        set -e

        if [ ~{has_service_account_file} = 'true' ]; then
            gsutil cp ~{service_account_json_path} local.service_account.json
            gcloud auth activate-service-account --key-file=local.service_account.json
            gcloud config set project ~{query_project_id}
        fi
        echo "project_id = ~{query_project_id}" > ~/.bigqueryrc

        bq query --nouse_legacy_sql --project_id=~{query_project_id} --format=csv 'SELECT
            vid,
            COUNT(vid) AS num_rows
        FROM
            ~{fq_vat_table}
        WHERE
            transcript_source is NULL AND
            transcript is NULL
        GROUP BY vid
        HAVING num_rows = 1' > bq_variant_count.csv

        # get number of lines in bq query output
        NUMRESULTS=$(awk 'END{print NR}' bq_variant_count.csv)

        # if the result of the query has any rows, that means there were vids will null transcripts and multiple
        # rows in the VAT, which should not be the case
        if [[ $NUMRESULTS = "0" ]]; then
            echo "PASS: The VAT table ~{fq_vat_table} only has 1 row per vid with a null transcript" > validation_results.txt
        else
            echo "FAIL: The VAT table ~{fq_vat_table} had at least one vid with a null transcript and more than one row: [csv output follows] " > validation_results.txt
            cat bq_variant_count.csv >> validation_results.txt
        fi
    >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
        memory: "1 GB"
        preemptible: 3
        cpu: "1"
        disks: "local-disk 100 HDD"
    }
    # ------------------------------------------------
    # Output: {"Name of validation rule": "PASS/FAIL plus additional validation results"}
    output {
        Map[String, String] result = {"SchemaOnlyOneRowPerNullTranscript": read_string('validation_results.txt')}
    }
}
