version 1.0

workflow GvsValidateVat {
    input {
        String query_project_id
        String default_dataset
        String vat_table_name
    }

    String fq_vat_table = "~{query_project_id}.~{default_dataset}.~{vat_table_name}"

    call GetBQTableLastModifiedDatetime {
        input:
            query_project = query_project_id,
            fq_table = fq_vat_table
    }

    call EnsureVatTableHasVariants {
        input:
            query_project_id = query_project_id,
            fq_vat_table = fq_vat_table,
            last_modified_timestamp = GetBQTableLastModifiedDatetime.last_modified_timestamp
    }

    call SpotCheckForExpectedTranscripts {
        input:
            query_project_id = query_project_id,
            fq_vat_table = fq_vat_table,
            last_modified_timestamp = GetBQTableLastModifiedDatetime.last_modified_timestamp
    }

    call SchemaOnlyOneRowPerNullTranscript {
        input:
            query_project_id = query_project_id,
            fq_vat_table = fq_vat_table,
            last_modified_timestamp = GetBQTableLastModifiedDatetime.last_modified_timestamp
    }

    call SchemaNullTranscriptsExist {
        input:
            query_project_id = query_project_id,
            fq_vat_table = fq_vat_table,
            last_modified_timestamp = GetBQTableLastModifiedDatetime.last_modified_timestamp
    }

    call SchemaNoNullRequiredFields {
        input:
            query_project_id = query_project_id,
            fq_vat_table = fq_vat_table,
            last_modified_timestamp = GetBQTableLastModifiedDatetime.last_modified_timestamp
    }

    call SchemaPrimaryKey {
        input:
            query_project_id = query_project_id,
            fq_vat_table = fq_vat_table,
            last_modified_timestamp = GetBQTableLastModifiedDatetime.last_modified_timestamp
    }

    call SchemaEnsemblTranscripts {
        input:
            query_project_id = query_project_id,
            fq_vat_table = fq_vat_table,
            last_modified_timestamp = GetBQTableLastModifiedDatetime.last_modified_timestamp
    }

    call SchemaNonzeroAcAn {
        input:
            query_project_id = query_project_id,
            fq_vat_table = fq_vat_table,
            last_modified_timestamp = GetBQTableLastModifiedDatetime.last_modified_timestamp
    }

    call SubpopulationMax {
        input:
            query_project_id = query_project_id,
            fq_vat_table = fq_vat_table,
            last_modified_timestamp = GetBQTableLastModifiedDatetime.last_modified_timestamp
    }

    call SubpopulationAlleleCount {
        input:
            query_project_id = query_project_id,
            fq_vat_table = fq_vat_table,
            last_modified_timestamp = GetBQTableLastModifiedDatetime.last_modified_timestamp
    }

    call SubpopulationAlleleNumber {
        input:
            query_project_id = query_project_id,
            fq_vat_table = fq_vat_table,
            last_modified_timestamp = GetBQTableLastModifiedDatetime.last_modified_timestamp
    }

    call ClinvarSignificance {
        input:
            query_project_id = query_project_id,
            fq_vat_table = fq_vat_table,
            last_modified_timestamp = GetBQTableLastModifiedDatetime.last_modified_timestamp
    }

    call SchemaAAChangeAndExonNumberConsistent {
        input:
            query_project_id = query_project_id,
            fq_vat_table = fq_vat_table,
            last_modified_timestamp = GetBQTableLastModifiedDatetime.last_modified_timestamp
    }

    call SpotCheckForAAChangeAndExonNumberConsistency {
         input:
            query_project_id = query_project_id,
            fq_vat_table = fq_vat_table,
            last_modified_timestamp = GetBQTableLastModifiedDatetime.last_modified_timestamp
    }

    call GenerateFinalReport {
        input:
            pf_flags = [
                       EnsureVatTableHasVariants.pass,
                       SpotCheckForExpectedTranscripts.pass,
                       SchemaOnlyOneRowPerNullTranscript.pass,
                       SchemaNullTranscriptsExist.pass,
                       SchemaNoNullRequiredFields.pass,
                       SchemaPrimaryKey.pass,
                       SchemaEnsemblTranscripts.pass,
                       SchemaNonzeroAcAn.pass,
                       SubpopulationMax.pass,
                       SubpopulationAlleleCount.pass,
                       SubpopulationAlleleNumber.pass,
                       ClinvarSignificance.pass,
                       SchemaAAChangeAndExonNumberConsistent.pass,
                       SpotCheckForAAChangeAndExonNumberConsistency.pass
                       ],
            validation_names = [
                               EnsureVatTableHasVariants.name,
                               SpotCheckForExpectedTranscripts.name,
                               SchemaOnlyOneRowPerNullTranscript.name,
                               SchemaNullTranscriptsExist.name,
                               SchemaNoNullRequiredFields.name,
                               SchemaPrimaryKey.name,
                               SchemaEnsemblTranscripts.name,
                               SchemaNonzeroAcAn.name,
                               SubpopulationMax.name,
                               SubpopulationAlleleCount.name,
                               SubpopulationAlleleNumber.name,
                               ClinvarSignificance.name,
                               SchemaAAChangeAndExonNumberConsistent.name,
                               SpotCheckForAAChangeAndExonNumberConsistency.name
                               ],
            validation_results = [
                                 EnsureVatTableHasVariants.result,
                                 SpotCheckForExpectedTranscripts.result,
                                 SchemaOnlyOneRowPerNullTranscript.result,
                                 SchemaNullTranscriptsExist.result,
                                 SchemaNoNullRequiredFields.result,
                                 SchemaPrimaryKey.result,
                                 SchemaEnsemblTranscripts.result,
                                 SchemaNonzeroAcAn.result,
                                 SubpopulationMax.result,
                                 SubpopulationAlleleCount.result,
                                 SubpopulationAlleleNumber.result,
                                 ClinvarSignificance.result,
                                 SchemaAAChangeAndExonNumberConsistent.result,
                                 SpotCheckForAAChangeAndExonNumberConsistency.result
                                 ]
    }

    output {
        String validation_results = GenerateFinalReport.results
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
    }


    # ------------------------------------------------
    # try to get the last modified date for the table in question; fail if something comes back from BigQuery
    # that isn't in the right format (e.g. an error)
    command <<<
        set -e

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
        String last_modified_timestamp
    }

    String pf_file = "pf.txt"
    String results_file = "results.txt"

    command <<<
        set -e
        echo "project_id = ~{query_project_id}" > ~/.bigqueryrc

        bq query --nouse_legacy_sql --project_id=~{query_project_id} --format=csv 'SELECT COUNT (DISTINCT vid) AS count FROM ~{fq_vat_table}' > bq_variant_count.csv

        NUMVARS=$(python3 -c "csvObj=open('bq_variant_count.csv','r');csvContents=csvObj.read();print(csvContents.split('\n')[1]);")

        echo "false" > ~{pf_file}
        # if the result of the bq call and the csv parsing is a series of digits, then check that it isn't 0
        if [[ $NUMVARS =~ ^[0-9]+$ ]]; then
            if [[ $NUMVARS = "0" ]]; then
                echo "The VAT table ~{fq_vat_table} has NO variants in it." > ~{results_file}
            else
                echo "The VAT table ~{fq_vat_table} has $NUMVARS variants in it." > ~{results_file}
                echo "true" > ~{pf_file}
            fi
        # otherwise, something is off, so return the output from the bq query call
        else
            echo "Something went wrong. The attempt to count the variants returned: " $(cat bq_variant_count.csv) >&2
            exit 1
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
        Boolean pass = read_boolean(pf_file)
        String name = "EnsureVatTableHasVariants"
        String result = read_string(results_file)
    }
}


task SpotCheckForExpectedTranscripts {
    input {
        String query_project_id
        String fq_vat_table
        String last_modified_timestamp
    }

    String pf_file = "pf.txt"
    String results_file = "results.txt"

    command <<<
        set -e

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

        echo "false" > ~{pf_file}
        # if the result of the query has any rows, that means there were unexpected transcripts at the
        # specified location, so report those back in the output
        if [[ $NUMRESULTS = "0" ]]; then
            echo "The VAT table ~{fq_vat_table} only has the expected transcripts at the tested location ('IGFLR1' and 'AD000671.2' in chromosome 19, between positions 35,740,407 - 35,740,469)." > ~{results_file}
            echo "true" > ~{pf_file}
        else
            echo "The VAT table ~{fq_vat_table} had unexpected transcripts at the tested location: [csv output follows] " > ~{results_file}
            cat bq_query_output.csv >> ~{results_file}
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
        Boolean pass = read_boolean(pf_file)
        String name = "SpotCheckForExpectedTranscripts"
        String result = read_string(results_file)
    }
}

task SchemaNoNullRequiredFields {
    input {
        String query_project_id
        String fq_vat_table
        String last_modified_timestamp
    }
    # No non-nullable fields contain null values

    String pf_file = "pf.txt"
    String results_file = "results.txt"

    command <<<
        echo "project_id = ~{query_project_id}" > ~/.bigqueryrc

        # non-nullable fields: vid, contig, position, ref_allele, alt_allele, gvs_all_ac, gvs_all_an, gvs_all_af, variant_type, genomic_location

        bq query --nouse_legacy_sql --project_id=~{query_project_id} --format=csv
        'SELECT
            contig,
            position,
            vid,
            concat(
              case(vid is null) when true then 'vid ' else '' end,
              case(contig is null) when true then 'contig ' else '' end,
              case(position is null) when true then 'position ' else '' end,
              case(ref_allele is null) when true then 'ref_allele ' else '' end,
              case(alt_allele is null) when true then 'alt_allele ' else '' end,
              case(gvs_all_ac is null) when true then 'gvs_all_ac ' else '' end,
              case(gvs_all_an is null) when true then 'gvs_all_an ' else '' end,
              case(gvs_all_af is null) when true then 'gvs_all_af ' else '' end,
              case(variant_type is null) when true then 'variant_type ' else '' end,
              case(genomic_location is null) when true then 'genomic_location ' else '' end
           ) AS null_fields
        FROM
            ~{fq_vat_table}
        WHERE
            vid IS NULL OR
            contig IS NULL OR
            position IS NULL OR
            ref_allele IS NULL OR
            alt_allele IS NULL OR
            gvs_all_ac IS NULL OR
            gvs_all_an IS NULL OR
            gvs_all_af IS NULL OR
            variant_type IS NULL OR
            genomic_location IS NULL' > bq_null_required_output.csv


            # get number of lines in bq query output
            NUMRESULTS=$(awk 'END{print NR}' bq_null_required_output.csv)

            echo "false" > ~{pf_file}
            # if the result of the query has any rows, that means there were unexpected null values in required fields, so report those back in the output
            if [[ $NUMRESULTS = "0" ]]; then
                echo "The VAT table ~{fq_vat_table} has no null values in required fields" > ~{results_file}
                echo "true" > ~{pf_file}
            else
                echo "The VAT table ~{fq_vat_table} had null values in required fields: [csv output follows] " > ~{results_file}
                cat bq_null_required_output.csv >> ~{results_file}
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
        Boolean pass = read_boolean(pf_file)
        String name = "SchemaNoNullRequiredFields"
        String result = read_string(results_file)
    }
}

task SchemaOnlyOneRowPerNullTranscript {
    input {
        String query_project_id
        String fq_vat_table
        String last_modified_timestamp
    }

    String pf_file = "pf.txt"
    String results_file = "results.txt"

    command <<<
        set -e

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
        HAVING num_rows > 1' > bq_variant_count.csv

        # get number of lines in bq query output
        NUMRESULTS=$(awk 'END{print NR}' bq_variant_count.csv)

        echo "false" > ~{pf_file}
        # if the result of the query has any rows, that means there were vids with null transcripts and multiple
        # rows in the VAT, which should not be the case
        if [[ $NUMRESULTS = "0" ]]; then
            echo "The VAT table ~{fq_vat_table} only has 1 row per vid with a null transcript" > ~{results_file}
            echo "true" > ~{pf_file}
        else
            echo "The VAT table ~{fq_vat_table} had at least one vid with a null transcript and more than one row: [csv output follows] " > ~{results_file}
            cat bq_variant_count.csv >> ~{results_file}
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
        Boolean pass = read_boolean(pf_file)
        String name = "SchemaOnlyOneRowPerNullTranscript"
        String result = read_string(results_file)
    }
}

task SchemaPrimaryKey {
    input {
        String query_project_id
        String fq_vat_table
        String last_modified_timestamp
    }
    # Each key combination (vid+transcript) is unique--confirms that primary key is enforced.

    String pf_file = "pf.txt"
    String results_file = "results.txt"

    command <<<
        echo "project_id = ~{query_project_id}" > ~/.bigqueryrc

        bq query --nouse_legacy_sql --project_id=~{query_project_id} --format=csv
        'SELECT
            vid,
            transcript,
            COUNT(vid) AS num_vids,
            COUNT(transcript) AS num_transcripts
        FROM
            ~{fq_vat_table}
        GROUP BY vid, transcript
        HAVING num_vids > 1 OR num_transcripts > 1' > bq_primary_key.csv

        # get number of lines in bq query output
        NUMRESULTS=$(awk 'END{print NR}' bq_primary_key.csv)

        echo "false" > ~{pf_file}
        # if the result of the query has any rows, that means not all key combinations (vid+transcript) are unique, so report those back in the output
        if [[ $NUMRESULTS = "0" ]]; then
          echo "The VAT table ~{fq_vat_table} has all unique key combinations (vid+transcript)" > ~{results_file}
          echo "true" > ~{pf_file}
        else
          echo "The VAT table ~{fq_vat_table} had repeating key combinations (vid+transcript): [csv output follows] " > ~{results_file}
          cat bq_primary_key.csv >> ~{results_file}
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
        Boolean pass = read_boolean(pf_file)
        String name = "SchemaPrimaryKey"
        String result = read_string(results_file)
    }
}

task SchemaEnsemblTranscripts {
    input {
        String query_project_id
        String fq_vat_table
        String last_modified_timestamp
    }
    # Every transcript_source is Ensembl or null

    String pf_file = "pf.txt"
    String results_file = "results.txt"

    command <<<
        echo "project_id = ~{query_project_id}" > ~/.bigqueryrc

        bq query --nouse_legacy_sql --project_id=~{query_project_id} --format=csv 'SELECT
            contig,
            position,
            vid,
            transcript,
            transcript_source
        FROM
            ~{fq_vat_table}
        WHERE
            transcript IS NOT NULL AND
            transcript_source != "Ensembl"' > bq_transcript_output.csv

        # get number of lines in bq query output
        NUMRESULTS=$(awk 'END{print NR}' bq_transcript_output.csv)

        echo "false" > ~{pf_file}
        # if the result of the query has any rows, that means there were unexpected transcripts (not from Ensembl), so report those back in the output
        if [[ $NUMRESULTS = "0" ]]; then
            echo "The VAT table ~{fq_vat_table} only has the expected Ensembl transcripts"  > ~{results_file}
            echo "true" > ~{pf_file}
        else
            echo "The VAT table ~{fq_vat_table} had unexpected transcripts (not from Ensembl): [csv output follows] " > ~{results_file}
            cat bq_transcript_output.csv >> ~{results_file}
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
        Boolean pass = read_boolean(pf_file)
        String name = "SchemaEnsemblTranscripts"
        String result = read_string(results_file)
    }
}

task SchemaNonzeroAcAn {
    input {
        String query_project_id
        String fq_vat_table
        String last_modified_timestamp
    }
    # No row has AC of zero or AN of zero.

    String pf_file = "pf.txt"
    String results_file = "results.txt"

    command <<<
        echo "project_id = ~{query_project_id}" > ~/.bigqueryrc

        bq query --nouse_legacy_sql --project_id=~{query_project_id} --format=csv 'SELECT
            contig,
            position,
            vid,
            gvs_all_ac,
            gvs_all_an
        FROM
            ~{fq_vat_table}
        WHERE
            gvs_all_ac IS NULL OR
            gvs_all_ac = 0 OR
            gvs_all_an IS NULL OR
            gvs_all_an = 0' > bq_ac_an_output.csv


        # get number of lines in bq query output
        NUMRESULTS=$(awk 'END{print NR}' bq_ac_an_output.csv)

        echo "false" > ~{pf_file}
        # if the result of the query has any rows, that means there were unexpected rows with either an AC of zero or AN of zero, so report those back in the output
        if [[ $NUMRESULTS = "0" ]]; then
            echo "The VAT table ~{fq_vat_table} only has no rows with AC of zero or AN of zero" > ~{results_file}
            echo "true" > ~{pf_file}
        else
            echo "The VAT table ~{fq_vat_table} had unexpected rows with AC of zero or AN of zero: [csv output follows] " > ~{results_file}
            cat bq_ac_an_output.csv >> ~{results_file}
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
        Boolean pass = read_boolean(pf_file)
        String name = "SchemaNonzeroAcAn"
        String result = read_string(results_file)
    }
}

task SchemaNullTranscriptsExist {
    input {
        String query_project_id
        String fq_vat_table
        String last_modified_timestamp
    }

    String pf_file = "pf.txt"
    String results_file = "results.txt"

    command <<<
        set -e

        echo "project_id = ~{query_project_id}" > ~/.bigqueryrc

        bq query --nouse_legacy_sql --project_id=~{query_project_id} --format=csv 'SELECT
            vid
        FROM
            ~{fq_vat_table}
        WHERE
            transcript_source is NULL AND
            transcript is NULL' > bq_variant_count.csv

        # get number of lines in bq query output
        NUMRESULTS=$(awk 'END{print NR}' bq_variant_count.csv)

        echo "false" > ~{pf_file}
        # if the result of the query has any rows, that means there were null transcripts as expected
        if [[ $NUMRESULTS != "0" ]]; then
           echo "The VAT table ~{fq_vat_table} has at least one null transcript" > ~{results_file}
           echo "true" > ~{pf_file}
        else
           echo "The VAT table ~{fq_vat_table} has no null transcripts" > ~{results_file}
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
        Boolean pass = read_boolean(pf_file)
        String name = "SchemaNullTranscriptsExist"
        String result = read_string(results_file)
    }
}

task SubpopulationMax {
    input {
        String query_project_id
        String fq_vat_table
        String last_modified_timestamp
    }
    # gvs_max_af is actually the max

    String pf_file = "pf.txt"
    String results_file = "results.txt"

    command <<<
        set -e

        echo "project_id = ~{query_project_id}" > ~/.bigqueryrc

        # gvs subpopulations:  [ "afr", "amr", "eas", "eur", "mid", "oth", "sas"]

        bq query --nouse_legacy_sql --project_id=~{query_project_id} --format=csv 'SELECT
            vid
        FROM
            ~{fq_vat_table}
        WHERE
            gvs_max_af < gvs_afr_af OR
            gvs_max_af < gvs_amr_af OR
            gvs_max_af < gvs_eas_af OR
            gvs_max_af < gvs_eur_af OR
            gvs_max_af < gvs_mid_af OR
            gvs_max_af < gvs_oth_af OR
            gvs_max_af < gvs_sas_af' > bq_query_output.csv

        # get number of lines in bq query output
        NUMRESULTS=$(awk 'END{print NR}' bq_query_output.csv)

        echo "false" > ~{pf_file}
        # if the result of the query has any rows, that means gvs_max_af is not in fact the max af
        if [[ $NUMRESULTS = "0" ]]; then
          echo "The VAT table ~{fq_vat_table} has a correct calculation for subpopulation" > ~{results_file}
          echo "true" > ~{pf_file}
        else
          echo "The VAT table ~{fq_vat_table} has an incorrect calculation for subpopulation" > ~{results_file}
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
        Boolean pass = read_boolean(pf_file)
        String name = "SubpopulationMax"
        String result = read_string(results_file)
    }
}

task SubpopulationAlleleCount {
    input {
        String query_project_id
        String fq_vat_table
        String last_modified_timestamp
    }
    # sum of subpop ACs equal the gvs_all ACs

    String pf_file = "pf.txt"
    String results_file = "results.txt"

    command <<<
        set -e

        echo "project_id = ~{query_project_id}" > ~/.bigqueryrc

        # gvs subpopulations:  [ "afr", "amr", "eas", "eur", "mid", "oth", "sas"]

        bq query --nouse_legacy_sql --project_id=~{query_project_id} --format=csv 'SELECT
            vid
        FROM
            ~{fq_vat_table}
        WHERE
            gvs_all_ac != gvs_afr_ac + gvs_amr_ac + gvs_eas_ac + gvs_eur_ac + gvs_mid_ac + gvs_oth_ac + gvs_sas_ac'  > bq_query_output.csv

        # get number of lines in bq query output
        NUMRESULTS=$(awk 'END{print NR}' bq_query_output.csv)

        echo "false" > ~{pf_file}
        # if the result of the query has any rows, that means gvs_all_ac has not been calculated correctly
        if [[ $NUMRESULTS = "0" ]]; then
            echo "The VAT table ~{fq_vat_table} has a correct calculation for AC and the AC of subpopulations" > ~{results_file}
            echo "true" > ~{pf_file}
        else
            echo "The VAT table ~{fq_vat_table} has an incorrect calculation for AC and the AC of subpopulations" > ~{results_file}
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
        Boolean pass = read_boolean(pf_file)
        String name = "SubpopulationAlleleCount"
        String result = read_string(results_file)
    }
}

task SubpopulationAlleleNumber {
    input {
        String query_project_id
        String fq_vat_table
        String last_modified_timestamp
    }
    # sum of subpop ACs equal the gvs_all ACs

    String pf_file = "pf.txt"
    String results_file = "results.txt"

    command <<<
        set -e

        echo "project_id = ~{query_project_id}" > ~/.bigqueryrc

        # gvs subpopulations:  [ "afr", "amr", "eas", "eur", "mid", "oth", "sas"]

        bq query --nouse_legacy_sql --project_id=~{query_project_id} --format=csv 'SELECT
        vid
        FROM
        ~{fq_vat_table}
        WHERE
        gvs_all_an != gvs_afr_an + gvs_amr_an + gvs_eas_an + gvs_eur_an + gvs_mid_an + gvs_oth_an + gvs_sas_an' > bq_an_output.csv

        # get number of lines in bq query output
        NUMRESULTS=$(awk 'END{print NR}' bq_an_output.csv)

        echo "false" > ~{pf_file}
        # if the result of the query has any rows, that means gvs_all_an has not been calculated correctly
        if [[ $NUMRESULTS = "0" ]]; then
          echo "The VAT table ~{fq_vat_table} has a correct calculation for AN and the AN of subpopulations" > ~{results_file}
          echo "true" > ~{pf_file}
        else
          echo "The VAT table ~{fq_vat_table} has an incorrect calculation for AN and the AN of subpopulations" > ~{results_file}
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
        Boolean pass = read_boolean(pf_file)
        String name = "SubpopulationAlleleNumber"
        String result = read_string(results_file)
    }
}

task ClinvarSignificance {
    input {
        String query_project_id
        String fq_vat_table
        String last_modified_timestamp
    }
    # check that all clinvar values are accounted for

    String pf_file = "pf.txt"
    String results_file = "results.txt"

    command <<<
        set -e

        echo "project_id = ~{query_project_id}" > ~/.bigqueryrc

        # clinvar significance values:  ["benign",
        #                                 "likely benign",
        #                                 "uncertain significance",
        #                                 "likely pathogenic",
        #                                 "pathogenic",
        #                                 "drug response",
        #                                 "association",
        #                                 "risk factor",
        #                                 "protective",
        #                                 "affects",
        #                                 "conflicting data from submitters",
        #                                 "other",
        #                                 "not provided"]

        bq query --nouse_legacy_sql --project_id=~{query_project_id} --format=csv 'SELECT
          distinct(unnested_clinvar_classification)
          FROM
        ~{fq_vat_table}, UNNEST(clinvar_classification) AS unnested_clinvar_classification' > bq_clinvar_classes.csv

        INCLUVALUES=$(awk -v RS='^$' 'END{print  !(index($0,"benign") && \
         index($0,"likely benign") && index($0,"uncertain significance") && \
         index($0,"likely pathogenic") && index($0,"pathogenic") && \
         index($0,"drug response") && index($0,"association") && \
         index($0,"risk factor") && index($0,"protective") && \
         index($0,"affects") && index($0,"conflicting data from submitters") && \
         index($0,"other") && \
         index($0,"not provided"))}'  bq_clinvar_classes.csv)

        NUMRESULTS=$( wc -l bq_clinvar_classes.csv | awk '{print $1;}' ) # we expect this to be 13+

        echo "false" > ~{pf_file}
        # if the result of the query has any rows, that means gvs_all_an has not been calculated correctly
        if [[ $NUMRESULTS -ge 13 && $INCLUVALUES = "0" ]]; then
          echo "The VAT table ~{fq_vat_table} has the correct values for clinvar classification" > ~{results_file}
          echo "true" > ~{pf_file}
        else
          echo "The VAT table ~{fq_vat_table} has incorrect values for clinvar classification" > ~{results_file}
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
        Boolean pass = read_boolean(pf_file)
        String name = "ClinvarSignificance"
        String result = read_string(results_file)
    }
}

task SchemaAAChangeAndExonNumberConsistent {
    input {
        String query_project_id
        String fq_vat_table
        String last_modified_timestamp
    }
    # Check that cases where (aa_change non-null AND exon_number null) OR (aa_change null AND exon_number non-null)
    # are all accounted for.
    # There are edge cases where aa_change is non-null AND exon_number is null AND intron_variant has a conseqence.
    # These cases are caught in this test, but too broadly, so we have an additional test:
    # SpotCheckForAAChangeAndExonNumberConsistency whose purpose is to spot check some of these edge cases AND normal cases
    # to validate our understanding (and that this doesn't change in future updates).

    String pf_file = "pf.txt"
    String results_file = "results.txt"

    command <<<
        set -e

        echo "project_id = ~{query_project_id}" > ~/.bigqueryrc

        bq query --nouse_legacy_sql --project_id=~{query_project_id} --format=csv 'SELECT
        COUNT (DISTINCT vid) AS count FROM
        (
            SELECT vid
            FROM ~{fq_vat_table}
            WHERE exon_number IS NULL
            AND aa_change IS NOT NULL
            AND NOT (("intron_variant" IN UNNEST(consequence) AND (variant_type IN ("insertion", "deletion"))))
            AND "splice_region_variant" NOT IN UNNEST(consequence)
            AND "splice_acceptor_variant" NOT IN UNNEST(consequence)
            AND "splice_donor_variant" NOT IN UNNEST(consequence)
            AND "NMD_transcript_variant" NOT IN UNNEST(consequence)
            AND "coding_sequence_variant" NOT IN UNNEST(consequence)
        UNION DISTINCT
            SELECT vid
            FROM ~{fq_vat_table}
            WHERE exon_number IS NOT NULL
            AND aa_change IS NULL
            AND "non_coding_transcript_exon_variant" NOT IN UNNEST(consequence)
            AND "3_prime_UTR_variant" NOT IN UNNEST(consequence)
            AND "5_prime_UTR_variant" NOT IN UNNEST(consequence)
            AND "splice_donor_variant" NOT IN UNNEST(consequence)
            AND "splice_acceptor_variant" NOT IN UNNEST(consequence)
            AND "splice_region_variant" NOT IN UNNEST(consequence)
            AND ("stop_gained" NOT IN UNNEST(consequence) AND "frameshift_variant" NOT IN UNNEST(consequence))
            AND ("stop_lost" NOT IN UNNEST(consequence) AND "frameshift_variant" NOT IN UNNEST(consequence))
            AND ("stop_lost" NOT IN UNNEST(consequence) AND "inframe_deletion" NOT IN UNNEST(consequence))
            AND "NMD_transcript_variant" NOT IN UNNEST(consequence)
            AND "stop_retained_variant" NOT IN UNNEST(consequence)
            AND "incomplete_terminal_codon_variant" NOT IN UNNEST(consequence)
            AND "transcript_variant" NOT IN UNNEST(consequence)
            AND "coding_sequence_variant" NOT IN UNNEST(consequence)
            AND "inframe_insertion" NOT IN UNNEST(consequence)
            AND "inframe_deletion" NOT IN UNNEST(consequence)
            AND "frameshift_variant" NOT IN UNNEST(consequence)
            AND "mature_miRNA_variant" NOT IN UNNEST(consequence)
        )' > bq_aachange_exonnumber.csv

        NUMVARS=$(python3 -c "csvObj=open('bq_aachange_exonnumber.csv','r');csvContents=csvObj.read();print(csvContents.split('\n')[1]);")

        echo "false" > ~{pf_file}
        # if the result of the bq call and the csv parsing is a series of digits, then check that it isn't 0
        if [[ $NUMVARS =~ ^[0-9]+$ ]]; then
            if [[ $NUMVARS = "0" ]]; then
                echo "The VAT table ~{fq_vat_table} has consistency across all aa_change and exon_number values in it." > ~{results_file}
                echo "true" > ~{pf_file}
            else
                echo "The VAT table ~{fq_vat_table} has $NUMVARS variants for which aa_change and exon_number are inconsistent." > ~{results_file}
            fi
        # otherwise, something is off, so return the output from the bq query call
        else
            echo "Something went wrong. The attempt to count the variants returned: " $(cat bq_aachange_exonnumber.csv) >&2
            exit 1
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
        Boolean pass = read_boolean(pf_file)
        String name = "SchemaAAChangeAndExonNumberConsistent"
        String result = read_string(results_file)
    }
}

task SpotCheckForAAChangeAndExonNumberConsistency {
    input {
        String query_project_id
        String fq_vat_table
        String last_modified_timestamp
    }

    String pf_file = "pf.txt"
    String results_file = "results.txt"

    # This test runs a spot check on the VAT table to verify that cases where (aa_change non-null AND exon_number null)
    # are understood.
    # The first four cases are cases where aa_change is non-null and exon_number is null
    # There are (currenly) 857 cases where this is allowed in our testing - these cases have a consequence of
    # intron_variant AND they are insertions or deletions. The test task SchemaAAChangeAndExonNumberConsistent (above)
    # allows these to pass.
    # BUT, since this is a fairly broad case we also spot-check that the normal case is still seen.
    # That is the purpose of the other four cases.
    # These cases test that for the 'normal' case (both aa_change is null AND exon_number is null)
    # we find cases which have a consequence of intron_variant AND they are insertions or deletionss

    command <<<
        set -e

        echo "project_id = ~{query_project_id}" > ~/.bigqueryrc

        bq query --nouse_legacy_sql --project_id=~{query_project_id} --format=csv 'SELECT
        COUNT (DISTINCT vid) FROM
        (
            SELECT vid, transcript, exon_number, aa_change, consequence, variant_type
            FROM ~{fq_vat_table}
            WHERE vid = "12-42446055-C-CGACCCTGCA" AND transcript = "ENST00000610488.4" AND exon_number IS NULL
            AND aa_change = "ENSP00000479913.1:p.(Ser375_Ala376insThrLeuGln)" AND "intron_variant" IN UNNEST(consequence) AND variant_type = "insertion"
            -- Another insertion
        UNION ALL
            SELECT vid, transcript, exon_number, aa_change, consequence, variant_type
            FROM ~{fq_vat_table}
            WHERE vid = "16-29810099-G-GGCGGCGGCAGCGGCA" AND transcript = "ENST00000566906.6" AND exon_number IS NULL
            AND aa_change = "ENSP00000461174.1:p.(Gly95_Ser99dup)" AND "intron_variant" IN UNNEST(consequence) AND variant_type = "insertion"
        UNION ALL
            SELECT vid, transcript, exon_number, aa_change, consequence, variant_type
            FROM ~{fq_vat_table}
            WHERE vid = "7-114631524-CCAG-C" AND transcript = "ENST00000393495.7" AND exon_number IS NULL
            AND aa_change = "ENSP00000377133.3:p.(Gln63del)" AND "intron_variant" IN UNNEST(consequence) AND variant_type = "deletion"
            --another
        UNION ALL
            SELECT vid, transcript, exon_number, aa_change, consequence, variant_type
            FROM ~{fq_vat_table}
            WHERE vid = "7-23198205-CTTTTTTT-C" and transcript = "ENST00000548367.2" AND exon_number IS NULL
            AND aa_change = "ENSP00000482736.1:p.(Phe25Ter)" AND "intron_variant" IN UNNEST(consequence) AND variant_type = "deletion"
        UNION ALL
            SELECT vid, transcript, exon_number, aa_change, consequence, variant_type
            FROM ~{fq_vat_table}
            WHERE vid = "1-40060064-T-TA" AND transcript = "ENST00000449311.5"
            AND exon_number IS NULL AND aa_change IS NULL
            AND "intron_variant" IN UNNEST(consequence) AND variant_type = "insertion"
        UNION ALL
            SELECT vid, transcript, exon_number, aa_change, consequence, variant_type
            FROM ~{fq_vat_table}
            WHERE vid = "17-1483061-T-TG" AND transcript = "ENST00000574790.5"
            AND exon_number IS NULL AND aa_change IS NULL
            AND "intron_variant" IN UNNEST(consequence) AND variant_type = "insertion"
        UNION ALL
            SELECT vid, transcript, exon_number, aa_change, consequence, variant_type
            FROM ~{fq_vat_table}
            WHERE vid = "1-46560284-TG-T" AND transcript = "ENST00000371945.8"
            AND exon_number IS NULL AND aa_change IS NULL
            AND "intron_variant" IN UNNEST(consequence) AND variant_type = "deletion"
        UNION ALL
            SELECT vid, transcript, exon_number, aa_change, consequence, variant_type
            FROM ~{fq_vat_table}
            WHERE vid = "6-111741100-CA-C" AND transcript = "ENST00000518295.5"
            AND exon_number IS NULL AND aa_change IS NULL
            AND "intron_variant" IN UNNEST(consequence) AND variant_type = "deletion"
        )' > output.csv

        NUMVARS=$(python3 -c "csvObj=open('output.csv','r');csvContents=csvObj.read();print(csvContents.split('\n')[1]);")

        echo "false" > ~{pf_file}
        # if the result of the bq call and the csv parsing is a series of digits, then check that it isn't 0
        if [[ $NUMVARS =~ ^[0-9]+$ ]]; then
            if [[ $NUMVARS = "8" ]]; then
                echo "true" > ~{pf_file}
                echo "The VAT table ~{fq_vat_table} has been spot checked for aa_change and exon_number consistency." > ~{results_file}
            else
                echo "The VAT table ~{fq_vat_table} has failed the spot check for aa_change and exon_number consistency." > ~{results_file}
        fi
        # otherwise, something is off, so return the output from the bq query call
        else
            echo "Something went wrong. The attempt to count the spot checked entries returned: " $(cat output.csv.csv) >&2
            exit 1
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
        Boolean pass = read_boolean(pf_file)
        String name = "SpotCheckForAAChangeAndExonNumberConsistency"
        String result = read_string(results_file)
    }
}

task GenerateFinalReport {
    input {
        Array[Boolean] pf_flags
        Array[String] validation_names
        Array[String] validation_results
    }

    String report_file = "report.txt"

    command <<<
        set -e

        declare -a PFS=(~{sep=' ' pf_flags})
        declare -a VNS=(~{sep=' ' validation_names})

        # the strings in validation_results may contain spaces, so we cannot directly sep it (by space) into an array
        # So, first convert it into a string with constituent strings sep'd by the '&' (chosen as it is unlikely to be used in results)
        # And then convert the string into a bash array using read with IFS set to '&'
        declare -a VRS_STRING="~{sep='&' validation_results}"
        IFS='&' read -r -a VRS <<< "$VRS_STRING"

        exit_code=0

        echo "VAT Table Validation Results" > ~{report_file}
        echo "" >> ~{report_file}
        echo "P/F   Name    Results" >> ~{report_file}
        for i in ${!PFS[@]};
        do
            if [ "${PFS[$i]}" = "true" ]; then
                echo "PASS  ${VNS[$i]}  ${VRS[$i]}" >> ~{report_file}
            else
                echo "FAIL  ${VNS[$i]}  ${VRS[$i]}" >> ~{report_file}
                exit_code=1
            fi
        done

        if [ $exit_code -ne 0 ]; then
            # cat the report to stderr as this task is a fail
            cat ~{report_file} >&2
        fi

        exit $exit_code
    >>>

    output {
        File results = report_file
    }
    runtime {
        docker: "python:3.8-slim-buster"
        memory: "1 GB"
        disks: "local-disk 10 HDD"
        preemptible: 3
        cpu: 1
    }
}


## TODO It would be great to spot check a few well known variants / genes