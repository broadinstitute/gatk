version 1.0

import "GvsCreateVATAnnotations.wdl" as Annotations

workflow GvsCreateVAT {
   input {
        #String data_project
        #String default_dataset

        File inputFileofFileNames
        File inputFileofIndexFileNames
        String project_id
        String dataset_name
        File nirvana_data_directory
        File vat_schema_json_file
        File variant_transcript_schema_json_file
        File genes_schema_json_file
        String output_path
        String table_suffix

        String? service_account_json_path
        File? gatk_override
        File AnAcAf_annotations_template
        File ancestry_file
        File reference
    }

    Array[String] contig_array = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"]


    call MakeSubpopulationFiles {
        input:
            input_ancestry_file = ancestry_file,
            service_account_json_path = service_account_json_path,
            inputFileofFileNames = inputFileofFileNames,
            inputFileofIndexFileNames = inputFileofIndexFileNames
    }

    ## Scatter across the shards from the GVS jointVCF
    scatter(i in range(length(MakeSubpopulationFiles.input_vcfs)) ) {
        ## Create a sites-only VCF from the original GVS jointVCF
        ## Calculate AC/AN/AF for subpopulations and extract them for custom annotations
        ## To prevent premature failures from this brittle step, the output will be saved to GCP
        String input_vcf_name = basename(MakeSubpopulationFiles.input_vcfs[i], ".vcf.gz")
        call Annotations.GvsCreateVATAnnotations {
            input:
              input_vcf = MakeSubpopulationFiles.input_vcfs[i],
              input_vcf_index = MakeSubpopulationFiles.input_vcf_indices[i],
              ancestry_mapping_list = MakeSubpopulationFiles.ancestry_mapping_list,
              nirvana_data_directory = nirvana_data_directory,
              output_path = output_path,
              service_account_json_path = service_account_json_path,
              custom_annotations_template = AnAcAf_annotations_template,
              ref = reference
        }
    }


    ## TODO move everything to a "done" folder so that stuff doesn't get loaded in twice

## presumably I want all of these above to finish running? And to be stored in a GCP bucket? Then how do I get all of the files from there?
## technically there are (approx) the same number as the loop before except for the failures -- how am I going to want to deal with the failures?

    call BigQueryLoadJson {
       input:
         nirvana_schema = vat_schema_json_file,
         vt_schema = variant_transcript_schema_json_file,
         genes_schema = genes_schema_json_file,
         project_id = project_id,
         dataset_name = dataset_name,
         output_path = output_path,
         table_suffix = table_suffix,
         service_account_json_path = service_account_json_path,
         prep_jsons_done = PrepAnnotationJson.done
  }

    call BigQuerySmokeTest {
       input:
         project_id = project_id,
         dataset_name = dataset_name,
         counts_variants = ExtractAnAcAfFromVCF.count_variants,
         track_dropped_variants = ExtractAnAcAfFromVCF.track_dropped,
         table_suffix = table_suffix,
         service_account_json_path = service_account_json_path,
         load_jsons_done = BigQueryLoadJson.done
    }

    scatter(i in range(length(contig_array)) ) {
      call BigQueryExportVat {
        input:
            contig = contig_array[i],
            project_id = project_id,
            dataset_name = dataset_name,
            output_path = output_path,
            table_suffix = table_suffix,
            service_account_json_path = service_account_json_path,
            validate_jsons_done = BigQuerySmokeTest.done
      }
    }
}

################################################################################

task MakeSubpopulationFiles {
    input {
        File input_ancestry_file
        String? service_account_json_path
        File inputFileofFileNames
        File inputFileofIndexFileNames
    }
    parameter_meta {
        input_ancestry_file: {
          localization_optional: true
        }
        inputFileofFileNames: {
          localization_optional: true
        }
        inputFileofIndexFileNames: {
          localization_optional: true
        }
    }
    String output_ancestry_filename =  "ancestry_mapping.tsv"
    String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'
    String updated_input_ancestry_file = if (defined(service_account_json_path)) then basename(input_ancestry_file) else input_ancestry_file
    String updated_input_vcfs_file = if (defined(service_account_json_path)) then basename(inputFileofFileNames) else inputFileofFileNames
    String updated_input_indices_file = if (defined(service_account_json_path)) then basename(inputFileofIndexFileNames) else inputFileofIndexFileNames

    command <<<
        set -e

        if [ ~{has_service_account_file} = 'true' ]; then
          gsutil cp ~{service_account_json_path} local.service_account.json
          export GOOGLE_APPLICATION_CREDENTIALS=local.service_account.json
          gcloud auth activate-service-account --key-file=local.service_account.json

          gsutil cp ~{input_ancestry_file} .
          gsutil cp ~{inputFileofFileNames} .
          gsutil cp ~{inputFileofIndexFileNames} .
       fi

        ## the ancestry file is processed down to a simple mapping from sample to subpopulation
        python3 /app/extract_subpop.py \
          --input_path ~{updated_input_ancestry_file} \
          --output_path ~{output_ancestry_filename}

    >>>

    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:ah_var_store_20211101"
        memory: "1 GB"
        preemptible: 3
        cpu: "1"
        disks: "local-disk 100 HDD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        File ancestry_mapping_list = "~{output_ancestry_filename}"
        Array[File] input_vcfs = read_lines(updated_input_vcfs_file)
        Array[File] input_vcf_indices = read_lines(updated_input_indices_file)
    }
}






task BigQueryLoadJson {
    meta { # since the WDL will not see the updated data (its getting put in a gcp bucket)
        volatile: true
    }

    input {
        File nirvana_schema
        File vt_schema
        File genes_schema
        String project_id
        String dataset_name
        String output_path
        String table_suffix
        String? service_account_json_path
        Array[String] prep_jsons_done
    }

    # There are two pre-vat tables. A variant table and a genes table. They are joined together for the vat table

    String vat_table = "vat_" + table_suffix
    String variant_transcript_table = "vat_vt_"  + table_suffix
    String genes_table = "vat_genes_" + table_suffix

    String vt_path = output_path + 'vt/*'
    String genes_path = output_path + 'genes/*'

    String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

    command <<<

       echo "project_id = ~{project_id}" > ~/.bigqueryrc

       DATE=86400 ## 24 hours in seconds


       if [ ~{has_service_account_file} = 'true' ]; then
            gsutil cp ~{service_account_json_path} local.service_account.json
            export GOOGLE_APPLICATION_CREDENTIALS=local.service_account.json
            gcloud auth activate-service-account --key-file=local.service_account.json
            gcloud config set project ~{project_id}
       fi

       set +e
       bq show --project_id ~{project_id} ~{dataset_name}.~{variant_transcript_table} > /dev/null
       BQ_SHOW_RC=$?
       set -e

       if [ $BQ_SHOW_RC -ne 0 ]; then
         echo "Creating a pre-vat table ~{dataset_name}.~{variant_transcript_table}"
         bq --location=US mk --expiration=$DATE --project_id=~{project_id}  ~{dataset_name}.~{variant_transcript_table} ~{vt_schema}
       fi

       echo "Loading data into a pre-vat table ~{dataset_name}.~{variant_transcript_table}"
       echo ~{vt_path}
       echo ~{genes_path}
       bq --location=US load --project_id=~{project_id} --source_format=NEWLINE_DELIMITED_JSON ~{dataset_name}.~{variant_transcript_table} ~{vt_path}

       set +e
       bq show --project_id ~{project_id} ~{dataset_name}.~{genes_table} > /dev/null
       BQ_SHOW_RC=$?
       set -e

       if [ $BQ_SHOW_RC -ne 0 ]; then
         echo "Creating a pre-vat table ~{dataset_name}.~{genes_table}"
         bq --location=US mk --expiration=$DATE --project_id=~{project_id}  ~{dataset_name}.~{genes_table} ~{genes_schema}
       fi

       echo "Loading data into a pre-vat table ~{dataset_name}.~{genes_table}"
       bq --location=US load  --project_id=~{project_id} --source_format=NEWLINE_DELIMITED_JSON  ~{dataset_name}.~{genes_table} ~{genes_path}

       set +e
       bq show --project_id ~{project_id} ~{dataset_name}.~{vat_table} > /dev/null
       BQ_SHOW_RC=$?
       set -e

       if [ $BQ_SHOW_RC -ne 0 ]; then
         echo "Creating the vat table ~{dataset_name}.~{vat_table}"
         bq --location=US mk --project_id=~{project_id} ~{dataset_name}.~{vat_table} ~{nirvana_schema}
       else
         bq rm -t -f --project_id=~{project_id} ~{dataset_name}.~{vat_table}
         bq --location=US mk --project_id=~{project_id} ~{dataset_name}.~{vat_table} ~{nirvana_schema}
       fi
       echo "And putting data into it"

       # Now we run a giant query in BQ to get this all in the right table and join the genes properly
       # Note the genes table join includes the group by to avoid the duplicates that get created from genes that span shards
       # Commented out columns in the query are to be added in the next release
       # We want the vat creation query to overwrite the destination table because if new data has been put into the pre-vat tables
       # and this workflow has been run an additional time, we dont want duplicates being appended from the original run

       bq query --nouse_legacy_sql --destination_table=~{dataset_name}.~{vat_table} --replace --project_id=~{project_id} \
        'SELECT
              v.vid,
              v.transcript,
              v.contig,
              v.position,
              v.ref_allele,
              v.alt_allele,
              v.gvs_all_ac,
              v.gvs_all_an,
              v.gvs_all_af,
              v.gvs_all_sc,
              v.gvs_max_af,
              v.gvs_max_ac,
              v.gvs_max_an,
              v.gvs_max_sc,
              v.gvs_max_subpop,
              v.gvs_afr_ac,
              v.gvs_afr_an,
              v.gvs_afr_af,
              v.gvs_afr_sc,
              v.gvs_amr_ac,
              v.gvs_amr_an,
              v.gvs_amr_af,
              v.gvs_amr_sc,
              v.gvs_eas_ac,
              v.gvs_eas_an,
              v.gvs_eas_af,
              v.gvs_eas_sc,
              v.gvs_eur_ac,
              v.gvs_eur_an,
              v.gvs_eur_af,
              v.gvs_eur_sc,
              v.gvs_mid_ac,
              v.gvs_mid_an,
              v.gvs_mid_af,
              v.gvs_mid_sc,
              v.gvs_oth_ac,
              v.gvs_oth_an,
              v.gvs_oth_af,
              v.gvs_oth_sc,
              v.gvs_sas_ac,
              v.gvs_sas_an,
              v.gvs_sas_af,
              v.gvs_sas_sc,
              v.gene_symbol,
              v.transcript_source,
              v.aa_change,
              v.consequence,
              v.dna_change_in_transcript,
              v.variant_type,
              v.exon_number,
              v.intron_number,
              v.genomic_location,
              # v.hgvsc AS splice_distance
              v.dbsnp_rsid,
              v.gene_id,
              # v.entrez_gene_id,
              # g.hgnc_gene_id,
              g.gene_omim_id,
              CASE WHEN ( v.transcript is not null and v.is_canonical_transcript is not True)
                THEN False WHEN ( v.transcript is not null and v.is_canonical_transcript is True) THEN True END AS is_canonical_transcript,
              v.gnomad_all_af,
              v.gnomad_all_ac,
              v.gnomad_all_an,
              v.gnomad_failed_filter,
              v.gnomad_max_af,
              v.gnomad_max_ac,
              v.gnomad_max_an,
              v.gnomad_max_subpop,
              v.gnomad_afr_ac,
              v.gnomad_afr_an,
              v.gnomad_afr_af,
              v.gnomad_amr_ac,
              v.gnomad_amr_an,
              v.gnomad_amr_af,
              v.gnomad_asj_ac,
              v.gnomad_asj_an,
              v.gnomad_asj_af,
              v.gnomad_eas_ac,
              v.gnomad_eas_an,
              v.gnomad_eas_af,
              v.gnomad_fin_ac,
              v.gnomad_fin_an,
              v.gnomad_fin_af,
              v.gnomad_nfr_ac,
              v.gnomad_nfr_an,
              v.gnomad_nfr_af,
              v.gnomad_sas_ac,
              v.gnomad_sas_an,
              v.gnomad_sas_af,
              v.gnomad_oth_ac,
              v.gnomad_oth_an,
              v.gnomad_oth_af,
              v.revel,
              v.splice_ai_acceptor_gain_score,
              v.splice_ai_acceptor_gain_distance,
              v.splice_ai_acceptor_loss_score,
              v.splice_ai_acceptor_loss_distance,
              v.splice_ai_donor_gain_score,
              v.splice_ai_donor_gain_distance,
              v.splice_ai_donor_loss_score,
              v.splice_ai_donor_loss_distance,
              g.omim_phenotypes_id,
              g.omim_phenotypes_name,
              v.clinvar_classification,
              v.clinvar_last_updated,
              v.clinvar_phenotype,
              FROM `~{dataset_name}.~{variant_transcript_table}` as v
              left join
              (SELECT gene_symbol, ANY_VALUE(gene_omim_id) AS gene_omim_id, ANY_VALUE(omim_phenotypes_id) AS omim_phenotypes_id, ANY_VALUE(omim_phenotypes_name) AS omim_phenotypes_name FROM `~{dataset_name}.~{genes_table}` group by gene_symbol) as g
              on v.gene_symbol = g.gene_symbol'

  >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "openbridge/ob_google-bigquery:latest"
        memory: "3 GB"
        preemptible: 3
        cpu: "1"
        disks: "local-disk 100 HDD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        Boolean done = true
    }
}

task BigQuerySmokeTest {
    input {
        String project_id
        String dataset_name
        Array[Int] counts_variants
        Array[File] track_dropped_variants
        String table_suffix
        String? service_account_json_path
        Boolean load_jsons_done
    }
    # Now query the final table for expected results
    # Compare the number of variants we expect from the input with the size of the output / VAT
    # The number of passing variants in GVS matches the number of variants in the VAT.
    # Please note that we are counting the number of variants in GVS, not the number of sites, which may add a difficulty to this task.

    String vat_table = "vat_" + table_suffix

    String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

    command <<<
        set -e
        if [ ~{has_service_account_file} = 'true' ]; then
            gsutil cp ~{service_account_json_path} local.service_account.json
            gcloud auth activate-service-account --key-file=local.service_account.json
            gcloud config set project ~{project_id}
        fi
        echo "project_id = ~{project_id}" > ~/.bigqueryrc

        # ------------------------------------------------
        # VALIDATION CALCULATION
        # sum all the initial input variants across the shards

        INITIAL_VARIANT_COUNT=$(python -c "print(sum([~{sep=', ' counts_variants}]))")
        ## here a list of files: ~{sep=', ' track_dropped_variants}
        ## I need to get the lines from each file
        awk 'FNR==1{print ""}1' ~{sep=' ' track_dropped_variants} > dropped_variants.txt

        echo "Dropped variants:"
        cat dropped_variants.txt

        # Count number of variants in the VAT
        bq query --nouse_legacy_sql --project_id=~{project_id} --format=csv 'SELECT COUNT (DISTINCT vid) AS count FROM `~{dataset_name}.~{vat_table}`' > bq_variant_count.csv
        VAT_COUNT=$(python3 -c "csvObj=open('bq_variant_count.csv','r');csvContents=csvObj.read();print(csvContents.split('\n')[1]);")
        # if the result of the bq call and the csv parsing is a series of digits, then check that it matches the input
        if [[ $VAT_COUNT =~ ^[0-9]+$ ]]; then
            if [[ $INITIAL_VARIANT_COUNT -ne $VAT_COUNT ]]; then
                echo "FAIL: The VAT table ~{vat_table} has $VAT_COUNT variants in it, and the input files had $INITIAL_VARIANT_COUNT."
            else
                echo "PASS: The VAT table ~{vat_table} has $VAT_COUNT variants in it, which is the expected number."
            fi
        # otherwise, something is off, so return the output from the bq query call
        else
            echo "Something went wrong. The attempt to count the variants returned: " $(cat bq_variant_count.csv)
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
    # Outputs:
    output {
        Boolean done = true
    }
}

task BigQueryExportVat {
    input {
        String contig
        String project_id
        String dataset_name
        String output_path
        String table_suffix
        String? service_account_json_path
        Boolean validate_jsons_done
    }

    String vat_table = "vat_" + table_suffix
    String export_path = output_path + "export/" + contig + "/*.tsv.gz"

    String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

    command <<<

        echo "project_id = ~{project_id}" > ~/.bigqueryrc

        if [ ~{has_service_account_file} = 'true' ]; then
          gsutil cp ~{service_account_json_path} local.service_account.json
          export GOOGLE_APPLICATION_CREDENTIALS=local.service_account.json
          gcloud auth activate-service-account --key-file=local.service_account.json
          gcloud config set project ~{project_id}
        fi

        bq query --nouse_legacy_sql --project_id=~{project_id} \
        'EXPORT DATA OPTIONS(
        uri="~{export_path}",
        format="CSV",
        compression="GZIP",
        overwrite=true,
        header=true,
        field_delimiter=" ") AS
        SELECT
        vid,
        transcript,
        contig,
        position,
        ref_allele,
        alt_allele,
        gvs_all_ac,
        gvs_all_an,
        gvs_all_af,
        gvs_all_sc,
        gvs_max_af,
        gvs_max_ac,
        gvs_max_an,
        gvs_max_sc,
        gvs_max_subpop,
        gvs_afr_ac,
        gvs_afr_an,
        gvs_afr_af,
        gvs_afr_sc,
        gvs_amr_ac,
        gvs_amr_an,
        gvs_amr_af,
        gvs_amr_sc,
        gvs_eas_ac,
        gvs_eas_an,
        gvs_eas_af,
        gvs_eas_sc,
        gvs_eur_ac,
        gvs_eur_an,
        gvs_eur_af,
        gvs_eur_sc,
        gvs_mid_ac,
        gvs_mid_an,
        gvs_mid_af,
        gvs_mid_sc,
        gvs_oth_ac,
        gvs_oth_an,
        gvs_oth_af,
        gvs_oth_sc,
        gvs_sas_ac,
        gvs_sas_an,
        gvs_sas_af,
        gvs_sas_sc,
        gene_symbol,
        transcript_source,
        aa_change,
        (SELECT STRING_AGG(c, ", ") FROM UNNEST(ARRAY(SELECT x FROM UNNEST(consequence) AS x ORDER BY x)) as c) AS consequence,
        dna_change_in_transcript,
        variant_type,
        exon_number,
        intron_number,
        genomic_location,
        (SELECT STRING_AGG(d, ", ") FROM UNNEST(ARRAY(SELECT x FROM UNNEST(dbsnp_rsid) AS x ORDER BY x)) as d) AS dbsnp_rsid,
        gene_id,
        gene_omim_id,
        is_canonical_transcript,
        gnomad_all_af,
        gnomad_all_ac,
        gnomad_all_an,
        gnomad_failed_filter,
        gnomad_max_af,
        gnomad_max_ac,
        gnomad_max_an,
        gnomad_max_subpop,
        gnomad_afr_ac,
        gnomad_afr_an,
        gnomad_afr_af,
        gnomad_amr_ac,
        gnomad_amr_an,
        gnomad_amr_af,
        gnomad_asj_ac,
        gnomad_asj_an,
        gnomad_asj_af,
        gnomad_eas_ac,
        gnomad_eas_an,
        gnomad_eas_af,
        gnomad_fin_ac,
        gnomad_fin_an,
        gnomad_fin_af,
        gnomad_nfr_ac,
        gnomad_nfr_an,
        gnomad_nfr_af,
        gnomad_sas_ac,
        gnomad_sas_an,
        gnomad_sas_af,
        gnomad_oth_ac,
        gnomad_oth_an,
        gnomad_oth_af,
        revel,
        splice_ai_acceptor_gain_score,
        splice_ai_acceptor_gain_distance,
        splice_ai_acceptor_loss_score,
        splice_ai_acceptor_loss_distance,
        splice_ai_donor_gain_score,
        splice_ai_donor_gain_distance,
        splice_ai_donor_loss_score,
        splice_ai_donor_loss_distance,
        (SELECT STRING_AGG(CAST(id AS STRING), ", ") FROM UNNEST(omim_phenotypes_id) id) as omim_phenotypes_id,
        ARRAY_TO_STRING(omim_phenotypes_name, ", ") as omim_phenotypes_name,
        ARRAY_TO_STRING(clinvar_classification, ", ") as clinvar_classification,
        clinvar_last_updated,
        ARRAY_TO_STRING(clinvar_phenotype, ", ") as clinvar_phenotype,
        FROM `~{dataset_name}.~{vat_table}`
        WHERE contig="~{contig}"
        ORDER BY position
        '

    >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "openbridge/ob_google-bigquery:latest"
        memory: "2 GB"
        preemptible: 3
        cpu: "1"
        disks: "local-disk 100 HDD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        Boolean done = true
    }
}

