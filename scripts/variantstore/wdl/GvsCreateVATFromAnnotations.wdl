version 1.0

workflow GvsCreateVATFromAnnotations {
   input {
        File inputFileofFileNames
        String project_id
        String dataset_name
        File? vat_schema_json_file = "gs://broad-dsp-spec-ops/scratch/rcremer/Nirvana/schemas/vat_schema.json"
        File? variant_transcript_schema_json_file = "gs://broad-dsp-spec-ops/scratch/rcremer/Nirvana/schemas/vt_schema.json"
        File? genes_schema_json_file = "gs://broad-dsp-spec-ops/scratch/rcremer/Nirvana/schemas/genes_schema.json"
        String output_path
        String table_suffix
        String? service_account_json_path
    }

    Array[String] contig_array = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"]
    File reference = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"

    call GetAnnotations {
        input:
            service_account_json_path = service_account_json_path,
            inputFileofFileNames = inputFileofFileNames
    }

    ## Scatter across the shards from the GVS jointVCF
    scatter(i in range(length(GetAnnotations.input_jsons)) ) {
        ## Parse the annotations
        call PrepAnnotationJson {
          input:
            annotation_json = GetAnnotations.input_jsons[i],
            output_file_suffix = basename(GetAnnotations.input_jsons[i], "_annotated.json.gz") + ".json.gz",
            output_path = output_path,
            service_account_json_path = service_account_json_path
       }
    }

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

    scatter(i in range(length(contig_array)) ) {
      call BigQueryExportVat {
        input:
            contig = contig_array[i],
            project_id = project_id,
            dataset_name = dataset_name,
            output_path = output_path,
            table_suffix = table_suffix,
            service_account_json_path = service_account_json_path,
            load_jsons_done = BigQueryLoadJson.done
      }
    }
}

################################################################################

task GetAnnotations {
    input {
        String? service_account_json_path
        File inputFileofFileNames
    }
    parameter_meta {
        inputFileofFileNames: {
          localization_optional: true
        }
    }
    String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'
    String updated_input_files = if (defined(service_account_json_path)) then basename(inputFileofFileNames) else inputFileofFileNames

    command <<<
        set -e

        if [ ~{has_service_account_file} = 'true' ]; then
          gsutil cp ~{service_account_json_path} local.service_account.json
          export GOOGLE_APPLICATION_CREDENTIALS=local.service_account.json
          gcloud auth activate-service-account --key-file=local.service_account.json

          gsutil cp ~{inputFileofFileNames} ~{updated_input_files}
       fi

    >>>

    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:ah_var_store_2022_07_08"
        memory: "1 GB"
        preemptible: 3
        cpu: "1"
        disks: "local-disk 100 HDD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        Array[File] input_jsons = read_lines(updated_input_files)
    }
}

task PrepAnnotationJson {
    input {
        File annotation_json
        String output_file_suffix
        String output_path
        String? service_account_json_path
    }
    parameter_meta {
        annotation_json: {
          localization_optional: true
        }
    }
    String output_vt_json = "vat_vt_bq_load" + output_file_suffix
    String output_genes_json = "vat_genes_bq_load" + output_file_suffix
    String output_vt_gcp_path = output_path + 'vt/'
    String output_genes_gcp_path = output_path + 'genes/'

    String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'
    String updated_annotation_json = if (defined(service_account_json_path)) then basename(annotation_json) else annotation_json



    command <<<
        set -e

        if [ ~{has_service_account_file} = 'true' ]; then
            gsutil cp ~{service_account_json_path} local.service_account.json
            export GOOGLE_APPLICATION_CREDENTIALS=local.service_account.json
            gcloud auth activate-service-account --key-file=local.service_account.json

            gsutil cp ~{annotation_json} ~{updated_annotation_json}
        fi

        ## the annotation jsons are split into the specific VAT schema
        python3 /app/create_variant_annotation_table.py \
          --annotated_json ~{updated_annotation_json} \
          --output_vt_json ~{output_vt_json} \
          --output_genes_json ~{output_genes_json}

        gsutil cp ~{output_vt_json} '~{output_vt_gcp_path}'
        gsutil cp ~{output_genes_json} '~{output_genes_gcp_path}'

     >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:ah_var_store_2022_07_08"
        memory: "8 GB"
        preemptible: 5
        cpu: "1"
        disks: "local-disk 250 SSD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        File vat_vt_json="~{output_vt_json}"
        File vat_genes_json="~{output_genes_json}"
        Boolean done = true
    }
}

task BigQueryLoadJson {
    meta { # since the WDL will not see the updated data (its getting put in a gcp bucket)
        volatile: true
    }

    input {
        File? nirvana_schema
        File? vt_schema
        File? genes_schema
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

       if [ $BQ_SHOW_RC -eq 0 ]; then
         bq rm -t -f --project_id=~{project_id} ~{dataset_name}.~{vat_table}
       fi
       echo "Creating the vat table ~{dataset_name}.~{vat_table}"
       bq --location=US mk --project_id=~{project_id} ~{dataset_name}.~{vat_table} ~{nirvana_schema}

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

task BigQueryExportVat {
    input {
        String contig
        String project_id
        String dataset_name
        String output_path
        String table_suffix
        String? service_account_json_path
        String load_jsons_done
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

        # note: tab delimiter and compression creates tsv.gz files
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

