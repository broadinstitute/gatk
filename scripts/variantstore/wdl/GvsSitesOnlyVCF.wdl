version 1.0
workflow GvsSitesOnlyVCF {
   input {
        Array[File] gvs_extract_cohort_filtered_vcfs
        Array[File] gvs_extract_cohort_filtered_vcf_indices
        String output_sites_only_file_name
        String output_annotated_file_name
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
    }

    scatter(i in range(length(gvs_extract_cohort_filtered_vcfs)) ) {
        call SitesOnlyVcf {
          input:
            input_vcf = gvs_extract_cohort_filtered_vcfs[i],
            input_vcf_index = gvs_extract_cohort_filtered_vcf_indices[i],
            service_account_json_path = service_account_json_path,
            output_filename = "${output_sites_only_file_name}_${i}.sites_only.vcf.gz",
        }

        call AnnotateVCF {
          input:
            input_vcf = SitesOnlyVcf.output_vcf,
            input_vcf_index = SitesOnlyVcf.output_vcf_idx,
            output_annotated_file_name = "${output_annotated_file_name}_${i}",
            nirvana_data_tar = nirvana_data_directory,
            custom_annotations_file = custom_annotations_file
        }

       call PrepAnnotationJson {
         input:
           annotation_json = AnnotateVCF.annotation_json,
           output_file_suffix = "${i}.json.gz",
           output_path = output_path,
           service_account_json_path = service_account_json_path,
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

     call BigQuerySmokeTest {
         input:
             project_id = project_id,
             dataset_name = dataset_name,
             table_suffix = table_suffix,
             service_account_json_path = service_account_json_path,
             annotation_jsons = AnnotateVCF.annotation_json,
             load_jsons_done = BigQueryLoadJson.done
         }
}

################################################################################
task SitesOnlyVcf {
    input {
        File input_vcf
        File input_vcf_index
        String? service_account_json_path
        String output_filename
    }
    String output_vcf_idx = basename(output_filename) + ".tbi"

    String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'
    String input_vcf_basename = basename(input_vcf)
    String updated_input_vcf = if (defined(service_account_json_path)) then input_vcf_basename else input_vcf

    parameter_meta {
        input_vcf: {
            localization_optional: true
        }
        input_vcf_index: {
            localization_optional: true
        }
    }
    command <<<
        set -e

        if [ ~{has_service_account_file} = 'true' ]; then
            gsutil cp ~{service_account_json_path} local.service_account.json
            export GOOGLE_APPLICATION_CREDENTIALS=local.service_account.json
            gcloud auth activate-service-account --key-file=local.service_account.json

        gsutil cp ~{input_vcf} .
            gsutil cp ~{input_vcf_index} .
        fi

        # Adding `--add-output-vcf-command-line false` so that the VCF header doesn't have a timestamp
        # in it so that downstream steps can call cache

        gatk --java-options "-Xmx2048m" \
            SelectVariants \
                -V ~{updated_input_vcf} \
                --add-output-vcf-command-line false \
                --exclude-filtered \
                --sites-only-vcf-output \
                -O ~{output_filename}
     >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "broadinstitute/gatk:4.2.0.0"
        memory: "3 GB"
        preemptible: 3
        cpu: "1"
        disks: "local-disk 100 HDD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        File output_vcf="~{output_filename}"
        File output_vcf_idx="~{output_vcf_idx}"
    }
}

task ExtractACANAF {
    input {
        File vcf_bgz_gts
        File vcf_index
        String output_filename
    }
    String output_vcf_idx = basename(output_filename) + ".tbi" # or will this be .idx if from .vcf.gz? or ".tbi" if a .vcf
    command <<<

        set -e
        gatk --java-options "-Xmx2048m" \
            SelectVariants \
                -V ~{vcf_bgz_gts} \
                --exclude-filtered \
                --sites-only-vcf-output \
                -O ~{output_filename}
     >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "broadinstitute/gatk:4.2.0.0"
        memory: "3 GB"
        cpu: "1"
        preemptible: 3
        disks: "local-disk 100 HDD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        File output_vcf="~{output_filename}"
        File output_vcf_idx="~{output_vcf_idx}"
    }
}

task AnnotateVCF {
    input {
        File input_vcf
        File input_vcf_index
        String output_annotated_file_name
        File nirvana_data_tar
        File custom_annotations_file
    }
    String annotation_json_name = output_annotated_file_name + ".json.gz"
    String annotation_json_name_jsi = annotation_json_name + ".jsi"
    String nirvana_location = "/opt/nirvana/Nirvana.dll"
    String custom_creation_location = "/opt/nirvana/SAUtils.dll"
    String path = "/Cache/GRCh38/Both"
    String path_supplementary_annotations = "/SupplementaryAnnotation/GRCh38"
    String path_reference = "/References/Homo_sapiens.GRCh38.Nirvana.dat"
    String custom_annotations_dir = "/CA"
    command <<<
        set -e
        # =======================================
        # Create custom annotations:
        echo "Creating custom annotations"
        dotnet ~{custom_creation_location} customvar\
        -r $DATA_SOURCES_FOLDER~{path_reference} \
        -i ~{custom_annotations_file} \
        -o ~{custom_annotations_dir}
        # =======================================
        # Handle our data sources:

        # Extract the tar.gz:
        echo "Extracting annotation data sources tar/gzip file..."
        mkdir datasources_dir
        tar zxvf ~{nirvana_data_tar} -C datasources_dir  --strip-components 2
        DATA_SOURCES_FOLDER="$PWD/datasources_dir"


        dotnet ~{nirvana_location} \
             -c $DATA_SOURCES_FOLDER~{path} \
             --sd $DATA_SOURCES_FOLDER~{path_supplementary_annotations} \
             -r $DATA_SOURCES_FOLDER~{path_reference} \
             -i ~{input_vcf} \
             -o ~{output_annotated_file_name}


    >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "annotation/nirvana:3.14"
        memory: "5 GB"
        cpu: "2"
        preemptible: 5
        disks: "local-disk 250 SSD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        File annotation_json = "~{annotation_json_name}"
        File annotation_json_jsi = "~{annotation_json_name_jsi}"
    }
}

task PrepAnnotationJson {
    input {
        File annotation_json
        String output_file_suffix
        String output_path
        String? service_account_json_path
    }

    String output_vt_json = "vat_vt_bq_load" + output_file_suffix
    String output_genes_json = "vat_genes_bq_load" + output_file_suffix
    String output_vt_gcp_path = output_path + 'vt/'
    String output_genes_gcp_path = output_path + 'genes/'

    String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

    command <<<
        set -e

        python3 /app/create_variant_annotation_table.py \
          --annotated_json ~{annotation_json} \
          --output_vt_json ~{output_vt_json} \
          --output_genes_json ~{output_genes_json}

        if [ ~{has_service_account_file} = 'true' ]; then
            gsutil cp ~{service_account_json_path} local.service_account.json
            export GOOGLE_APPLICATION_CREDENTIALS=local.service_account.json
            gcloud auth activate-service-account --key-file=local.service_account.json
        fi

        gsutil cp ~{output_vt_json} '~{output_vt_gcp_path}'
        gsutil cp ~{output_genes_json} '~{output_genes_gcp_path}'

     >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:ah_var_store_20200709"
        memory: "3 GB"
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
         bq --location=US mk --project_id=~{project_id}  ~{dataset_name}.~{variant_transcript_table} ~{vt_schema}
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
         bq --location=US mk --project_id=~{project_id}  ~{dataset_name}.~{genes_table} ~{genes_schema}
       fi

       echo "Loading data into a pre-vat table ~{dataset_name}.~{genes_table}"
       bq --location=US load  --project_id=~{project_id} --source_format=NEWLINE_DELIMITED_JSON  ~{dataset_name}.~{genes_table} ~{genes_path}

       # create the final vat table with the correct fields
       PARTITION_FIELD="position"
       CLUSTERING_FIELD="vid"
       PARTITION_STRING="" #--range_partitioning=$PARTITION_FIELD,0,4000,4000"
       CLUSTERING_STRING="" #--clustering_fields=$CLUSTERING_FIELD"

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
              null AS gvs_all_ac,
              null AS gvs_all_an,
              null AS  gvs_all_af,
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
              # v.gnomad_max_af,
              # v.gnomad_max_ac,
              # v.gnomad_max_an,
              # null AS gnomad_max_subpop,
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

task BigQuerySmokeTest { # TO BE BROKEN UP AND POTENTIALLY RUN SEPARATELY
    input {
        String project_id
        String dataset_name
        Array[File] annotation_jsons
        Boolean load_jsons_done
        String? service_account_json_path
        String table_suffix
    }

    # Now query the final table for expected results

    String vat_table = "vat_" + table_suffix

    String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

    command <<<
        set +e

        if [ ~{has_service_account_file} = 'true' ]; then
            gsutil cp ~{service_account_json_path} local.service_account.json
            export GOOGLE_APPLICATION_CREDENTIALS=local.service_account.json
            gcloud auth activate-service-account --key-file=local.service_account.json
            gcloud config set project ~{project_id}
        fi

        # For the validation rules that use count in SQL--the results are parsed to grab all digits in the response, which will also include error responses

        # ------------------------------------------------
        # VALIDATION #1
        echo  "VALIDATION #1"
        # The pipeline completed without an error message
        # The VAT has data in it

        bq query --nouse_legacy_sql --project_id=~{project_id} 'SELECT COUNT (DISTINCT vid) FROM `~{dataset_name}.~{vat_table}`'
        BQ_VAT_VARIANT_COUNT=$(bq query --nouse_legacy_sql --project_id=~{project_id} 'SELECT COUNT (DISTINCT vid) AS count FROM `~{dataset_name}.~{vat_table}`'| tr -dc '0-9')
        echo $BQ_VAT_VARIANT_COUNT

        if [[ $BQ_VAT_VARIANT_COUNT -le 0 ]]
        then
          echo "The VAT has no data in it"
          echo  "Validation has failed"
        else
          echo "The VAT has been updated"
          echo "Validation #1 has passed"
        fi

        # ------------------------------------------------
        # VALIDATION #2
        # The number of passing variants in GVS matches the number of variants in the VAT.
        echo  "VALIDATION #2"
        # Please note that we are counting the number of variants in GVS, not the number of sites, which may add a difficulty to this task.
        ANNOTATE_JSON_VARIANT_COUNT=$(gunzip -dc ~{sep=" " annotation_jsons} | grep -o -i '"vid":' | wc -l)

        echo $ANNOTATE_JSON_VARIANT_COUNT

        if [[ $ANNOTATE_JSON_VARIANT_COUNT -ne $BQ_VAT_VARIANT_COUNT ]]
        then
          echo "The number of variants is incorrect"
          echo  "Validation has failed"
        else
          echo "The number of passing variants in GVS matches the number of variants in the VAT"
          echo "Validation #2 has passed"
        fi

        # ------------------------------------------------
        # VALIDATION #3
        # All variants in the TESK2 gene region (chr1:45,343,883-45,491,163) list multiple genes and those genes are always TESK2 and AL451136.1.
        echo  "VALIDATION #3"

        # Note: Lee has swapped this valdation to chr19
        # SELECT gene_symbol, consequence  FROM `~{dataset_name}.~{vat_table}` WHERE contig = "chr19" and position >= 35740407 and position <= 35740469 and gene_symbol!="IGFLR1" and gene_symbol!="AD000671.2"

        POSITIONAL_COUNT=$(bq query --nouse_legacy_sql --project_id=~{project_id} 'SELECT COUNT (DISTINCT vid) AS distinct_vid_count FROM `~{dataset_name}.~{vat_table}` WHERE contig = "chr1" and position >= 45343883 and position <= 45491163'| tr -dc '0-9')
        TESK2_COUNT=$(bq query --nouse_legacy_sql --project_id=~{project_id} 'SELECT COUNT (DISTINCT vid) AS distinct_vid_count FROM `~{dataset_name}.~{vat_table}` WHERE contig = "chr1" and position >= 45343883 and position <= 45491163 and gene_symbol="TESK2"'| tr -dc '0-9')
        AL451136_COUNT=$(bq query --nouse_legacy_sql --project_id=~{project_id} 'SELECT COUNT (DISTINCT vid) AS distinct_vid_count FROM `~{dataset_name}.~{vat_table}` WHERE contig = "chr1" and position >= 45343883 and position <= 45491163 and gene_symbol="AL451136.1"'| tr -dc '0-9')
        GENE_COUNT_SUM=$(( $TESK2_COUNT + $AL451136_COUNT ))

        if [[ $POSITIONAL_COUNT -ne $GENE_COUNT_SUM ]]
          then echo "There are unexpected genes in the TESK2 region"
          echo  "Validation has failed"
        else
          echo "All variants in the TESK2 gene region list multiple genes and those genes are always TESK2 and AL451136.1"
          echo "Validation #4 has passed"
        fi

        # ------------------------------------------------
        # VALIDATION #4
        # If a vid has a null transcript, then the vid is only in one row of the VAT.
        echo  "VALIDATION #4"
        # Get a count of all the rows in the vat with no transcript
        # Get a count of all distinct VID in vat with no transcript
        # Make sure they are the same
        BQ_VAT_ROWS_NO_TRANSCRIPT_COUNT=$(bq query --nouse_legacy_sql --project_id=~{project_id} 'SELECT COUNT (DISTINCT vid) AS distinct_vid_count FROM `~{dataset_name}.~{vat_table}` WHERE transcript IS NULL'| tr -dc '0-9')
        BQ_VAT_VARIANT_NO_TRANSCRIPT_COUNT=$(bq query --nouse_legacy_sql --project_id=~{project_id} 'SELECT COUNT(*) AS count FROM `~{dataset_name}.~{vat_table}` WHERE transcript IS NULL'| tr -dc '0-9')

        if [[ $BQ_VAT_ROWS_NO_TRANSCRIPT_COUNT -ne $BQ_VAT_VARIANT_NO_TRANSCRIPT_COUNT ]]
        then
          echo "The number of rows for variants with no transcripts is incorrect"
          echo  "Validation has failed"
        else
          echo "If a vid has a null transcript, then the vid is only in one row of the VAT"
          echo "Validation #4 has passed"
        fi

        # ------------------------------------------------
        # VALIDATION #5
        # There is a non-zero number of transcript fields with null values in the VAT.
        echo  "VALIDATION #5"

        if [[ $BQ_VAT_VARIANT_NO_TRANSCRIPT_COUNT -gt 0 ]]
        then
          echo  "Validation has failed"
        else
          echo "Validation #5 has passed"
        fi

        # ------------------------------------------------
        # VALIDATION #6
        # No non-nullable fields contain null values.
        # non-nullable fields: vid, contig, position, ref_allele, alt_allele, gvs_all_ac, gvs_all_an, gvs_all_af, variant_type, genomic_location
        echo  "VALIDATION #6"

        BQ_VAT_NULL=$(bq query --nouse_legacy_sql --project_id=~{project_id} 'SELECT COUNT (DISTINCT vid) as count FROM `~{dataset_name}.~{vat_table}` WHERE vid IS NULL OR contig IS NULL OR position IS NULL OR ref_allele IS NULL OR alt_allele IS NULL OR variant_type IS NULL OR genomic_location IS NULL')
        echo $BQ_VAT_NULL

        if [[ $BQ_VAT_NULL -gt 0 ]]
        then
          echo "A required value is null"
          echo  "Validation has failed"
        else
          echo "No required values are null"
          echo "Validation #6 has passed"
        fi

        # ------------------------------------------------
        # VALIDATION #7
        # Each key combination (vid+transcript) is unique.
        echo  "VALIDATION #7"
        # get the sum of all the distinct vids where transcript is null and all the distinct transcript where transcript is not null
        BQ_VAT_UNIQUE_IDS_COUNT=$(bq query --nouse_legacy_sql --project_id=~{project_id} 'SELECT COUNT(*) AS count FROM (SELECT vid, transcript FROM `~{dataset_name}.~{vat_table}` group by vid, transcript)')
        BQ_VAT_ROW_COUNT=$(bq query --nouse_legacy_sql --project_id=~{project_id} 'SELECT COUNT(*) AS count FROM `~{dataset_name}.~{vat_table}` ')

        if [ $BQ_VAT_UNIQUE_IDS_COUNT -ne $BQ_VAT_ROW_COUNT ];
        then echo "There are duplicate variant - transcript rows"
          echo  "Validation has failed"
        else:
          echo "Each key combination is unique"
          echo "Validation #7 has passed"
        fi

        # ------------------------------------------------

        # VALIDATION #8
        # If a vid has any non-null transcripts then one transcript must be Ensembl (transcript_source). Every transcript_source is Ensembl or null.
        echo  "VALIDATION #8"
        BQ_VAT_ENSEMBL_COUNT=$(bq query --nouse_legacy_sql --project_id=~{project_id} 'SELECT COUNT(*) AS count FROM `~{dataset_name}.~{vat_table}` where transcript is not null and transcript_source="Ensembl"')
        BQ_VAT_TRANSCRIPT_COUNT=$(bq query --nouse_legacy_sql --project_id=~{project_id} 'SELECT COUNT(*) AS count FROM `~{dataset_name}.~{vat_table}` where transcript is not null')

        if [[ $BQ_VAT_ENSEMBL_COUNT -ne $BQ_VAT_TRANSCRIPT_COUNT ]]
        then
          echo "All transcripts should be from Ensembl"
          echo  "Validation has failed"
        else
          echo "If a vid has any non-null transcripts then one transcript must be Ensembl"
          echo "Validation #8 has passed"
        fi

    >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "openbridge/ob_google-bigquery:latest"
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
