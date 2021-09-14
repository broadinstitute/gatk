version 1.0
workflow GvsSitesOnlyVCF {
   input {
        File inputFileofFileNames
        File inputFileofIndexFileNames
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
        File AnAcAf_annotations_template
        File ancestry_file
    }

    call MakeSubpopulationFiles {
        input:
            input_ancestry_file = ancestry_file,
            service_account_json_path = service_account_json_path,
            inputFileofFileNames = inputFileofFileNames,
            inputFileofIndexFileNames = inputFileofIndexFileNames
    }

    ## Scatter across the shards from the GVS jointVCF
    scatter(i in range(length(MakeSubpopulationFiles.input_vcfs)) ) {
        call ExtractAnAcAfFromVCF {
            input:
              input_vcf = MakeSubpopulationFiles.input_vcfs[i],
              input_vcf_index = MakeSubpopulationFiles.input_vcf_indices[i],
              service_account_json_path = service_account_json_path,
              subpopulation_sample_list = MakeSubpopulationFiles.ancestry_mapping_list,
              custom_annotations_template = AnAcAf_annotations_template
        }

      ## Create a sites-only VCF from the original GVS jointVCF
        call SitesOnlyVcf {
          input:
            input_vcf = ExtractAnAcAfFromVCF.output_vcf,
            input_vcf_index = ExtractAnAcAfFromVCF.output_vcf_index,
            service_account_json_path = service_account_json_path,
            output_filename = "${output_sites_only_file_name}_${i}.sites_only.vcf.gz"
        }

      ## Use Nirvana to annotate the sites-only VCF and include the AC/AN/AF calculations as custom annotations
        call AnnotateVCF {
          input:
            input_vcf = SitesOnlyVcf.output_vcf,
            input_vcf_index = SitesOnlyVcf.output_vcf_idx,
            output_annotated_file_name = "${output_annotated_file_name}_${i}",
            nirvana_data_tar = nirvana_data_directory,
            custom_annotations_file = ExtractAnAcAfFromVCF.annotations_file,
        }

        call PrepAnnotationJson {
          input:
            annotation_json = AnnotateVCF.annotation_json,
            output_file_suffix = "${i}.json.gz",
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

    call BigQuerySmokeTest {
       input:
         project_id = project_id,
         dataset_name = dataset_name,
         counts_variants = ExtractAnAcAfFromVCF.count_variants,
         table_suffix = table_suffix,
         service_account_json_path = service_account_json_path,
         load_jsons_done = BigQueryLoadJson.done
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


        python3 /app/extract_subpop.py \
          --input_path ~{updated_input_ancestry_file} \
          --output_path ~{output_ancestry_filename}

    >>>

    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:ah_var_store_20210910"
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


task ExtractAnAcAfFromVCF {
    input {
        File input_vcf
        File input_vcf_index
        String? service_account_json_path
        File subpopulation_sample_list
        File custom_annotations_template
    }
    parameter_meta {
        input_vcf: {
          localization_optional: true
        }
        input_vcf_index: {
          localization_optional: true
        }
    }

    String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'
    String custom_annotations_file_name = "an_ac_af.tsv"
    String local_input_vcf = basename(input_vcf)
    String local_input_vcf_index = basename(input_vcf_index)

    # separate multi-allelic sites into their own lines, remove deletions and extract the an/ac/af & sc
    command <<<
        set -e

        if [ ~{has_service_account_file} = 'true' ]; then
          gsutil cp ~{service_account_json_path} local.service_account.json
          export GOOGLE_APPLICATION_CREDENTIALS=local.service_account.json
          gcloud auth activate-service-account --key-file=local.service_account.json

        fi

        cp ~{custom_annotations_template} ~{custom_annotations_file_name}

        gsutil cp ~{input_vcf} ~{local_input_vcf}
        gsutil cp ~{input_vcf_index} ~{local_input_vcf_index}

        awk '{print $2}' ~{subpopulation_sample_list} | tail -n +2 | sort -u > collected_subpopulations.txt

        # expected_subpopulations = [
        # "afr",
        # "amr",
        # "eas",
        # "eur",
        # "mid",
        # "oth",
        # "sas"
        #]

        bcftools norm -m- ~{local_input_vcf} | bcftools plugin fill-tags -- -S ~{subpopulation_sample_list} \
        | bcftools query -f \
        '%CHROM\t%POS\t%REF\t%ALT\t%AC\t%AN\t%AF\t%AC_Hom\t%AC_Het\t%AC_afr\t%AN_afr\t%AF_afr\t%AC_Hom_afr\t%AC_Het_afr\t%AC_amr\t%AN_amr\t%AF_amr\t%AC_Hom_amr\t%AC_Het_amr\t%AC_eas\t%AN_eas\t%AF_eas\t%AC_Hom_eas\t%AC_Het_eas\t%AC_eur\t%AN_eur\t%AF_eur\t%AC_Hom_eur\t%AC_Het_eur\t%AC_mid\t%AN_mid\t%AF_mid\t%AC_Hom_mid\t%AC_Het_mid\t%AC_oth\t%AN_oth\t%AF_oth\t%AC_Hom_oth\t%AC_Het_oth\t%AC_sas\t%AN_sas\t%AF_sas\t%AC_Hom_sas\t%AC_Het_sas\n' \
        | grep -v "*" >> ~{custom_annotations_file_name}

        ### for validation of the pipeline
        tail -n +7 ~{custom_annotations_file_name} | grep -v "AC=0;" | grep "AC=" | grep "AN=" | grep "AF=" | grep -v "*" | wc -l > count.txt
        # I find this ^ clearer, but could also do a regex like:  grep "AC=[1-9][0-9]*;A[N|F]=[.0-9]*;A[N|F]=[.0-9]*"
        # Should this be where we do the filtering of the AC/AN/AF values rather than in the python?
    >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:ah_var_store_20210910"
        memory: "12 GB"
        preemptible: 3
        cpu: "2"
        disks: "local-disk 500 SSD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        File annotations_file = "~{custom_annotations_file_name}"
        Int count_variants = read_int("count.txt")
        File output_vcf = local_input_vcf
        File output_vcf_index = local_input_vcf_index
    }
}

task SitesOnlyVcf {
    input {
        File input_vcf
        File input_vcf_index
        String? service_account_json_path
        String output_filename
    }
    String output_vcf_idx = basename(output_filename) + ".tbi"

    command <<<
        set -e

        # Adding `--add-output-vcf-command-line false` so that the VCF header doesn't have a timestamp
        # in it so that downstream steps can call cache

        gatk --java-options "-Xmx12288m" \
            SelectVariants \
                -V ~{input_vcf} \
                --add-output-vcf-command-line false \
                --exclude-filtered \
                --sites-only-vcf-output \
                -O ~{output_filename}

     >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "broadinstitute/gatk:4.2.0.0"
        memory: "15 GB"
        preemptible: 3
        cpu: "2"
        disks: "local-disk 250 HDD"
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

    command <<<
        set -e

        # =======================================
        # Handle our data sources:

        # Extract the tar.gz:
        echo "Extracting annotation data sources tar/gzip file..."
        mkdir datasources_dir
        tar zxvf ~{nirvana_data_tar} -C datasources_dir  --strip-components 2
        DATA_SOURCES_FOLDER="$PWD/datasources_dir"


        # =======================================
        # Create custom annotations:
        echo "Creating custom annotations"
        mkdir customannotations_dir
        CUSTOM_ANNOTATIONS_FOLDER="$PWD/customannotations_dir"

        # Add AC/AN/AF as custom annotations
        dotnet ~{custom_creation_location} customvar\
             -r $DATA_SOURCES_FOLDER~{path_reference} \
             -i ~{custom_annotations_file} \
             -o $CUSTOM_ANNOTATIONS_FOLDER

        # =======================================
        # Create Nirvana annotations:


        dotnet ~{nirvana_location} \
             -c $DATA_SOURCES_FOLDER~{path} \
             --sd $DATA_SOURCES_FOLDER~{path_supplementary_annotations} \
             --sd $CUSTOM_ANNOTATIONS_FOLDER \
             -r $DATA_SOURCES_FOLDER~{path_reference} \
             -i ~{input_vcf} \
             -o ~{output_annotated_file_name}

    >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "annotation/nirvana:3.14"
        memory: "20 GB"
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
    String output_annotations_gcp_path = output_path + 'annotations/'

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

        # for debugging only
        gsutil cp ~{annotation_json} '~{output_annotations_gcp_path}'

     >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:ah_var_store_20210910"
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
