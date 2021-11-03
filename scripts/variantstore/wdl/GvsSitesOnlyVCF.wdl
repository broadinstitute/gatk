version 1.0
workflow GvsSitesOnlyVCF {
   input {
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
        String input_vcf_name = basename(MakeSubpopulationFiles.input_vcfs[i], ".vcf.gz")
        call ExtractAnAcAfFromVCF {
            input:
              input_vcf = MakeSubpopulationFiles.input_vcfs[i],
              input_vcf_index = MakeSubpopulationFiles.input_vcf_indices[i],
              service_account_json_path = service_account_json_path,
              subpopulation_sample_list = MakeSubpopulationFiles.ancestry_mapping_list,
              custom_annotations_template = AnAcAf_annotations_template,
              ref = reference
        }

      ## Use Nirvana to annotate the sites-only VCF and include the AC/AN/AF calculations as custom annotations
        call AnnotateVCF {
          input:
            input_vcf = ExtractAnAcAfFromVCF.output_vcf,
            input_vcf_index = ExtractAnAcAfFromVCF.output_vcf_index,
            output_annotated_file_name = "${input_vcf_name}_annotated",
            nirvana_data_tar = nirvana_data_directory,
            custom_annotations_file = ExtractAnAcAfFromVCF.annotations_file,
        }

        call PrepAnnotationJson {
          input:
            annotation_json = AnnotateVCF.annotation_json,
            output_file_suffix = "${input_vcf_name}.json.gz",
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


        python3 /app/extract_subpop.py \
          --input_path ~{updated_input_ancestry_file} \
          --output_path ~{output_ancestry_filename}

    >>>

    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:ah_var_store_20211007"
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
        File ref
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
    String custom_annotations_file_name = "ac_an_af.tsv"
    String local_input_vcf = basename(input_vcf)
    String local_input_vcf_index = basename(input_vcf_index)
    String normalized_vcf = "normalized.vcf"
    String normalized_vcf_compressed = "normalized.vcf.gz"
    String normalized_vcf_indexed = "normalized.vcf.gz.tbi"

    # separate multi-allelic sites into their own lines, remove deletions and filtered sites and make a sites only vcf
    # while extracting the an/ac/af & sc by subpopulation into a tsv
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
        gsutil cp ~{ref} Homo_sapiens_assembly38.fasta

        ## compare the ancestry sample list with the vcf sample list
        # TODO throw an error (but dont fail the job) if there are samples that are in one, but not the other. Throw two different errors.
        # Currently commented out to save on io with the AoU beta VAT creation
        # awk '{print $1}' ~{subpopulation_sample_list} | tail -n +2 | sort -u > collected_subpopulation_samples.txt
        # bcftools query --list-samples ~{local_input_vcf} | sort -u > collected_samples.txt
        # diff collected_subpopulation_samples.txt collected_samples.txt

        # expected_subpopulations = [
        # "afr",
        # "amr",
        # "eas",
        # "eur",
        # "mid",
        # "oth",
        # "sas"
        #]

        ## track the dropped variants with +500 alt alleles or N's in the reference (Since Nirvana cant handle N as a base, drop them for now)
        bcftools view -i 'N_ALT>500 || REF~"N"' ~{local_input_vcf} | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' > track_dropped.tsv

        ## filter out sites with too many alt alleles
        bcftools view -e 'N_ALT>500 || REF~"N"' --no-update  ~{local_input_vcf} | \
        ## filter out the non-passing sites
        bcftools view  -f 'PASS,.' --no-update | \
        ## normalize, left align and split multi allelic sites to new lines, remove duplicate lines
        bcftools norm -m- --check-ref w -f Homo_sapiens_assembly38.fasta | \
        ## filter out spanning deletions and variants with an AC of 0
        bcftools view  -e 'ALT[0]="*" || AC=0' --no-update | \
        ## ensure that we respect the FT tag
        bcftools filter -i "FORMAT/FT='PASS,.'" --set-GTs . > ~{normalized_vcf}

        ## clean up unneeded file
        rm ~{local_input_vcf}

        ## make a file of just the first 5 columns of the tsv
        bcftools query ~{normalized_vcf} -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' > check_duplicates.tsv
        ## check it for duplicates and put them in a new file
        sort check_duplicates.tsv | uniq -d | cut -f1,2,3,4,5  > duplicates.tsv
        rm check_duplicates.tsv ## clean up
        ## remove those rows (that match up to the first 5 cols)
        grep -v -wFf duplicates.tsv ~{normalized_vcf} > deduplicated.vcf
        rm ~{normalized_vcf} ## clean up

        ## add duplicates to the file tracking dropped variants
        cat duplicates.tsv >> track_dropped.tsv
        rm duplicates.tsv ## clean up unneeded file

        ## calculate annotations for all subpopulations
        bcftools plugin fill-tags  -- deduplicated.vcf -S ~{subpopulation_sample_list} -t AC,AF,AN,AC_het,AC_hom,AC_Hemi | bcftools query -f \
        '%CHROM\t%POS\t%REF\t%ALT\t%AC\t%AN\t%AF\t%AC_Hom\t%AC_Het\t%AC_Hemi\t%AC_afr\t%AN_afr\t%AF_afr\t%AC_Hom_afr\t%AC_Het_afr\t%AC_Hemi_afr\t%AC_amr\t%AN_amr\t%AF_amr\t%AC_Hom_amr\t%AC_Het_amr\t%AC_Hemi_amr\t%AC_eas\t%AN_eas\t%AF_eas\t%AC_Hom_eas\t%AC_Het_eas\t%AC_Hemi_eas\t%AC_eur\t%AN_eur\t%AF_eur\t%AC_Hom_eur\t%AC_Het_eur\t%AC_Hemi_eur\t%AC_mid\t%AN_mid\t%AF_mid\t%AC_Hom_mid\t%AC_Het_mid\t%AC_Hemi_mid\t%AC_oth\t%AN_oth\t%AF_oth\t%AC_Hom_oth\t%AC_Het_oth\t%AC_Hemi_oth\t%AC_sas\t%AN_sas\t%AF_sas\t%AC_Hom_sas\t%AC_Het_sas\t%AC_Hemi_sas\n' \
        >> ~{custom_annotations_file_name}

        ## for validation of the pipeline
        wc -l ~{custom_annotations_file_name} | awk '{print $1 -7}'  > count.txt

        ## compress the vcf and index it, make it sites-only
        bcftools view --no-update --drop-genotypes deduplicated.vcf -Oz -o ~{normalized_vcf_compressed}
        ## if we can spare the IO and want to pass a smaller file we can also drop the info field w bcftools annotate -x INFO
        bcftools index --tbi  ~{normalized_vcf_compressed}

    >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:ah_var_store_20211007"
        memory: "32 GB"
        preemptible: 3
        cpu: "2"
        disks: "local-disk 500 SSD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        File annotations_file = "~{custom_annotations_file_name}"
        Int count_variants = read_int("count.txt")
        File track_dropped = "track_dropped.tsv"
        File output_vcf = "~{normalized_vcf_compressed}"
        File output_vcf_index = "~{normalized_vcf_indexed}"
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
        memory: "32 GB"
        cpu: "2"
        preemptible: 5
        disks: "local-disk 500 SSD"
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

    ## TODO these temp files do not currently get cleaned up. Some of them may be helpful for recovery.

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
        docker: "us.gcr.io/broad-dsde-methods/variantstore:ah_var_store_20211007"
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

    # There are two pre-vat tables. A variant table and a genes table. They are joined together for the vat table
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
