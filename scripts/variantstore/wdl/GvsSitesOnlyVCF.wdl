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
        File nirvana_schema_json_file
        File vat_vt_schema_json_file
        File vat_genes_schema_json_file
        String output_path # TODO Is there a Path wdl type?
        String table_id

        File? service_account_json
        File? gatk_override
    }

    # TODO need to decide where to specify the name of the VAT table now that it is needed in 2 steps: BigQueryLoadJson and BigQuerySmokeTest <-- should these be one step?

    scatter(i in range(length(gvs_extract_cohort_filtered_vcfs)) ) {
        call SitesOnlyVcf {
          input:
            input_vcf = gvs_extract_cohort_filtered_vcfs[i],
            input_vcf_index = gvs_extract_cohort_filtered_vcf_indices[i],
            service_account_json = service_account_json,
            output_filename = "${output_sites_only_file_name}_${i}.sites_only.vcf.gz",
        }

        # we want a step here that uses bcftools to create a TSV for ac, an and af each and then we will nirvana-ize it and then use it for annotations

        call AnnotateShardedVCF {
          input:
            input_vcf = SitesOnlyVcf.output_vcf,
            input_vcf_index = SitesOnlyVcf.output_vcf_idx,
            output_annotated_file_name = "${output_annotated_file_name}_${i}",
            nirvana_data_tar = nirvana_data_directory
        }

       call PrepAnnotationJson {
         input:
           annotation_json = AnnotateShardedVCF.annotation_json,
           output_name = "${i}.json.gz",
           output_path = output_path
       }
    }
     call BigQueryLoadJson {
         input:
             nirvana_schema = nirvana_schema_json_file,
             vt_schema = vat_vt_schema_json_file,
             genes_schema = vat_genes_schema_json_file,
             project_id = project_id,
             dataset_name = dataset_name,
             output_path = output_path,
             table_id = table_id,
             prep_jsons_done = PrepAnnotationJson.done
         }

     call BigQuerySmokeTest {
         input:
             project_id = project_id,
             dataset_name = dataset_name,
             table_id = table_id,
             annotation_jsons = AnnotateShardedVCF.annotation_json,
             load_jsons_done = BigQueryLoadJson.done
         }
}

################################################################################
task SitesOnlyVcf {
    input {
        File input_vcf
        File input_vcf_index
        File? service_account_json
        String output_filename
    }
    String output_vcf_idx = basename(output_filename) + ".tbi" # or will this be .idx if from .vcf.gz? or ".tbi" if a .vcf

    String has_service_account_file = if (defined(service_account_json)) then 'true' else 'false'
    String input_vcf_basename = basename(input_vcf)
    String updated_input_vcf = if (defined(service_account_json)) then input_vcf_basename else input_vcf

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
            export GOOGLE_APPLICATION_CREDENTIALS=~{service_account_json}
            gcloud auth activate-service-account --key-file='~{service_account_json}'

            gsutil cp ~{input_vcf} .
            gsutil cp ~{input_vcf_index} .
        fi

        gatk --java-options "-Xmx2048m" \
            SelectVariants \
                -V ~{updated_input_vcf} \
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

task AnnotateShardedVCF {
    input {
        File input_vcf
        File input_vcf_index
        String output_annotated_file_name
        File nirvana_data_tar
    }
    String annotation_json_name = output_annotated_file_name + ".json.gz"
    String annotation_json_name_jsi = annotation_json_name + ".jsi"

    String nirvana_location = "/opt/nirvana/Nirvana.dll"
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
        docker: "annotation/nirvana:3.14" # this download is too slow---can we beef this up?
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
        String output_name
        String output_path
    }

    String output_vt_json = "vat_vt_bq_load" + output_name
    String output_genes_json = "vat_genes_bq_load" + output_name
    String output_vt_gcp_path = output_path + 'vt/'
    String output_genes_gcp_path = output_path + 'genes/'
    String output_ant_gcp_path = output_path + 'annotations/'

    command <<<
        set -e

        python3 /app/create_variant_annotation_table.py \
          --annotated_json ~{annotation_json} \
          --output_vt_json ~{output_vt_json} \
          --output_genes_json ~{output_genes_json}

        gsutil cp ~{output_vt_json} '~{output_vt_gcp_path}'
        gsutil cp ~{output_genes_json} '~{output_genes_gcp_path}'

     >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:ah_var_store_20200709"
        memory: "10 GB"
        cpu: "2"
        disks: "local-disk 100 HDD"
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
    input {
        File nirvana_schema
        File vt_schema
        File genes_schema
        String project_id
        String dataset_name
        String output_path
        String table_id
        File? service_account_json
        Array[String] prep_jsons_done
    }

    # I am going to want to have two pre-vat tables. A variant table and a genes table. They will be joined together for the vat table
    # See if we can grab the annotations json directly from the gcp bucket (so pull it in as a string so it wont)

    String vat_table = "vat_" + table_id
    String variant_transcript_table = "vat_vt_"  + table_id
    String genes_table = "vat_genes_" + table_id
    String vt_path = output_path + 'vt/*'
    String genes_path = output_path + 'genes/*'

    String has_service_account_file = if (defined(service_account_json)) then 'true' else 'false'

    command <<<

       if [ ~{has_service_account_file} = 'true' ]; then
           export GOOGLE_APPLICATION_CREDENTIALS=~{service_account_json}
           gcloud auth activate-service-account --key-file='~{service_account_json}'
           gcloud config set project ~{project_id}
       fi

       bq show --project_id ~{project_id} ~{dataset_name}.~{variant_transcript_table} > /dev/null
       BQ_SHOW_RC=$?

       set -e

       if [ $BQ_SHOW_RC -ne 0 ]; then # do we want to delete if it's still there?  no shards on this part of the workflow
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

       # TODO check the below logic. Seems to be appending.....?
       # TODO how do we want to be doing this? With a hash of some sort?

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

       bq query --nouse_legacy_sql --destination_table=~{dataset_name}.~{vat_table} --project_id=~{project_id} \
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
              # this is the first value in spliceAI (need to validate that there will only ever be one)
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
              v.clinvar_classification, # validate that we only grab lines WHERE id LIKE "RCV%"
              v.clinvar_last_updated,
              v.clinvar_phenotype,
              FROM `~{dataset_name}.~{variant_transcript_table}` as v
              left join
              (SELECT gene_symbol, ANY_VALUE(gene_omim_id) AS gene_omim_id, ANY_VALUE(omim_phenotypes_id) AS omim_phenotypes_id, ANY_VALUE(omim_phenotypes_name) AS omim_phenotypes_name FROM `~{dataset_name}.~{genes_table}` group by gene_symbol) as g
              on v.gene_symbol = g.gene_symbol'

 # TODO why do I sometimes hit an error above, but still make it to the smoke test?!?!??! Seems like the auth errors dont fail the workflow
  >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "openbridge/ob_google-bigquery:latest"
        memory: "3 GB"
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
        Array[File] annotation_jsons
        Boolean load_jsons_done
        String table_id
    }

    # What I want to do here is query the final table for my expected results
    # This will be hardcoded for now, but in the future I may want to pull a line out of the annotations json and use thats

    String vat_table = "vat_" + table_id

    command <<<
        set +e

        # validate the VAT table
        # We will pass, or fail, a pipeline run by checking the following

        # TODO do we want to bork the pipeline if one of these following rules fail? Or simply collect the errors?

        # ------------------------------------------------
        # VALIDATION #1
        # TODO sometimes when there's an error msg, this is stripping out the #s and doing math with that like a dummy
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
          # exit 1  <-- collect the errors? or fail the pipeline?
        else
          echo "The VAT has been updated"
          echo "Validation #1 has passed"
        fi

        # ------------------------------------------------
        # VALIDATION #2
        # The number of passing variants in GVS matches the number of variants in the VAT.
        # TODO this tests the python---what other qs should we ask here?
        echo  "VALIDATION #2"
        # Please note that we are counting the number of variants in GVS, not the number of sites, which may add a difficulty to this task.
        # TODO should these get broken down more so as not to test my sloppy bash over testing the data?
        # TODO should I write this test in python instead?!?!? YES! I think I should!!!!
        ANNOTATE_JSON_VARIANT_COUNT=$(gunzip -dc ~{sep=" " annotation_jsons} | grep -o -i '"vid":' | wc -l)
        echo $ANNOTATE_JSON_VARIANT_COUNT

        if [[ $ANNOTATE_JSON_VARIANT_COUNT -ne $BQ_VAT_VARIANT_COUNT ]]
        then
          echo "The number of variants is incorrect"
          echo  "Validation has failed"
          # exit 1  <-- collect the errors? or fail the pipeline?
        else
          echo "The number of passing variants in GVS matches the number of variants in the VAT"
          echo "Validation #2 has passed"
        fi

        # ------------------------------------------------
        # VALIDATION #3
        # Less than 5% of the variants have a non-null value in the gene field
        # TODO CHECKING WITH LEE --> We dont think 5% is correct
        echo  "VALIDATION #3"
        # Get the number of variants which have been joined (TODO is it the actual joining? or just the potential?)
        bq query --nouse_legacy_sql --project_id=~{project_id} \
          'SELECT COUNT (DISTINCT vid) AS count FROM `~{dataset_name}.~{vat_table}` WHERE gene_symbol IS NOT NULL'
        BQ_VAT_GENE_COUNT=$(bq query --nouse_legacy_sql --project_id=~{project_id} 'SELECT COUNT (DISTINCT vid) AS count FROM `~{dataset_name}.~{vat_table}` WHERE gene_symbol IS NOT NULL'| tr -dc '0-9')
        echo "Get the percent"
        PERCENT=$((BQ_VAT_GENE_COUNT * 100 / BQ_VAT_VARIANT_COUNT))
        echo $PERCENT

        if [[ $PERCENT -gt 5 ]]
          then echo "There are too many genes"
          echo  "Validation has failed"
          # exit 1  <-- collect the errors? or fail the pipeline?
        else
          echo "Less than 5% of the variants have a non-null value in the gene field"
          echo "Validation #3 has passed"
        fi

        # ------------------------------------------------
        # VALIDATION #4
        # All variants in the TESK2 gene region (chr1:45,343,883-45,491,163) list multiple genes and those genes are always TESK2 and AL451136.1.
        echo  "VALIDATION #4"
        # TODO ask Lee for help here. I'm not sure how to do this one?
        # 'SELECT COUNT (DISTINCT vid) AS distinct_vid_count FROM `~{dataset_name}.~{vat_table}` WHERE contig = "chr1" and position >= 45343883 and position <= 45491163'
        POSITIONAL_COUNT=$('SELECT COUNT (DISTINCT vid) AS distinct_vid_count FROM `~{dataset_name}.~{vat_table}` WHERE contig = "chr1" and position >= 45343883 and position <= 45491163')
        TESK2_COUNT=$('SELECT COUNT (DISTINCT vid) AS distinct_vid_count FROM `~{dataset_name}.~{vat_table}` WHERE contig = "chr1" and position >= 45343883 and position <= 45491163 and gene_symbol="TESK2"')
        AL451136_COUNT=$('SELECT COUNT (DISTINCT vid) AS distinct_vid_count FROM `~{dataset_name}.~{vat_table}` WHERE contig = "chr1" and position >= 45343883 and position <= 45491163 and gene_symbol="AL451136.1"')
        echo  "Still need to validate chr1"


        if [[ $POSITIONAL_COUNT -ne ($TESK2_COUNT + $AL451136_COUNT) ]]
          then echo "There are unexpected genes in the TESK2 region"
          echo  "Validation has failed"
          # exit 1  <-- collect the errors? or fail the pipeline?
        else
          echo "All variants in the TESK2 gene region list multiple genes and those genes are always TESK2 and AL451136.1"
          echo "Validation #4 has passed"
        fi

        # ------------------------------------------------
        # VALIDATION #5
        # If a vid has a null transcript, then the vid is only in one row of the VAT.
        echo  "VALIDATION #5"
        # Get a count of all the rows in the vat with no transcript
        # Get a count of all distinct VID in vat with no transcript
        # Make sure they are the same
        # TODO we could in the future specify which VIDs are not distinct
        bq query --nouse_legacy_sql --project_id=~{project_id} \
          'SELECT COUNT (DISTINCT vid) AS distinct_vid_count FROM `~{dataset_name}.~{vat_table}` WHERE transcript IS NULL'
        BQ_VAT_ROWS_NO_TRANSCRIPT_COUNT=$(bq query --nouse_legacy_sql --project_id=~{project_id} 'SELECT COUNT (DISTINCT vid) AS distinct_vid_count FROM `~{dataset_name}.~{vat_table}` WHERE transcript IS NULL'| tr -dc '0-9')
        bq query --nouse_legacy_sql --project_id=~{project_id} \
          'SELECT COUNT(*) AS count FROM `~{dataset_name}.~{vat_table}` WHERE transcript IS NULL'
        BQ_VAT_VARIANT_NO_TRANSCRIPT_COUNT=$(bq query --nouse_legacy_sql --project_id=~{project_id} 'SELECT COUNT(*) AS count FROM `~{dataset_name}.~{vat_table}` WHERE transcript IS NULL'| tr -dc '0-9')

        if [[ $BQ_VAT_ROWS_NO_TRANSCRIPT_COUNT -ne $BQ_VAT_VARIANT_NO_TRANSCRIPT_COUNT ]]
        then
          echo "The number of rows for variants with no transcripts is incorrect"
          echo  "Validation has failed"
          # exit 1  <-- collect the errors? or fail the pipeline?
        else
          echo "If a vid has a null transcript, then the vid is only in one row of the VAT"
          echo "Validation #5 has passed"
        fi

        # ------------------------------------------------
        # VALIDATION #6
        # If a vid has any non-null transcripts then one transcript must be Ensembl (transcript_source) and canonical (is_canonical).
        echo  "VALIDATION #6"
        # Lee is removing this for now.
        echo "Lee is removing this for now"
        # SELECT vid, string_agg(DISTINCT is_canonical_transcript) as canonical FROM `spec-ops-aou.anvil_100_for_testing.aou_shard_223_vat` where transcript is not null group by vid order by canonical <--there are plenty that are false
        BQ_VAT_ENSEMBL_COUNT=$(bq query --nouse_legacy_sql --project_id=~{project_id} 'SELECT COUNT(*) AS count FROM `~{dataset_name}.~{vat_table}` where transcript is not null and transcript_source="Ensembl"')
        BQ_VAT_TRANSCRIPT_COUNT=$(bq query --nouse_legacy_sql --project_id=~{project_id} 'SELECT COUNT(*) AS count FROM `~{dataset_name}.~{vat_table}` where transcript is not null')

        if [[ $BQ_VAT_ENSEMBL_COUNT -ne $BQ_VAT_TRANSCRIPT_COUNT ]]
        then
          echo "All transcripts should be from Ensembl"
          echo  "Validation has failed"
          # exit 1  <-- collect the errors? or fail the pipeline?
        else
          echo "If a vid has any non-null transcripts then one transcript must be Ensembl"
          echo "Validation #6 has passed"
        fi

        # ------------------------------------------------
        # VALIDATION #7
        # No vid may have a mix of non-null and null transcripts.
        echo  "VALIDATION #7"
        # Get a list of all distinct vids with non-null transcripts
        # Get a list of all distinct vids with no transcripts
        # Make sure those lists have no intersection

        BQ_VAT_VID_OVERLAP=$(bq query --nouse_legacy_sql --project_id="spec-ops-aou" 'SELECT COUNT(vids_with_transcript_table.vid) AS count FROM (SELECT DISTINCT vid FROM `~{dataset_name}.~{vat_table}` WHERE transcript IS NULL) AS vids_no_transcript_table inner join (SELECT DISTINCT vid FROM `~{dataset_name}.~{vat_table}` WHERE transcript IS NOT NULL) AS vids_with_transcript_table on vids_with_transcript_table.vid = vids_no_transcript_table.vid')

        if [[ $BQ_VAT_VID_OVERLAP -gt 0 ]]
          then echo "A vid has a mix of non-null and null transcripts"
          echo  "Validation has failed"
          # exit 1  <-- collect the errors? or fail the pipeline?
        else
          echo "No vids have a mix of non-null and null transcripts"
          echo "Validation #7 has passed"
        fi

        # ------------------------------------------------
        # VALIDATION #8
        # No non-nullable fields contain null values.
        echo  "VALIDATION #8"
        # Hmmm--- the use of the vat_schema should do that for us
        # Could get a count of each where _ is null
        # note that the below returns nothing
        BQ_VAT_NULL=$(bq query --nouse_legacy_sql --project_id=~{project_id} 'SELECT COUNT (DISTINCT vid) as count FROM `~{dataset_name}.~{vat_table}` WHERE vid IS NULL OR contig IS NULL OR position IS NULL OR ref_allele IS NULL OR alt_allele IS NULL OR variant_type IS NULL OR genomic_location IS NULL'
        # vid, contig, position, ref_allele, alt_allele, gvs_all_ac, gvs_all_an, gvs_all_af, variant_type, genomic_location + the aou max stuff Lee still has to figure out
        echo $BQ_VAT_NULL

        if [[ $BQ_VAT_NULL -gt 0 ]]
          then echo "A required value is null"
          echo  "Validation has failed"
          # exit 1  <-- collect the errors? or fail the pipeline?
        else
          echo "No required values are null"
          echo "Validation #8 has passed"
        fi

        # ------------------------------------------------
        # VALIDATION #9
        # Each key combination is unique.
        echo  "VALIDATION #9"
        # get the sum of all the distinct vids where transcript is null and all the distinct transcript where transcript is not null
        BQ_VAT_UNIQUE_IDS_COUNT=$(bq query --nouse_legacy_sql --project_id=~{project_id} 'SELECT COUNT(*) FROM (SELECT vid, transcript FROM `spec-ops-aou.anvil_100_for_testing.vat_jul_6_chr1` group by vid, transcript)')
        BQ_VAT_ROW_COUNT=$(bq query --nouse_legacy_sql --project_id=~{project_id} 'SELECT COUNT (*) FROM `~{dataset_name}.~{vat_table}` ')

        if [ $BQ_VAT_UNIQUE_IDS_COUNT -ne $BQ_VAT_ROW_COUNT ];
           then echo "There are duplicate variant - transcript rows"
           echo  "Validation has failed"
          # exit 1  <-- collect the errors? or fail the pipeline?
        else:
          echo "Each key combination is unique"
          echo "Validation #9 has passed"
        fi

        # ------------------------------------------------
        # FURTHER VALIDATION ?!??!
        # Do I want to add additional checks that validate the mapping from the annotations json--ie count instances of 'Ensembl' in the annotations json ?
        # RCA rule in clinvar data?
        # check that vt and vat are the same length

    >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "openbridge/ob_google-bigquery:latest"
        memory: "1 GB"
        cpu: "1"
        disks: "local-disk 100 HDD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        Boolean done = true
    }
}
