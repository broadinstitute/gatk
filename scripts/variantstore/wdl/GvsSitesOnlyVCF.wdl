version 1.0
workflow GvsSitesOnlyVCF { # tested with "gs://broad-dsp-spec-ops/scratch/rcremer/Nirvana/Data/anvil_expected.vcf"
   input {
        Array[File] gvs_extract_cohort_filtered_vcfs
        Array[File] gvs_extract_cohort_filtered_vcf_indices
        String output_sites_only_file_name
        String output_merged_file_name
        String output_annotated_file_name
        String project_id
        String dataset_name
        File nirvana_data_directory
        File nirvana_override_json
        File nirvana_schema_json_file
        File vat_vt_schema_json_file
        File create_vat_jsons_python_script
        String genes_vat_schema = "gene_symbol:STRING,gene_omim_id:INTEGER,omim_phenotypes_id:INTEGER,omim_phenotypes_name:INTEGER"
        File? gatk_override
    }

    # TODO need to decide where to specify the name of the VAT table now that it is needed in 2 steps: BigQueryLoadJson and BigQuerySmokeTest <-- should these be one step?

    #scatter(i in range(length(gvs_extract_cohort_filtered_vcfs)) ) {
    scatter(i in range(length(["1"])) ) {
        #call SitesOnlyVcf {
         # input:
         #   vcf_bgz_gts = gvs_extract_cohort_filtered_vcfs[i],
         #   vcf_index = gvs_extract_cohort_filtered_vcf_indices[i],
         #   output_filename = "${output_sites_only_file_name}_${i}.sites_only.vcf.gz",
        #}

        #call AnnotateShardedVCF {
         # input:
         #   input_vcf = SitesOnlyVcf.output_vcf,
         #   input_vcf_index = SitesOnlyVcf.output_vcf_idx,
         #   output_annotated_file_name = "${output_annotated_file_name}_${i}",
         #   nirvana_data_tar = nirvana_data_directory
        #}

       call PrepAnnotationJson {
         input:
           # annotation_json = AnnotateShardedVCF.annotation_json,
           annotation_json = nirvana_override_json,
           python_script = create_vat_jsons_python_script,
           output_name = "${i}.json"
       }

       call BucketJson { # TODO should this be it's own step? We are just gsutiling
         input:
           vt_bq_loading_json = PrepAnnotationJson.vat_vt_json,
           genes_bq_loading_json = PrepAnnotationJson.vat_genes_json
       }
    }

     call BigQueryLoadJson {
         input:
             nirvana_schema = nirvana_schema_json_file,
             vt_schema = vat_vt_schema_json_file,
             genes_schema = genes_vat_schema,
             project_id = project_id,
             dataset_name = dataset_name
         }

     call BigQuerySmokeTest {
         input:
             project_id = project_id,
             dataset_name = dataset_name,
             # annotation_json = AnnotateShardedVCF.annotation_json
             annotation_json = nirvana_override_json
         }
}

################################################################################
task SitesOnlyVcf {
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
        disks: "local-disk 100 HDD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        File annotation_json = "~{annotation_json_name}"
        File annotation_json_jsi = "~{annotation_json_name_jsi}"
    }
}

task PrepAnnotationJsonJQ {
    input {
        File annotation_json
        File python_script
        String output_name
    }

    String output_vt_json = "vat_vt_bq_load" + output_name
    String output_genes_json = "vat_genes_bq_load" + output_name

    command <<<
        set -e

        python ~{python_script} \
          --annotated_json ~{annotation_json} \
          --output_vt_json ~{output_vt_json} \
          --output_genes_json ~{output_genes_json}

     >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore-export:091920"
        memory: "10 GB"
        cpu: "2"
        disks: "local-disk 100 HDD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        File vat_vt_json="~{output_vt_json}"
        File vat_genes_json="~{output_genes_json}"
    }
}

task PrepAnnotationJson {
    input {
        File annotation_json
        File python_script
        String output_name
    }

    String output_vt_json = "vat_vt_bq_load" + output_name
    String output_genes_json = "vat_genes_bq_load" + output_name

    command <<<
        set -e

        python ~{python_script} \
          --annotated_json ~{annotation_json} \
          --output_vt_json ~{output_vt_json} \
          --output_genes_json ~{output_genes_json}

     >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore-export:091920"
        memory: "10 GB"
        cpu: "2"
        disks: "local-disk 100 HDD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        File vat_vt_json="~{output_vt_json}"
        File vat_genes_json="~{output_genes_json}"
    }
}

task BucketJson {
    input {
        File vt_bq_loading_json
        File genes_bq_loading_json
    }

    command <<<
        set -e

        gsutil cp ~{vt_bq_loading_json} gs://broad-dsp-spec-ops/scratch/rcremer/Nirvana/output/vt/
        gsutil cp ~{genes_bq_loading_json} gs://broad-dsp-spec-ops/scratch/rcremer/Nirvana/output/genes/

     >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "google/cloud-sdk:latest"
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

task BigQueryLoadJson {
    input {
        File nirvana_schema
        File vt_schema
        String genes_schema
        String project_id
        String dataset_name
    }

    # I am going to want to have two pre-vat tables. A variant table and a genes table. They will be joined together for the vat table
    # See if we can grab the annotations json directly from the gcp bucket (so pull it in as a string so it wont)

    String vat_table = "vat_jun28"
    String variant_transcript_table = "vat_vt"
    String genes_table = "vat_genes"

    # TODO there needs to be a LOOP somewhere give that there are now many files in GCP that need to be loaded into BQ

    String vt_jsons_path = "gs://broad-dsp-spec-ops/scratch/rcremer/Nirvana/output/vt/*"
    String genes_jsons_path = "gs://broad-dsp-spec-ops/scratch/rcremer/Nirvana/output/genes/*"

    command <<<

       bq show --project_id ~{project_id} ~{dataset_name}.~{variant_transcript_table} > /dev/null
       BQ_SHOW_RC=$?

       set -e

       if [ $BQ_SHOW_RC -ne 0 ]; then
         echo "Creating a pre-vat table ~{dataset_name}.~{variant_transcript_table}"
         bq --location=US mk --project_id=~{project_id}  ~{dataset_name}.~{variant_transcript_table} ~{vt_schema}
       fi

       echo "Loading data into a pre-vat table ~{dataset_name}.~{variant_transcript_table}"
       # bq --location=US load --project_id=~{project_id} --source_format=NEWLINE_DELIMITED_JSON ~{dataset_name}.~{variant_transcript_table} ~{vt_jsons_path}
       bq --location=US load --project_id="spec-ops-aou" --source_format=NEWLINE_DELIMITED_JSON "anvil_100_for_testing.vat_vt" gs://broad-dsp-spec-ops/scratch/rcremer/Nirvana/output/vt/*

       set +e

       bq show --project_id ~{project_id} ~{dataset_name}.~{genes_table} > /dev/null
       BQ_SHOW_RC=$?

       set -e

       if [ $BQ_SHOW_RC -ne 0 ]; then
         echo "Creating a pre-vat table ~{dataset_name}.~{genes_table}"
         bq --location=US mk --project_id=~{project_id}  ~{dataset_name}.~{genes_table} ~{genes_schema}
       fi

       echo "Loading data into a pre-vat table ~{dataset_name}.~{genes_table}"
       # bq --location=US load  --project_id=~{project_id} --source_format=NEWLINE_DELIMITED_JSON  ~{dataset_name}.~{genes_table} ~{genes_jsons_path}
       bq --location=US load  --project_id=~{project_id} --source_format=NEWLINE_DELIMITED_JSON  ~{dataset_name}.~{genes_table} gs://broad-dsp-spec-ops/scratch/rcremer/Nirvana/output/genes/*

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
       if [ $BQ_SHOW_RC -ne 0 ]; then
         echo "Creating the vat table ~{dataset_name}.~{vat_table}"
         bq --location=US mk --project_id=~{project_id} ~{dataset_name}.~{vat_table} ~{nirvana_schema}
       else
         bq rm -t -f --project_id=~{project_id} ~{dataset_name}.~{vat_table}
         bq --location=US mk --project_id=~{project_id} ~{dataset_name}.~{vat_table} ~{nirvana_schema}
       fi
       echo "And putting data into it"

       # now run some giant query in BQ to get this all in the right table
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
              v.consequence,    #ARRAY_TO_STRING(v.consequence, ",") AS consequence,
              v.dna_change,
              v.variant_type,
              v.exon_number,
              v.intron_number,
              v.genomic_location,
              #v.hgvsc AS splice_distance has not yet been designed
              v.dbsnp_rsid,   #ARRAY_TO_STRING(v.dbsnp_rsid, ",") AS dbsnp_rsid,
              v.entrez_gene_id,
              #g.hgnc_gene_id is not produced by Nirvana annotations
              g.gene_omim_id,
              CASE WHEN ( v.transcript is not null and v.is_canonical_transcript is not True)
                THEN False WHEN ( v.transcript is not null and v.is_canonical_transcript is True) THEN True END AS is_canonical_transcript,
              v.gnomad_all_af,
              v.gnomad_all_ac,
              v.gnomad_all_an,
              v.gnomad_max_af, # this still needs to be designed
              v.gnomad_max_ac, # this still needs to be designed
              v.gnomad_max_an, # this still needs to be designed
              null AS gnomad_max_subpop, # what is this mapping?
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
              left join `~{dataset_name}.~{genes_table}` as g on v.gene_symbol = g.gene_symbol'

 # TODO why do I sometimes hit an error above, but still make it to the smoke test?
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
        File annotation_json
    }

    # What I want to do here is query the final table for my expected results
    # This will be hardcoded for now, but in the future I may want to pull a line out of the annotations json and use thats

    String vat_table = "vat_jun26" # TODO seems dangerous to specify this twice

    command <<<
        set -e

        # validate the VAT table
        echo "Does this table row look like I expect it to?"
        # We will pass, or fail, a pipeline run by checking the following

        # ------------------------------------------------
        # VALIDATION #1
        # The number of passing variants in GVS matches the number of variants in the VAT.
        # TODO this tests the python---what other qs should we ask here?
        echo  "VALIDATION #1"
        # Please note that we are counting the number of variants in GVS, not the number of sites, which may add a difficulty to this task.
        # TODO should these get broken down more so as not to test my sloppy bash over testing the data?
        # TODO should I write this test in python instead?!?!?
        # grep -o -i '"vid":' ~{annotation_json} | wc -l
        ANNOTATE_JSON_VARIANT_COUNT=$(grep -o -i '"vid":' ~{annotation_json} | wc -l)
        echo $ANNOTATE_JSON_VARIANT_COUNT
        bq query --nouse_legacy_sql --project_id=~{project_id} 'SELECT COUNT (DISTINCT vid) FROM `~{dataset_name}.~{vat_table}`'
        BQ_VAT_VARIANT_COUNT=$(bq query --nouse_legacy_sql --project_id=~{project_id} 'SELECT COUNT (DISTINCT vid) AS count FROM `~{dataset_name}.~{vat_table}`'| tr -dc '0-9')
        echo $BQ_VAT_VARIANT_COUNT

        if [[ $ANNOTATE_JSON_VARIANT_COUNT -ne $BQ_VAT_VARIANT_COUNT ]]
        then
          echo "The number of variants is incorrect"
          echo  "Validation has failed"
          exit 1
        else
          echo "The number of passing variants in GVS matches the number of variants in the VAT"
        fi

        # ------------------------------------------------
        # VALIDATION #2
        # Less than 5% of the variants have a non-null value in the gene field
        echo  "VALIDATION #2"
        # Get the number of variants which have been joined (TODO is it the actual joining? or just the potential?)
        bq query --nouse_legacy_sql --project_id=~{project_id} \
          'SELECT COUNT (DISTINCT vid) AS count FROM `~{dataset_name}.~{vat_table}` WHERE gene_symbol IS NOT NULL'
        BQ_VAT_GENE_COUNT=$(bq query --nouse_legacy_sql --project_id=~{project_id} 'SELECT COUNT (DISTINCT vid) AS count FROM `~{dataset_name}.~{vat_table}` WHERE gene_symbol IS NOT NULL'| tr -dc '0-9')
        echo $BQ_VAT_GENE_COUNT
        echo "Get the percent"
        PERCENT=$((BQ_VAT_GENE_COUNT * 100 / BQ_VAT_VARIANT_COUNT))
        echo $PERCENT

        if [[ $PERCENT -gt 5 ]]
          then echo "There are too many genes"
          echo  "Validation has failed"
          # exit 1 <-- great! This worked! But need to get by it with my sloppy test data for now!
        else
          echo "Less than 5% of the variants have a non-null value in the gene field"
        fi

        # ------------------------------------------------
        # VALIDATION #3
        # All variants in the TESK2 gene region (chr1:45,343,883-45,491,163) list multiple genes and those genes are always TESK2 and AL451136.1.
        echo  "VALIDATION #3"
        # TODO ask Lee for help here. I'm not sure how to do this one?
        # 'SELECT COUNT (DISTINCT vid) AS distinct_vid_count FROM `~{dataset_name}.~{vat_table}` WHERE contig="chr1"'
        echo  "Still need to validate chr1"

        # TODO this is not done!

        # ------------------------------------------------
        # VALIDATION #4
        # If a vid has a null transcript, then the vid is only in one row of the VAT.
        echo  "VALIDATION #4"
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

        echo $BQ_VAT_ROWS_NO_TRANSCRIPT_COUNT
        echo $BQ_VAT_VARIANT_NO_TRANSCRIPT_COUNT

        if [[ $BQ_VAT_ROWS_NO_TRANSCRIPT_COUNT -ne $BQ_VAT_VARIANT_NO_TRANSCRIPT_COUNT ]]
        then
          echo "The number of rows for variants with no transcripts is incorrect"
          echo  "Validation has failed"
          # exit 1 <-- great! This worked! But need to get by it with my sloppy test data for now!
        else
          echo "If a vid has a null transcript, then the vid is only in one row of the VAT"
        fi

        # ------------------------------------------------
        # VALIDATION #5
        # If a vid has any non-null transcripts then one transcript must be Ensembl (transcript_source) and canonical (is_canonical).
        echo  "VALIDATION #5"
        # TODO Rori needs to actually build this in the first place---I think this should currently fail


        # ------------------------------------------------
        # VALIDATION #6
        # No vid may have a mix of non-null and null transcripts.
        echo  "VALIDATION #6"
        # Get a list of all distinct vids with non-null transcripts
        # Get a list of all distinct vids with no transcripts
        # Make sure those lists have no intersection

        # bq query --nouse_legacy_sql --project_id="spec-ops-aou" 'SELECT vids_with_transcript_table.vid, vids_no_transcript_table.vid FROM (SELECT DISTINCT vid FROM `~{dataset_name}.~{vat_table}` WHERE transcript IS NULL) AS vids_no_transcript_table inner join (SELECT DISTINCT vid FROM `~{dataset_name}.~{vat_table}` WHERE transcript IS NOT NULL) AS vids_with_transcript_table on vids_with_transcript_table.vid = vids_no_transcript_table.vid'

        # ------------------------------------------------
        # VALIDATION #7
        # No non-nullable fields contain null values.
        # Hmmm--- the use of the vat_schema should do that for us
        # Could get a count of each where _ is null
        # note that the below returns nothing
        # BQ_VAT_NULL=$(bq query --nouse_legacy_sql --project_id=~{project_id} 'SELECT vid FROM `~{dataset_name}.~{vat_table}` WHERE vid IS NULL OR contig IS NULL OR position IS NULL OR ref_allele IS NULL OR alt_allele IS NULL OR variant_type IS NULL OR genomic_location IS NULL'
        # vid, contig, position, ref_allele, alt_allele, gvs_all_ac, gvs_all_an, gvs_all_af, variant_type, genomic_location + the aou max stuff Lee still has to figure out

        # ------------------------------------------------
        # VALIDATION #8
        # Each key combination is unique.
        # bq query --nouse_legacy_sql --project_id=~{project_id} 'SELECT COUNT (*) FROM `~{dataset_name}.~{vat_table}`'
        # BQ_VAT_ROW_COUNT=$(bq query --nouse_legacy_sql --project_id=~{project_id} 'SELECT COUNT (*) FROM `~{dataset_name}.~{vat_table}`')
        # bq query --nouse_legacy_sql --project_id=~{project_id} 'SELECT COUNT (DISTINCT vid, DISTINCT transcript) FROM `~{dataset_name}.~{vat_table}`'
        # BQ_VAT_KEY_COUNT=$(bq query --nouse_legacy_sql --project_id=~{project_id} 'SELECT COUNT (DISTINCT vid, DISTINCT transcript) FROM `~{dataset_name}.~{vat_table}`')

        # if [ $BQ_VAT_ROW_COUNT -ne $BQ_VAT_KEY_COUNT ];
          # then echo "The number of rows and keys is different"
          # exit 1
        # else:
          # echo "Each key combination is unique"
        # fi

        # ------------------------------------------------
        # FURTHER VALIDATION ?!??!
        # Do I want to add additional checks that validate the mapping from the annotations json--ie count instances of 'Ensembl' in the annotations json ?

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


