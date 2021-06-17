version 1.0
workflow GvsSitesOnlyVCF {
   input {
        Array[File] gvs_extract_cohort_filtered_vcfs
        String output_sites_only_file_name
        String output_merged_file_name
        String output_annotated_file_name
        String project_id
        String dataset_name
        File nirvana_data_directory
        String nirvana_schema = "position:INTEGER,vid:String,contig:STRING,ref_allele:STRING,alt_allele:STRING,variant_type:STRING,genomic_location:STRING,dbsnp_rsid:STRING,transcript:STRING,gene_symbol:STRING,transcript_source:STRING,aa_change:STRING,consequence:STRING,dna_change:String,exon_number:STRING,intron_number:STRING,splice_distance:STRING,entrez_gene_id:STRING,is_canonical_transcript:STRING,gvs_all_ac:INTEGER,gvs_all_an:INTEGER,gvs_all_af:INTEGER,revel:FLOAT,splice_ai_acceptor_gain_score:INTEGER,splice_ai_acceptor_gain_distance:INTEGER,splice_ai_acceptor_loss_score:INTEGER,splice_ai_acceptor_loss_distance:INTEGER,splice_ai_donor_gain_score:INTEGER,splice_ai_donor_gain_distance:INTEGER,splice_ai_donor_loss_score:INTEGER,splice_ai_donor_loss_distance:INTEGER,clinvar_classification:STRING,clinvar_last_updated:DATE,clinvar_phenotype:STRING,gnomad_all_af:STRING,gnomad_all_ac:STRING,gnomad_all_an:STRING,gnomad_max_af:STRING,gnomad_max_ac:STRING,gnomad_max_an:STRING,gnomad_max_subpop:STRING,gene_omim_id:STRING,omim_phenotypes_id:STRING,omim_phenotypes_name:STRING"
        File? gatk_override
    }

    scatter(i in range(length(gvs_extract_cohort_filtered_vcfs)) ) {
        call SitesOnlyVcf {
                    input:
                        vcf_bgz_gts                 = gvs_extract_cohort_filtered_vcfs[i],
                        output_filename             = "${output_sites_only_file_name}_${i}.sites_only.vcf.gz",
                }
    }

    call MergeVCFs { # why are we merging before annotating? couldn't we keep them separate?
    # I guess if we create the table in one single thread, we can load data w multi (we would be downloading the annotations many times perhaps...)
      input:
          input_vcfs = SitesOnlyVcf.output_vcf,
          input_vcf_indices = SitesOnlyVcf.output_vcf_idx,
          output_merged_file_name = output_merged_file_name,
    }

    call AnnotateVCF {
      input:
          input_vcf = MergeVCFs.merged_vcf,
          output_annotated_file_name = output_annotated_file_name,
          nirvana_data_tar = nirvana_data_directory
    }

    call PrepAnnotationJson {
      input:
          annotation_json = AnnotateVCF.annotation_json,
          annotation_json_jsi = AnnotateVCF.annotation_json_jsi

    }

    call BigQueryLoadJson {
      input:
          annotation_json = PrepAnnotationJson.annotations_edited_file,
          nirvana_schema = nirvana_schema,
          project_id = project_id,
          dataset_name = dataset_name
    }
}



################################################################################
task SitesOnlyVcf {
    input {
        File vcf_bgz_gts
        String output_filename
    }
    String output_vcf_idx = basename(output_filename) + ".tbi" # or will this be .idx if from .vcf.gz?
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

task MergeVCFs {
    input {
        Array[File] input_vcfs
        Array[File] input_vcf_indices
        String output_merged_file_name
    }
    String output_vcf = basename(output_merged_file_name) + ".vcf.gz"
    String output_vcf_idx = basename(output_vcf) + ".tbi"
    command <<<
        set -e
        gatk --java-options "-Xmx2048m" \
            MergeVcfs -I ~{sep=' -I ' input_vcfs} -O ~{output_vcf}

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
        File merged_vcf = "~{output_vcf}"
        File merged_vcf_idx = "~{output_vcf_idx}"
    }
}

task AnnotateVCF {
    input {
        File input_vcf
        String output_annotated_file_name
        File nirvana_data_tar
    }
    String annotation_json_name = basename(output_annotated_file_name) + ".json.gz"
    String annotation_json_name_jsi = annotation_json_name + ".jsi"

    String nirvana_location = "/opt/nirvana/Nirvana.dll"
    String path = "/Cache/GRCh38/Both"
    String path_supplementary_annotations = "/SupplementaryAnnotation/GRCh38"
    String path_reference = "/References/Homo_sapiens.GRCh38.Nirvana.dat"

    command <<<
        set -e

        # NOTE: Validate a lil so that we don't waste time copying down the data sources if there's an error.
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
        memory: "3 GB"
        cpu: "1"
        disks: "local-disk 100 HDD"
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
        File annotation_json_jsi
    }

    String unzipped_json = "unzipped.json"
    String newline_delimited_json = "load_into_bq.json"

    command <<<
        set -e

        # prepare the json file for loading into BQ by making it a new line delimited json
        gunzip -c ~{annotation_json} >> ~{unzipped_json}
        jq -rc '.positions | .[]' ~{unzipped_json} >> ~{newline_delimited_json}

    >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "stedolan/jq:latest"
        memory: "3 GB"
        cpu: "1"
        disks: "local-disk 100 HDD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        File annotations_edited_file = newline_delimited_json
    }
}

task BigQueryLoadJson {
    input {
        File annotation_json
        String nirvana_schema
        String project_id
        String dataset_name
    }

    command <<<
        set -e

        # load the annotation json in as a temp interim BQ vat table
        TEMP_TABLE="~{dataset_name}.pre-vat"
        bq --location=US load --project_id=~{project_id} --source_format=NEWLINE_DELIMITED_JSON --autodetect $TEMP_TABLE ~{annotation_json}

        # create the final vat table with the correct fields
        # TODO this is a hacky set of fields---would ideally use the vat_schema.json
        PARTITION_FIELD="position"
        CLUSTERING_FIELD="vid"
        PARTITION_STRING="" #--range_partitioning=$PARTITION_FIELD,0,4000,4000"
        CLUSTERING_STRING="" #--clustering_fields=$CLUSTERING_FIELD"
        TABLE="~{dataset_name}.vatter"
        SCHEMA="~{nirvana_schema}"
        PROJECT="~{project_id}"
        echo "Creating a vat table..."
        bq --location=US mk ${PARTITION_STRING} ${CLUSTERING_STRING} --project_id=~{project_id} $TABLE $SCHEMA > status_bq_submission
        echo "And putting data into it"

        # now run some giant query in BQ to get this all in the right table
        bq query --destination_table="anvil_100_for_testing.vatter" --project_id="spec-ops-aou" 'SELECT
              v.position,
              v.vid,
              v.chromosome AS contig,
              v.refAllele AS ref_allele,
              v.altAllele AS alt_allele,
              v.variantType AS variant_type,
              v.hgvsg AS genomic_location,
              ARRAY_TO_STRING(v.dbsnp, ",") AS dbsnp_rsid,
              t.transcript,
              t.hgnc AS gene_symbol,
              t.source AS transcript_source,
              t.hgvsp AS aa_change,
              ARRAY_TO_STRING(t.consequence, ",") AS consequence,
              t.hgvsc AS dna_change,
              t.exons AS exon_number,
              t.introns AS intron_number,
              t.hgvsc AS splice_distance,
              t.geneId AS entrez_gene_id,
        CASE WHEN ( t.transcript is not null and t.isCanonical is not True) THEN "false" WHEN ( t.transcript is not null and t.isCanonical is True) THEN "true" END AS is_canonical_transcript,
        null AS gvs_all_ac,
        null AS gvs_all_an,
        null AS  gvs_all_af,
        v.revel.score  AS revel,
              # we just grab the first value in spliceAI (need to validate that there will only ever be one)
              CASE WHEN (select array_length(v.spliceAI)) > 0
                 THEN v.spliceAI[offset(0)].acceptorGainScore END AS splice_ai_acceptor_gain_score,
              CASE WHEN (select array_length(v.spliceAI)) > 0
                 THEN v.spliceAI[offset(0)].acceptorGainDistance END AS splice_ai_acceptor_gain_distance,
              CASE WHEN (select array_length(v.spliceAI)) > 0
                 THEN v.spliceAI[offset(0)].acceptorLossScore END AS splice_ai_acceptor_loss_score,
              CASE WHEN (select array_length(v.spliceAI)) > 0
                 THEN v.spliceAI[offset(0)].acceptorLossDistance END AS splice_ai_acceptor_loss_distance,
              CASE WHEN (select array_length(v.spliceAI)) > 0
                 THEN v.spliceAI[offset(0)].donorGainScore END AS splice_ai_donor_gain_score,
              CASE WHEN (select array_length(v.spliceAI)) > 0
                 THEN v.spliceAI[offset(0)].donorGainDistance END AS splice_ai_donor_gain_distance,
              CASE WHEN (select array_length(v.spliceAI)) > 0
                 THEN v.spliceAI[offset(0)].donorLossScore END AS splice_ai_donor_loss_score,
              CASE WHEN (select array_length(v.spliceAI)) > 0
                 THEN v.spliceAI[offset(0)].donorLossDistance END AS splice_ai_donor_loss_distance,
              (SELECT ARRAY_TO_STRING(significance, ",") FROM v.clinvar WHERE id LIKE "RCV%") AS clinvar_classification,
              (SELECT lastUpdatedDate FROM v.clinvar WHERE id LIKE "RCV%") AS clinvar_last_updated,
              (SELECT ARRAY_TO_STRING(phenotypes, ",") FROM v.clinvar WHERE id LIKE "RCV%") AS clinvar_phenotype,
              v.gnomad.allAf AS gnomad_all_af,
              v.gnomad.allAc AS gnomad_all_ac,
              v.gnomad.allAn AS gnomad_all_an,
              v.gnomad.afrAf AS gnomad_max_af,
              v.gnomad.afrAc AS gnomad_max_ac,
              v.gnomad.afrAn AS gnomad_max_an,
              null AS gnomad_max_subpop, # what is this mapping?
              null AS gene_omim_id,
              null AS omim_phenotypes_id,
              null AS omim_phenotypes_name,
              from (SELECT position, variantline.* FROM `spec-ops-aou.anvil_100_for_testing.pre-vat`, UNNEST(variants) as variantline) as v left join
              (SELECT position, variantline.vid, transcriptline.* FROM `spec-ops-aou.anvil_100_for_testing.pre-vat`, UNNEST(variants) as variantline, UNNEST(variantline.transcripts) as transcriptline) as t on v.vid = t.vid'
       # cat status_bq_submission | tail -n 1 > status_bq_submission_last_line
       # bq_job_id=$(sed 's/.*://' status_bq_submission_last_line)
        #echo $bq_job_id

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


