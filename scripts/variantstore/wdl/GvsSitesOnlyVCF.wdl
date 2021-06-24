version 1.0
workflow GvsSitesOnlyVCF { # tested with "gs://broad-dsp-spec-ops/scratch/rcremer/Nirvana/Data/anvil_expected.vcf"
   input {
        Array[File] gvs_extract_cohort_filtered_vcfs # gonna try this one next: gs://broad-dsp-spec-ops/kcibul/ft/gvs.chr20.vcf.gz,
        Array[File] gvs_extract_cohort_filtered_vcf_indices # gonna try this one next: gs://broad-dsp-spec-ops/kcibul/ft/gvs.chr20.vcf.gz,
        String output_sites_only_file_name
        String output_merged_file_name
        String output_annotated_file_name
        String project_id
        String dataset_name
        File nirvana_data_directory
        File nirvana_schema_json_file
        String? genes_vat_schema = "gene_symbol:STRING,gene_omim_id:INTEGER,omim_phenotypes_id:INTEGER,omim_phenotypes_name:STRING"
        File? gatk_override
    }

    scatter(i in range(length(gvs_extract_cohort_filtered_vcfs)) ) {
        call SitesOnlyVcf {
                    input:
                        vcf_bgz_gts = gvs_extract_cohort_filtered_vcfs[i],
                        vcf_index = gvs_extract_cohort_filtered_vcf_indices[i],
                        output_filename = "${output_sites_only_file_name}_${i}.sites_only.vcf.gz",
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
    #call AnnotateShardedVCF {
    #    input:
    #       input_vcfs = SitesOnlyVcf.output_vcf,
    #       output_annotated_file_name = output_annotated_file_name,
    #       nirvana_data_tar = nirvana_data_directory
    #}

  #  call PrepAnnotationJson { # TODO should this be it's own step? We need to run a python script to prep the jsons for bq
  #    input:
  #        annotation_json = AnnotateVCF.annotation_json,
   # }

   # call BigQueryLoadJson { # If we dont use the above step to prep w a python script-- we are assuming BQ docker has python
    #  input:
     #     annotation_json = AnnotateVCF.annotation_json,
       #   #annotation_json_jsi = AnnotateVCF.annotation_json_jsi
      #    nirvana_schema = nirvana_schema_json_file,
       #   project_id = project_id,
        #  dataset_name = dataset_name
   # }
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

task MergeVCFs {
    input {
        Array[File] input_vcfs
        Array[File] input_vcf_indices
        String output_merged_file_name
    }
    # TODO Ideally we wouldn't merge the vcfs at all, but keep them separated by position
    # and then annotate them separately and then load each of the annotation.jsons into BQ (may involve multiple temp tables that each get mapped to the final vat)

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

task AnnotateShardedVCF {
    input {
        Array[File] input_vcfs
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

        # =======================================
        # Handle our data sources:

        # Extract the tar.gz:
        echo "Extracting annotation data sources tar/gzip file..."
        mkdir datasources_dir
        tar zxvf ~{nirvana_data_tar} -C datasources_dir  --strip-components 2
        DATA_SOURCES_FOLDER="$PWD/datasources_dir"

        # do a bash loop here to single thread annotate each of the files in the array

        for input_vcf in input_vcfs
        do
          dotnet ~{nirvana_location} \
             -c $DATA_SOURCES_FOLDER~{path} \
             --sd $DATA_SOURCES_FOLDER~{path_supplementary_annotations} \
             -r $DATA_SOURCES_FOLDER~{path_reference} \
             -i ~{input_vcfs} \
             -o ~{output_annotated_file_name+"input_vcf letter or Ii or somethinig?"}
        done


    >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "annotation/nirvana:3.14" # this download is too slow---can we beef this up?
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

task AnnotateVCF { # can we add the tar to the docker container?
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
        docker: "annotation/nirvana:3.14" # this download is too slow---can we beef this up?
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

task BigQueryLoadJson {
    input {
        File annotation_json
        File nirvana_schema
        String project_id
        String dataset_name
    }

    # I am going to want to have two pre-vat tables. A variant table and a genes table. They will be joined together for the vat table
    # See if we can grab the annotations json directly from the gcp bucket (so pull it in as a string so it wont)

    String vat_table = "vatter"
    String variant_transcript_table = "vat_vt"
    String genes_table = "vat_genes"

    # instead of the annotation_json we now need two jsons which are created from the annotation_json with a python script
    # TODO run that python script on the annotation_json
    # we will also want specific schemas for each of these tables---right now we are over-using one. Let's get Lee's sign off on the final one
    File variant_transcript_bq_load.json = annotation_json
    File genes_bq_load.json = annotation_json

    command <<<
        set -e

        # load the annotation json in as a temp interim BQ vat table
        # TODO make sure when you load---if you split in some way, that you dont drop "header lines" from secondary shards that have no headers

        # since we expect to add to these tables, --replace=true is not a flag we want. Verify that this will not be loaded with a bq *
        echo "Creating a pre-vat table ~{dataset_name}.~{variant_transcript_table}"
        bq --location=US load --project_id=~{project_id} --source_format=NEWLINE_DELIMITED_JSON --autodetect ~{dataset_name}.~{variant_transcript_table} ~{variant_transcript_bq_load.json}
        # bq --location=US load --replace=true --project_id="spec-ops-aou" --source_format=NEWLINE_DELIMITED_JSON --autodetect "anvil_100_for_testing.temp" /Users/aurora/Desktop/repositories/gatk/hello_did_I_annotate.json

        echo "Creating a pre-vat table ~{dataset_name}.~{genes_table}"
        bq --location=US load  --project_id=~{project_id} --source_format=NEWLINE_DELIMITED_JSON --autodetect ~{dataset_name}.~{genes_table} ~{genes_bq_load.json}

        set +e
        bq show --project_id ~{project_id} ~{dataset_name}.~{vat_table} > /dev/null
        BQ_SHOW_RC=$?
        set -e
        if [ $BQ_SHOW_RC -ne 0 ]; then
          echo "Creating the vat table ~{dataset_name}.~{vat_table}"
          bq --location=US mk --project_id=~{project_id} ~{dataset_name}.~{vat_table} ~{nirvana_schema}
        fi


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
              ARRAY_TO_STRING(v.consequence, ",") AS consequence,
              v.dna_change,
              v.variant_type,
              v.exon_number,
              v.intron_number,
              v.genomic_location,
              #v.hgvsc AS splice_distance has not yet been designed
              ARRAY_TO_STRING(v.dbsnp_rsid, ",") AS dbsnp_rsid,
              v.entrez_gene_id,
              #g.hgnc_gene_id is not produced by Nirvana annotations
              g.gene_omim_id,
              CASE WHEN ( v.transcript is not null and v.is_canonical_transcript is not True)
                THEN "false" WHEN ( v.transcript is not null and v.is_canonical_transcript is True) THEN "true" END AS is_canonical_transcript,
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
              from `~{dataset_name}.~{variant_transcript_table}` as v
              left join `~{dataset_name}.~{genes_table}` as g on v.gene_symbol = g.gene_symbol'

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


