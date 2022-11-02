version 1.0

workflow GvsCreateVATfromVDS {
    input {
        File input_sites_only_vcf
        File custom_annotations_file
        File ancestry_file

        String project_id
        String dataset_name
        String filter_set_name
        String? vat_version

        String output_path
        Int? merge_vcfs_disk_size_override
    }

    Array[String] contig_array = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"]
    File reference = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
    File nirvana_data_directory = "gs://gvs-internal/bigquery-jointcalling/VAT/Nirvana/Nirvana-references-2022-10-07.tgz"

    ## TODO # do we need to shard the sites only VCF?
    # String input_vds_name = basename(input_VDS_file)
    ## TODO: where do we need to validate that there are no hemis?
    String input_vcf_name = "gvs"

    call MakeSubpopulationFilesAndReadSchemaFiles {
        input:
            input_ancestry_file = ancestry_file
    }

    call CreateCustomAnnotationsFile {
        input:
            custom_annotations_header = MakeSubpopulationFilesAndReadSchemaFiles.custom_annotations_template_file,
            custom_annotations_body = custom_annotations_file,
    }

    ## Use Nirvana to annotate the sites-only VCF and include the AC/AN/AF calculations as custom annotations
    call AnnotateVCF {
        input:
            input_vcf = input_sites_only_vcf,
            output_annotated_file_name = "${input_vcf_name}_annotated",
            nirvana_data_tar = nirvana_data_directory,
            custom_annotations_file = CreateCustomAnnotationsFile.custom_annotations,
    }

    call PrepAnnotationJson {
        input:
            annotation_json = AnnotateVCF.annotation_json,
            output_file_suffix = "${input_vcf_name}.json.gz",
            output_path = output_path,
    }



    call BigQueryLoadJson {
        input:
            nirvana_schema = MakeSubpopulationFilesAndReadSchemaFiles.vat_schema_json_file,
            vt_schema = MakeSubpopulationFilesAndReadSchemaFiles.variant_transcript_schema_json_file,
            genes_schema = MakeSubpopulationFilesAndReadSchemaFiles.genes_schema_json_file,
            project_id = project_id,
            dataset_name = dataset_name,
            output_path = output_path,
            filter_set_name = filter_set_name,
            vat_version = vat_version,
            prep_jsons_done = PrepAnnotationJson.done,
    }


    scatter(i in range(length(contig_array)) ) {
        call BigQueryExportVat {
            input:
                contig = contig_array[i],
                project_id = project_id,
                dataset_name = dataset_name,
                output_path = output_path,
                vat_table = BigQueryLoadJson.vat_table_name,
                load_jsons_done = BigQueryLoadJson.done
        }
    }

    call MergeVatTSVs {
        input:
            export_done = BigQueryExportVat.done,
            contig_array = contig_array,
            output_path = output_path,
            project_id = project_id,
            merge_vcfs_disk_size_override = merge_vcfs_disk_size_override
    }

    output {
        File final_tsv_file = MergeVatTSVs.tsv_file
    }
}

################################################################################

task MakeSubpopulationFilesAndReadSchemaFiles {
    input {
        File input_ancestry_file

        String schema_filepath = "/data/variant_annotation_table/schema/"
        String vat_schema_json_filename = "vat_schema.json"
        String variant_transcript_schema_json_filename = "vt_schema.json"
        String genes_schema_json_filename = "genes_schema.json"
    }
    String output_ancestry_filename =  "ancestry_mapping.tsv"
    String custom_annotations_template_filename =  "custom_annotations_template.tsv"

    command <<<
        set -e

        cp ~{schema_filepath}* .

        ## the ancestry file is processed down to a simple mapping from sample to subpopulation
        python3 /app/extract_subpop.py \
            --input_path ~{input_ancestry_file} \
            --output_path ~{output_ancestry_filename} \
            --custom_annotations_template_path ~{custom_annotations_template_filename}
    >>>

    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:rc_616_var_store_2022_10_25"
        memory: "1 GB"
        preemptible: 3
        cpu: "1"
        disks: "local-disk 100 HDD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        File vat_schema_json_file = vat_schema_json_filename
        File variant_transcript_schema_json_file = variant_transcript_schema_json_filename
        File genes_schema_json_file = genes_schema_json_filename

        File ancestry_mapping_list = output_ancestry_filename
        File custom_annotations_template_file = custom_annotations_template_filename
    }
}


task CreateCustomAnnotationsFile {
    input {
        File custom_annotations_header
        File custom_annotations_body
    }

    String custom_annotations_file_name = "ac_an_af.tsv"

    command <<<
        set -e

        # I know this is gross, but I just want to get this working:
        gsutil cp "gs://fc-d5e319d4-b044-4376-afde-22ef0afc4088/submissions/8a5b026d-950a-49dd-8122-26e20e255996/GvsCreateVATfromVDS/54d7f7be-06d7-4bb9-bc65-72631185c10e/call-CreateCustomAnnotationsFile/mod_custom_annotations_template.tsv" .

        # cat ~{custom_annotations_body} | cut -f3- >> ~{custom_annotations_header}
        cat mod_custom_annotations_template.tsv > ~{custom_annotations_file_name}
        cat ~{custom_annotations_body} | cut -f3- | sed '1d'  >> ~{custom_annotations_file_name}


    >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:rc_616_var_store_2022_10_25"
        memory: "1 GB"
        cpu: "1"
        preemptible: 3
        disks: "local-disk 100 HDD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        File custom_annotations = "~{custom_annotations_file_name}"
    }
}


task AnnotateVCF {
    input {
        File input_vcf ## TODO do we need a sites only index file?
        String output_annotated_file_name
        File nirvana_data_tar
        File custom_annotations_file
    }
    String annotation_json_name = output_annotated_file_name + ".json.gz"
    String annotation_json_name_jsi = annotation_json_name + ".jsi"
    String nirvana_location = "/Nirvana/Nirvana.dll"
    String custom_creation_location = "/Nirvana/SAUtils.dll"
    String path = "/Cache/GRCh38/Both"
    String path_supplementary_annotations = "/SupplementaryAnnotation/GRCh38"
    String path_reference = "/References/Homo_sapiens.GRCh38.Nirvana.dat"

    command <<<
        # set -e

        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        # =======================================
        # Handle our data sources:

        echo "Extracting annotation data sources tar/gzip file..."
        mkdir datasources_dir
        tar zxvf ~{nirvana_data_tar} -C datasources_dir  ## --strip-components 2
        DATA_SOURCES_FOLDER="$PWD/datasources_dir/references"


        # =======================================
        echo "Creating custom annotations"
        mkdir customannotations_dir
        CUSTOM_ANNOTATIONS_FOLDER="$PWD/customannotations_dir"

        # Add AC/AN/AF as custom annotations
        ## use --skip-ref once you are on a version of nirvana later than 3.14 (once they have created a docker image for it)
        dotnet ~{custom_creation_location} customvar \
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
        docker: "us.gcr.io/broad-dsde-methods/variantstore:nirvana_2022_10_19"
        memory: "32 GB"
        cpu: "2"
        preemptible: 5
        disks: "local-disk 250 HDD"
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
    }

    String output_vt_json = "vat_vt_bq_load" + output_file_suffix
    String output_genes_json = "vat_genes_bq_load" + output_file_suffix
    String output_vt_gcp_path = output_path + 'vt/'
    String output_genes_gcp_path = output_path + 'genes/'
    String output_annotations_gcp_path = output_path + 'annotations/'

    ## note: these temp files do not currently get cleaned up as some of them may be helpful for recovery.

    command <<<
        # set -e

        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        # for debugging purposes only
        gsutil cp ~{annotation_json} '~{output_annotations_gcp_path}' #TODO: Can we run this task without cp-ing down this file?

        ## the annotation jsons are split into the specific VAT schema
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
        docker: "us.gcr.io/broad-dsde-methods/variantstore:rc_616_var_store_2022_10_25"
        memory: "60 GB"
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
    meta {
        # since the WDL will not see the updated data (its getting put in a gcp bucket)
        volatile: true
    }

    input {
        String filter_set_name
        String? vat_version
        File? nirvana_schema
        File? vt_schema
        File? genes_schema
        String project_id
        String dataset_name
        String output_path
        Boolean prep_jsons_done
    }

    # If the vat version is undefined or v1 then the vat tables would be named like filter_vat, otherwise filter_vat_v2.
    String effective_vat_version = if (defined(vat_version) && vat_version != "v1") then "_" + vat_version else ""

    # There are two pre-vat tables. A variant table and a genes table. They are joined together for the vat table
    String vat_table = filter_set_name + "_vat" + effective_vat_version
    String variant_transcript_table = filter_set_name + "_vat_vt" + effective_vat_version
    String genes_table = filter_set_name + "_vat_genes" + effective_vat_version

    String vt_path = output_path + 'vt/*'
    String genes_path = output_path + 'genes/*'

    command <<<
        echo "project_id = ~{project_id}" > ~/.bigqueryrc

        DATE=86400 ## 24 hours in seconds

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
        String vat_table_name = vat_table
        Boolean done = true
    }
}

task BigQueryExportVat {
    input {
        String contig
        String project_id
        String dataset_name
        String vat_table
        String output_path
        Boolean load_jsons_done
    }

    String export_path = output_path + "export/" + contig + "/*.tsv.gz"

    command <<<
        echo "project_id = ~{project_id}" > ~/.bigqueryrc

        # note: tab delimiter and compression creates tsv.gz files
        bq query --nouse_legacy_sql --project_id=~{project_id} \
        'EXPORT DATA OPTIONS(
        uri="~{export_path}",
        format="CSV",
        compression="GZIP",
        overwrite=true,
        header=false,
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

task MergeVatTSVs {
    input {
        Array[Boolean] export_done
        Array[String] contig_array
        String project_id
        String output_path

        Int? merge_vcfs_disk_size_override
    }

    # going large with the default to make gsutil -m cp really zippy
    Int disk_size = if (defined(merge_vcfs_disk_size_override)) then merge_vcfs_disk_size_override else 250

    command <<<
        apt-get update
        apt-get install tabix

        # custom function to prepend the current datetime to an echo statement "borrowed" from ExtractAnAcAfFromVCF
        echo_date () { echo "`date "+%Y/%m/%d %H:%M:%S"` $1"; }

        mkdir TSVs
        contigs=( ~{sep=' ' contig_array} )
        files="header.gz"

        echo_date "looping over contigs: $contigs"
        for i in "${contigs[@]}"
        do
            echo_date "copying files from ~{output_path}export/$i/*.tsv.gz"
            gsutil -m cp ~{output_path}export/$i/*.tsv.gz TSVs/
            echo_date "concatenating local tsv.gz files"

            # the internet says that * is deterministic, see https://serverfault.com/questions/122737/in-bash-are-wildcard-expansions-guaranteed-to-be-in-order
            cat TSVs/*.tsv.gz > vat_$i.tsv.gz

            echo_date "removing now concatenated files"
            rm TSVs/*.tsv.gz
            files="$files vat_$i.tsv.gz"
        done

        echo_date "making header.gz"
        echo "vid transcript contig position ref_allele alt_allele gvs_all_ac gvs_all_an gvs_all_af gvs_all_sc gvs_max_af gvs_max_ac gvs_max_an gvs_max_sc gvs_max_subpop gvs_afr_ac gvs_afr_an gvs_afr_af gvs_afr_sc gvs_amr_ac gvs_amr_an gvs_amr_af gvs_amr_sc gvs_eas_ac gvs_eas_an gvs_eas_af gvs_eas_sc gvs_eur_ac gvs_eur_an gvs_eur_af gvs_eur_sc gvs_mid_ac gvs_mid_an gvs_mid_af gvs_mid_sc gvs_oth_ac gvs_oth_an gvs_oth_af gvs_oth_sc gvs_sas_ac gvs_sas_an gvs_sas_af gvs_sas_sc gene_symbol transcript_source aa_change consequence dna_change_in_transcript variant_type exon_number intron_number genomic_location dbsnp_rsid gene_id gene_omim_id is_canonical_transcript gnomad_all_af gnomad_all_ac gnomad_all_an gnomad_failed_filter gnomad_max_af gnomad_max_ac gnomad_max_an gnomad_max_subpop gnomad_afr_ac gnomad_afr_an gnomad_afr_af gnomad_amr_ac gnomad_amr_an gnomad_amr_af gnomad_asj_ac gnomad_asj_an gnomad_asj_af gnomad_eas_ac gnomad_eas_an gnomad_eas_af gnomad_fin_ac gnomad_fin_an gnomad_fin_af gnomad_nfr_ac gnomad_nfr_an gnomad_nfr_af gnomad_sas_ac gnomad_sas_an gnomad_sas_af gnomad_oth_ac gnomad_oth_an gnomad_oth_af revel splice_ai_acceptor_gain_score splice_ai_acceptor_gain_distance splice_ai_acceptor_loss_score splice_ai_acceptor_loss_distance splice_ai_donor_gain_score splice_ai_donor_gain_distance splice_ai_donor_loss_score splice_ai_donor_loss_distance omim_phenotypes_id omim_phenotypes_name clinvar_classification clinvar_last_updated clinvar_phenotype" | gzip > header.gz

        echo_date "concatenating $files"
        cat $(echo $files) > vat_complete.tsv.gz
        echo_date "bgzipping concatenated file"
        cat vat_complete.tsv.gz | gunzip | bgzip > vat_complete.bgz.tsv.gz
        echo_date "copying bgzipped file to ~{output_path}"
        gsutil -m cp vat_complete.bgz.tsv.gz ~{output_path}
    >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/broad-gatk/gatk:4.2.6.1"
        memory: "2 GB"
        preemptible: 3
        cpu: "1"
        disks: "local-disk ~{disk_size} HDD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        File tsv_file = "vat_complete.bgz.tsv.gz"
    }
}

