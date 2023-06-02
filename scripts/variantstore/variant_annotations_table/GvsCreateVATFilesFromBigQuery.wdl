version 1.0

workflow GvsCreateVATFilesFromBigQuery {
    input {
        String project_id
        String dataset_name
        String vat_table_name

        String output_path
        Int? merge_vcfs_disk_size_override
        Boolean precondition_met = true
    }

    Array[String] contig_array = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"]

    scatter(i in range(length(contig_array)) ) {
        call BigQueryExportVat {
            input:
                contig = contig_array[i],
                project_id = project_id,
                dataset_name = dataset_name,
                output_path = output_path,
                vat_table = vat_table_name,
                load_jsons_done = precondition_met
        }
    }

    call MergeVatTSVs {
        input:
            export_done = BigQueryExportVat.done,
            contig_array = contig_array,
            output_path = output_path,
            merge_vcfs_disk_size_override = merge_vcfs_disk_size_override
    }

    output {
        File final_tsv_file = MergeVatTSVs.tsv_file
    }
}

################################################################################


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
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace
        
        echo "project_id = ~{project_id}" > ~/.bigqueryrc

        # note: tab delimiter and compression creates tsv.gz files
        bq --apilog=false query --nouse_legacy_sql --project_id=~{project_id} \
        'EXPORT DATA OPTIONS(
        uri="~{export_path}",
        format="CSV",
        compression="GZIP",
        overwrite=true,
        header=false,
        field_delimiter="\t") AS
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
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:426.0.0-alpine"
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
        String output_path

        Int? merge_vcfs_disk_size_override
    }

    File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"

    # going large with the default to make gsutil -m cp really zippy
    Int disk_size = if (defined(merge_vcfs_disk_size_override)) then select_first([merge_vcfs_disk_size_override]) else 500

    command <<<
        # Kick off the monitoring script
        bash ~{monitoring_script} > monitoring.log &
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace
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
            gcloud storage cp ~{output_path}export/$i/*.tsv.gz TSVs/
            echo_date "concatenating local tsv.gz files"

            # the internet says that * is deterministic, see https://serverfault.com/questions/122737/in-bash-are-wildcard-expansions-guaranteed-to-be-in-order
            cat TSVs/*.tsv.gz > vat_$i.tsv.gz

            echo_date "removing now concatenated files"
            rm TSVs/*.tsv.gz
            files="$files vat_$i.tsv.gz"
        done

        echo_date "making header.gz"
        # NOTE: Contents of tsvs exported from BigQuery are tab-separated, the header must also be tab-separated!
        echo -e "vid\ttranscript\tcontig\tposition\tref_allele\talt_allele\tgvs_all_ac\tgvs_all_an\tgvs_all_af\tgvs_all_sc\tgvs_max_af\tgvs_max_ac\tgvs_max_an\tgvs_max_sc\tgvs_max_subpop\tgvs_afr_ac\tgvs_afr_an\tgvs_afr_af\tgvs_afr_sc\tgvs_amr_ac\tgvs_amr_an\tgvs_amr_af\tgvs_amr_sc\tgvs_eas_ac\tgvs_eas_an\tgvs_eas_af\tgvs_eas_sc\tgvs_eur_ac\tgvs_eur_an\tgvs_eur_af\tgvs_eur_sc\tgvs_mid_ac\tgvs_mid_an\tgvs_mid_af\tgvs_mid_sc\tgvs_oth_ac\tgvs_oth_an\tgvs_oth_af\tgvs_oth_sc\tgvs_sas_ac\tgvs_sas_an\tgvs_sas_af\tgvs_sas_sc\tgene_symbol\ttranscript_source\taa_change\tconsequence\tdna_change_in_transcript\tvariant_type\texon_number\tintron_number\tgenomic_location\tdbsnp_rsid\tgene_id\tgene_omim_id\tis_canonical_transcript\tgnomad_all_af\tgnomad_all_ac\tgnomad_all_an\tgnomad_failed_filter\tgnomad_max_af\tgnomad_max_ac\tgnomad_max_an\tgnomad_max_subpop\tgnomad_afr_ac\tgnomad_afr_an\tgnomad_afr_af\tgnomad_amr_ac\tgnomad_amr_an\tgnomad_amr_af\tgnomad_asj_ac\tgnomad_asj_an\tgnomad_asj_af\tgnomad_eas_ac\tgnomad_eas_an\tgnomad_eas_af\tgnomad_fin_ac\tgnomad_fin_an\tgnomad_fin_af\tgnomad_nfr_ac\tgnomad_nfr_an\tgnomad_nfr_af\tgnomad_sas_ac\tgnomad_sas_an\tgnomad_sas_af\tgnomad_oth_ac\tgnomad_oth_an\tgnomad_oth_af\trevel\tsplice_ai_acceptor_gain_score\tsplice_ai_acceptor_gain_distance\tsplice_ai_acceptor_loss_score\tsplice_ai_acceptor_loss_distance\tsplice_ai_donor_gain_score\tsplice_ai_donor_gain_distance\tsplice_ai_donor_loss_score\tsplice_ai_donor_loss_distance\tomim_phenotypes_id\tomim_phenotypes_name\tclinvar_classification\tclinvar_last_updated\tclinvar_phenotype" | gzip > header.gz

        echo_date "concatenating $files"
        cat $(echo $files) > vat_complete.tsv.gz
        echo_date "bgzipping concatenated file"
        cat vat_complete.tsv.gz | gunzip | bgzip > vat_complete.bgz.tsv.gz
        echo_date "copying bgzipped file to ~{output_path}"
        gcloud storage cp vat_complete.bgz.tsv.gz ~{output_path}
    >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:426.0.0-slim"
        memory: "4 GB"
        preemptible: 3
        cpu: "2"
        disks: "local-disk ~{disk_size} HDD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        File tsv_file = "vat_complete.bgz.tsv.gz"
        File monitoring_log = "monitoring.log"
    }
}
