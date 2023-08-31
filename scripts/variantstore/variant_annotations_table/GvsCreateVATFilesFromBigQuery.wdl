version 1.0

import "../wdl/GvsUtils.wdl" as Utils

workflow GvsCreateVATFilesFromBigQuery {
    input {
        String project_id
        String? git_branch_or_tag
        String? git_hash
        String dataset_name
        String vat_table_name

        String output_path
        Int? merge_vcfs_disk_size_override
        Boolean precondition_met = true
        String? cloud_sdk_docker
        String? cloud_sdk_slim_docker  # TODO - probably unused.
    }

    Array[String] contig_array = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"]
    String small_contig_for_header = "chrY"

    if (!defined(git_hash) || !defined(cloud_sdk_docker) || !defined(cloud_sdk_slim_docker)) {
        call Utils.GetToolVersions {
            input:
                git_branch_or_tag = git_branch_or_tag,
        }
    }

    String effective_cloud_sdk_docker = select_first([cloud_sdk_docker, GetToolVersions.cloud_sdk_docker])
    String effective_cloud_sdk_slim_docker = select_first([cloud_sdk_slim_docker, GetToolVersions.cloud_sdk_slim_docker])
    String effective_git_hash = select_first([git_hash, GetToolVersions.git_hash])

    call BigQueryExportVat as BigQueryExportVatForHeader {
        input:
            contig = small_contig_for_header,
            project_id = project_id,
            dataset_name = dataset_name,
            output_path = output_path + "export_for_header/",
            add_header = true,
            vat_table = vat_table_name,
            load_jsons_done = precondition_met,
            cloud_sdk_docker = effective_cloud_sdk_docker,
    }

    scatter(i in range(length(contig_array)) ) {
        call BigQueryExportVat {
            input:
                contig = contig_array[i],
                project_id = project_id,
                dataset_name = dataset_name,
                output_path = output_path + "export/",
                vat_table = vat_table_name,
                load_jsons_done = precondition_met,
                cloud_sdk_docker = effective_cloud_sdk_docker,
        }
    }

    call GetHeader {
        input:
            export_done = BigQueryExportVatForHeader.done,
            header_contig = small_contig_for_header,
            output_path = output_path + "export_for_header/",
            merge_vcfs_disk_size_override = merge_vcfs_disk_size_override,
            cloud_sdk_docker = effective_cloud_sdk_docker,
    }

    call MergeVatTSVs {
        input:
            export_done = BigQueryExportVat.done,
            contig_array = contig_array,
            header_tsv_file = GetHeader.header_tsv_file,
            output_path = output_path + "export/",
            merge_vcfs_disk_size_override = merge_vcfs_disk_size_override,
            cloud_sdk_docker = effective_cloud_sdk_docker,
    }

    output {
        File final_tsv_file = MergeVatTSVs.tsv_file
        String recorded_git_hash = effective_git_hash
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
        Boolean add_header = false
        Boolean load_jsons_done
        String cloud_sdk_docker
    }

    String export_path = output_path + contig + "/*.tsv.gz"

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
        header=~{add_header},
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
        docker: cloud_sdk_docker
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

task GetHeader {
    input {
        Boolean export_done
        String header_contig
        String output_path

        Int? merge_vcfs_disk_size_override
        String cloud_sdk_docker
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

        # custom function to prepend the current datetime to an echo statement "borrowed" from ExtractAnAcAfFromVCF
        echo_date () { echo "`date "+%Y/%m/%d %H:%M:%S"` $1"; }

        mkdir TSVs
        echo_date "copying files from ~{output_path}~{header_contig}/*.tsv.gz"
        gcloud storage cp ~{output_path}~{header_contig}/*.tsv.gz TSVs/
        tsv_files=(TSVs/*.tsv.gz)

        # Get the first file and strip off its header.
        # We can't safely pipe to `head -1` because while `head` will exit successfully after reading the first line, the
        # pipeline will continue trying to write data to the `head` process. If this happens we'll get a 141 exit code and
        # with `set -o pipefail` turned on this will fail our task. As a workaround use this `<(...)` temp file construct.
        # https://news.ycombinator.com/item?id=9255830
        head -1 <(gzip -cd ${tsv_files[0]}) | gzip > header.tsv.gz
    >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: cloud_sdk_docker
        memory: "4 GB"
        preemptible: 3
        cpu: "2"
        disks: "local-disk ~{disk_size} HDD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        File header_tsv_file = "header.tsv.gz"
        File monitoring_log = "monitoring.log"
    }
}

task MergeVatTSVs {
    input {
        Array[Boolean] export_done
        Array[String] contig_array
        File header_tsv_file
        String output_path

        Int? merge_vcfs_disk_size_override
        String cloud_sdk_docker
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

        # custom function to prepend the current datetime to an echo statement "borrowed" from ExtractAnAcAfFromVCF
        echo_date () { echo "`date "+%Y/%m/%d %H:%M:%S"` $1"; }

        mkdir TSVs
        contigs=( ~{sep=' ' contig_array} )
        cp ~{header_tsv_file} header.gz
        files="header.gz"

        echo_date "looping over contigs: $contigs"
        for i in "${contigs[@]}"
        do
            echo_date "copying files from ~{output_path}$i/*.tsv.gz"
            gcloud storage cp ~{output_path}$i/*.tsv.gz TSVs/
            echo_date "concatenating local tsv.gz files"

            # the internet says that * is deterministic, see https://serverfault.com/questions/122737/in-bash-are-wildcard-expansions-guaranteed-to-be-in-order
            cat TSVs/*.tsv.gz > vat_$i.tsv.gz

            echo_date "removing now concatenated files"
            rm TSVs/*.tsv.gz
            files="$files vat_$i.tsv.gz"
        done

        echo_date "concatenating $files"
        cat $(echo $files) > vat_complete.tsv.gz
        echo_date "copying output file to ~{output_path}"
        gcloud storage cp vat_complete.tsv.gz ~{output_path}
    >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: cloud_sdk_docker
        memory: "4 GB"
        preemptible: 3
        cpu: "2"
        disks: "local-disk ~{disk_size} HDD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        File tsv_file = "vat_complete.tsv.gz"
        File monitoring_log = "monitoring.log"
    }
}
