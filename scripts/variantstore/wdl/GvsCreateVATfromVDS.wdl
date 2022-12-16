version 1.0

import "GvsUtils.wdl" as Utils

workflow GvsCreateVATfromVDS {
    input {
        File input_sites_only_vcf
        File ancestry_file

        String project_id
        String dataset_name
        String filter_set_name
        String? vat_version

        Int effective_scatter_count = 10

        String output_path
        Int? split_intervals_disk_size_override
        Int? split_intervals_mem_override
        Int? merge_vcfs_disk_size_override
    }

    File interval_list = "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list"

    Array[String] contig_array = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"]
    File reference = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
    File reference_dict = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict"
    File reference_index = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"

    File nirvana_data_directory = "gs://gvs_quickstart_storage/Nirvana/Nirvana-references-2022-10-07.tgz"

    ## Flag E

    call MakeSubpopulationFilesAndReadSchemaFiles {
        input:
            input_ancestry_file = ancestry_file
    }

    call RemoveDuplicatesFromSitesOnlyVCF {
        input:
            sites_only_vcf = input_sites_only_vcf,
            ref = reference
    }

    call Utils.IndexVcf {
        input:
            input_vcf = RemoveDuplicatesFromSitesOnlyVCF.output_vcf
    }

    call Utils.SplitIntervals {
        input:
            intervals = interval_list,
            ref_fasta = reference,
            ref_fai = reference_index,
            ref_dict = reference_dict,
            scatter_count = effective_scatter_count,
            output_gcs_dir = output_path + "intervals",
            split_intervals_disk_size_override = split_intervals_disk_size_override,
            split_intervals_mem_override = split_intervals_mem_override,
    }

    String sites_only_vcf_basename = basename(basename(input_sites_only_vcf, ".gz"), ".vcf")

    scatter(i in range(length(SplitIntervals.interval_files))) {
        String interval_file_basename = basename(SplitIntervals.interval_files[i], ".interval_list")
        String vcf_filename = interval_file_basename + "." + sites_only_vcf_basename

        call Utils.SelectVariants {
            input:
                input_vcf = IndexVcf.output_vcf,
                input_vcf_index = IndexVcf.output_vcf_index,
                interval_list = SplitIntervals.interval_files[i],
                output_basename = vcf_filename
        }

        call StripCustomAnnotationsFromSitesOnlyVCF {
            input:
                input_vcf = SelectVariants.output_vcf,
                custom_annotations_header = MakeSubpopulationFilesAndReadSchemaFiles.custom_annotations_template_file,
                output_vcf_name = "${vcf_filename}.unannotated.sites_only.vcf",
                output_custom_annotations_filename = "${vcf_filename}.custom_annotations.tsv"
        }

        ## Use Nirvana to annotate the sites-only VCF and include the AC/AN/AF calculations as custom annotations
        call AnnotateVCF {
            input:
                input_vcf = StripCustomAnnotationsFromSitesOnlyVCF.output_vcf,
                output_annotated_file_name = "${vcf_filename}_annotated",
                nirvana_data_tar = nirvana_data_directory,
                custom_annotations_file = StripCustomAnnotationsFromSitesOnlyVCF.output_custom_annotations_file,
        }


        call PrepVtAnnotationJson {
            input:
                positions_annotation_json = AnnotateVCF.positions_annotation_json,
                output_file_suffix = "${vcf_filename}.json.gz",
                output_path = output_path,
        }

        call PrepGenesAnnotationJson {
            input:
                genes_annotation_json = AnnotateVCF.genes_annotation_json,
                output_file_suffix = "${vcf_filename}.json.gz",
                output_path = output_path,
        }

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
            prep_vt_json_done = PrepVtAnnotationJson.done,
            prep_genes_json_done = PrepGenesAnnotationJson.done
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
        docker: "us.gcr.io/broad-dsde-methods/variantstore:gg_VS-757_var_store_2022_12_15d"
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


task StripCustomAnnotationsFromSitesOnlyVCF {
    input {
        File input_vcf
        File custom_annotations_header
        String output_vcf_name
        String output_custom_annotations_filename
        File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"
    }

    Int disk_size = ceil((size(input_vcf, "GB") + size(custom_annotations_header, "GB")) * 4) + 100

    command <<<
        set -o errexit -o nounset -o pipefail -o xtrace

        bash ~{monitoring_script} > monitoring.log &

        python3 /app/strip_custom_annotations_from_sites_only_vcf.py \
        --input_vcf ~{input_vcf} \
        --input_custom_annotations_tsv ~{custom_annotations_header} \
        --output_vcf ~{output_vcf_name} \
        --output_custom_annotations_tsv ~{output_custom_annotations_filename}

    >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:gg_VS-757_var_store_2022_12_15d"
        memory: "7 GiB"
        cpu: "2"
        preemptible: 3
        disks: "local-disk " + disk_size + " HDD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        File output_vcf = output_vcf_name
        File output_custom_annotations_file = output_custom_annotations_filename
        File monitoring_log = "monitoring.log"
    }
}


task RemoveDuplicatesFromSitesOnlyVCF {
    input {
        File sites_only_vcf
        File ref
        File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"
    }

    Int disk_size = ceil(size(sites_only_vcf, "GB") * 5) + 100

    # separate multi-allelic sites into their own lines, remove deletions and filtered sites and make a sites only vcf
    # while extracting and calculating the an/ac/af & sc by subpopulation into a tsv
    command <<<
        set -e

        bash ~{monitoring_script} > monitoring.log &

        # custom function to prepend the current datetime to an echo statement
        echo_date () { echo "`date "+%Y/%m/%d %H:%M:%S"` $1"; }

        echo_date "VAT: Convert input to BCF format"
        bcftools convert --threads 4 -O b -o sites_only.bcf ~{sites_only_vcf}

        echo_date "VAT: Calculating number of sites with Ns"

        ## track the dropped variants with N's in the reference (Since Nirvana cant handle N as a base, drop them for now)
        bcftools view --threads 4 -i 'REF~"N"' -O u sites_only.bcf | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' > track_dropped.tsv

        echo_date "VAT: filter out sites with N's in the reference AND sites with AC=0"
        # TODO - NOTE - we are NOT tracking the sites with AC=0 - should we? For Lee (Rori is checking)
        bcftools view --threads 4 -e 'REF~"N" || AC=0' -O b sites_only.bcf -o filtered_sites_only.bcf
        rm sites_only.bcf

        echo_date "VAT: normalize, left align and split multi allelic sites to new lines, remove duplicate lines"
        ## note that normalization may create sites with more than 50 alt alleles
        bcftools norm --threads 4 -m- --check-ref w -f ~{ref} filtered_sites_only.bcf -O b -o normalized.bcf
        rm filtered_sites_only.bcf

        echo_date "VAT: detecting and removing duplicate rows from sites-only VCF"

        ## During normalization, sometimes duplicate variants appear but with different calculations. This seems to be a bug in bcftools. For now we are dropping all duplicate variants
        ## to locate the duplicates, we first make a file of just the first 5 columns
        bcftools query normalized.bcf -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' | sort | uniq -d > duplicates.tsv

        echo_date "VAT: done with duplicate detection"
        wc -l duplicates.tsv
        echo_date "VAT: Duplicates may have been found"

        # If there ARE dupes to remove
        if [ -s duplicates.tsv ]; then
            ## remove those rows (that match up to the first 5 cols)
            echo_date "VAT: Removing those rows"
            bcftools view --threads 4 normalized.bcf | grep -v -wFf duplicates.tsv > deduplicated.vcf
        else
            # There are no duplicates to remove
            echo_date "VAT: No duplicates found"
            bcftools view --threads 4 normalized.bcf -o deduplicated.vcf
        fi
        rm normalized.bcf

        ## add duplicates to the file that's tracking dropped variants
        cat duplicates.tsv >> track_dropped.tsv
        rm duplicates.tsv ## clean up unneeded file

        echo_date "VAT: finished"
    >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:2022-10-25-alpine"
        maxRetries: 3
        memory: "16 GB"
        preemptible: 3
        cpu: "8"
        disks: "local-disk " + disk_size + " HDD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        File track_dropped = "track_dropped.tsv"
        File output_vcf = "deduplicated.vcf"
        File monitoring_log = "monitoring.log"
    }
}


task AnnotateVCF {
    input {
        File input_vcf
        String output_annotated_file_name
        File nirvana_data_tar
        File custom_annotations_file
        File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"
    }
    String annotation_json_name = output_annotated_file_name + ".json.gz"
    String gene_annotation_json_name = output_annotated_file_name + ".genes.json.gz"
    String positions_annotation_json_name = output_annotated_file_name + ".positions.json.gz"
    String nirvana_location = "/Nirvana/Nirvana.dll"
    String custom_creation_location = "/Nirvana/SAUtils.dll"
    String jasix_location = "/Nirvana/Jasix.dll"
    String path = "/Cache/GRCh38/Both"
    String path_supplementary_annotations = "/SupplementaryAnnotation/GRCh38"
    String path_reference = "/References/Homo_sapiens.GRCh38.Nirvana.dat"

    command <<<
        # set -e

        bash ~{monitoring_script} > monitoring.log &

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
            -i ~{input_vcf} \
            -c $DATA_SOURCES_FOLDER~{path} \
            --sd $DATA_SOURCES_FOLDER~{path_supplementary_annotations} \
            --sd $CUSTOM_ANNOTATIONS_FOLDER \
            -r $DATA_SOURCES_FOLDER~{path_reference} \
            -o ~{output_annotated_file_name}

        # https://illumina.github.io/NirvanaDocumentation/introduction/parsing-json#jasix
        # Parse out the Genes section into a separate annotated json
        dotnet  ~{jasix_location} \
            --in ~{annotation_json_name} \
            --section genes \
            --out ~{gene_annotation_json_name}

        # Parse out the Positions section into a separate annotated json
        dotnet  ~{jasix_location} \
        --in ~{annotation_json_name} \
        --section positions \
        --out ~{positions_annotation_json_name}

    >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:nirvana_2022_10_19"
        memory: "64 GB"
        cpu: "4"
        preemptible: 3
        disks: "local-disk 2000 HDD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        File genes_annotation_json = "~{gene_annotation_json_name}"
        File positions_annotation_json = "~{positions_annotation_json_name}"
        File monitoring_log = "monitoring.log"
    }
}

task PrepVtAnnotationJson {
    input {
        File positions_annotation_json
        String output_file_suffix
        String output_path
    }

    File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"

    String output_vt_json = "vat_vt_bq_load" + output_file_suffix
    String output_vt_gcp_path = output_path + 'vt/'

    command <<<
        set -o errexit -o nounset -o pipefail -o xtrace

        # Kick off the monitoring script
        bash ~{monitoring_script} > monitoring.log &

        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '

        ## the annotation jsons are split into the specific VAT schema
        python3 /app/create_vt_bqloadjson_from_annotations.py \
            --annotated_json ~{positions_annotation_json} \
            --output_vt_json ~{output_vt_json}

        gsutil cp ~{output_vt_json} '~{output_vt_gcp_path}'

    >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:gg_VS-757_var_store_2022_12_15d"
        memory: "7 GB"
        preemptible: 3
        cpu: "1"
        disks: "local-disk 500 HDD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        File vat_vt_json="~{output_vt_json}"
        Boolean done = true
        File monitoring_log = "monitoring.log"
    }
}

task PrepGenesAnnotationJson {
    input {
        File genes_annotation_json
        String output_file_suffix
        String output_path
    }

    # Kick off the monitoring script
    File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"

    String output_genes_json = "vat_genes_bq_load" + output_file_suffix
    String output_genes_gcp_path = output_path + 'genes/'

    command <<<
        set -o errexit -o nounset -o pipefail -o xtrace

        bash ~{monitoring_script} > monitoring.log &

        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '

        ## the annotation jsons are split into the specific VAT schema
        python3 /app/create_genes_bqloadjson_from_annotations.py \
            --annotated_json ~{genes_annotation_json} \
            --output_genes_json ~{output_genes_json}

        gsutil cp ~{output_genes_json} '~{output_genes_gcp_path}'

    >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:gg_VS-757_var_store_2022_12_15d"
        memory: "7 GB"
        preemptible: 3
        cpu: "1"
        disks: "local-disk 500 HDD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        File vat_genes_json="~{output_genes_json}"
        Boolean done = true
        File monitoring_log = "monitoring.log"
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
        Array[Boolean] prep_vt_json_done
        Array[Boolean] prep_genes_json_done
    }

    # If the vat version is undefined or v1 then the vat tables would be named like filter_vat, otherwise filter_vat_v2.
    String effective_vat_version = if (defined(vat_version) && select_first([vat_version]) != "v1") then "_" + select_first([vat_version]) else ""

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
        String output_path

        Int? merge_vcfs_disk_size_override
    }

    # going large with the default to make gsutil -m cp really zippy
    Int disk_size = if (defined(merge_vcfs_disk_size_override)) then select_first([merge_vcfs_disk_size_override]) else 250

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
