version 1.0

import "../wdl/GvsUtils.wdl" as Utils

workflow GvsAnnotateAnyVcf {
    input {
        String project_id
        String dataset_name
        String reference_name
        File sites_only_vcf
        File sites_only_vcf_index
        String output_path

        String? variants_nirvana_docker
        String? git_branch_or_tag
        Boolean use_reference_disk = true
    }

    parameter_meta {
        project_id: {
            help: "Google project ID for the GVS BigQuery dataset"
        }
        dataset_name: {
            help: "BigQuery dataset name for GVS"
        }
        sites_only_vcf: {
            help: "Sites-only VCF file for annotation"
        }
        output_path: {
            help: "GCS location (with a trailing '/') to put temporary and output files for the VAT pipeline"
        }

    }

    call Utils.GetToolVersions {
        input:
            git_branch_or_tag = git_branch_or_tag,
    }

    String effective_basic_docker = GetToolVersions.basic_docker
    String effective_cloud_sdk_docker = GetToolVersions.cloud_sdk_docker
    String effective_cloud_sdk_slim_docker = GetToolVersions.cloud_sdk_slim_docker
    String effective_variants_docker = GetToolVersions.variants_docker
    String effective_gatk_docker = GetToolVersions.gatk_docker
    String effective_variants_nirvana_docker = GetToolVersions.variants_nirvana_docker
    String effective_google_project = GetToolVersions.google_project

    String output_path_without_a_trailing_slash = sub(output_path, "/$", "")
    String effective_output_path = if (output_path == output_path_without_a_trailing_slash) then output_path + "/" else output_path


    if (!defined(sites_only_vcf)) {
        call Utils.TerminateWorkflow as MustSetSitesOnlyVcfCreationParameters {
            input:
                basic_docker = effective_basic_docker,
        }
    }

    call Utils.GetReference {
        input:
            reference_name = reference_name,
            basic_docker = effective_basic_docker,
    }

    String interval_list = GetReference.reference.wgs_calling_interval_list
    String reference_fasta = GetReference.reference.reference_fasta

    if (defined(sites_only_vcf)) {

        call Utils.CopyFile as CopySitesOnlyVcf {
            input:
                input_file = sites_only_vcf,
                output_gcs_dir = effective_output_path + "sites_only_vcf",
                cloud_sdk_docker = effective_cloud_sdk_docker,
        }

        if (!defined(sites_only_vcf_index)) {
            call Utils.IndexVcf {
                input:
                    input_vcf = CopySitesOnlyVcf.output_file_path,
                    gatk_docker = effective_gatk_docker,
            }
        }

        call Utils.CopyFile as CopySitesOnlyVcfIndex {
            input:
                input_file = select_first([sites_only_vcf_index, IndexVcf.output_vcf_index]),
                output_gcs_dir = effective_output_path + "sites_only_vcf",
                cloud_sdk_docker = effective_cloud_sdk_docker,
        }

#        call Utils.SplitIntervals { ## if I dont know how many samples---I prob dont want to split intervals
#            input:
#                intervals = interval_list,
#                ref_fasta = reference_fasta,
#                scatter_count = effective_scatter_count,
#                output_gcs_dir = effective_output_path + "intervals",
#                split_intervals_disk_size_override = split_intervals_disk_size_override,
#                split_intervals_mem_override = split_intervals_mem_override,
#                gatk_docker = effective_gatk_docker,
#        }

        String sites_only_vcf_basename = basename(CopySitesOnlyVcf.output_file_path, ".sites-only.vcf")

#        scatter(i in range(length(SplitIntervals.interval_files))) {
#            String interval_file_basename = basename(SplitIntervals.interval_files[i], ".interval_list")
#            String vcf_filename = interval_file_basename + "." + sites_only_vcf_basename

#            call Utils.SelectVariants {
#                input:
#                    input_vcf = CopySitesOnlyVcf.output_file_path,
#                    input_vcf_index = CopySitesOnlyVcfIndex.output_file_path,
#                    interval_list = SplitIntervals.interval_files[i],
#                    output_basename = vcf_filename,
#                    gatk_docker = effective_gatk_docker,
#            }

#            call RemoveDuplicatesFromSitesOnlyVCF {
#                input:
#                    sites_only_vcf = SelectVariants.output_vcf,
#                    ref = reference_fasta,
#                    variants_docker = effective_variants_docker,
#            }

            ## Use Nirvana to annotate the sites-only VCF and include the AC/AN/AF calculations as custom annotations
            call AnnotateVCF {
                input:
                    cromwell_root = GetToolVersions.cromwell_root,
                    input_vcf = CopySitesOnlyVcf.output_file_path,
                    output_annotated_file_name = "${vcf_filename}_annotated",
                    variants_nirvana_docker = effective_variants_nirvana_docker,
                    use_reference_disk = use_reference_disk,
            }

#            call PrepVtAnnotationJson {
#                input:
#                    positions_annotation_json = AnnotateVCF.positions_annotation_json,
#                    output_file_suffix = "${vcf_filename}.json.gz",
#                    output_path = effective_output_path,
#                    variants_docker = effective_variants_docker,
#            }
#
#            call PrepGenesAnnotationJson {
#                input:
#                    genes_annotation_json = AnnotateVCF.genes_annotation_json,
#                    output_file_suffix = "${vcf_filename}.json.gz",
#                    output_path = effective_output_path,
#                    variants_docker = effective_variants_docker,
#            }

#        }
    }

    output {
        File genes_annotation_json = AnnotateVCF.genes_annotation_json
        File positions_annotation_json = AnnotateVCF.positions_annotation_json
        File monitoring_log = AnnotateVCF.monitoring_log
    }
}



task RemoveDuplicatesFromSitesOnlyVCF {
    input {
        File sites_only_vcf
        File ref
        String variants_docker
    }

    File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"
    Int disk_size = ceil(size(sites_only_vcf, "GB") * 5) + 100

    # separate multi-allelic sites into their own lines, remove deletions and filtered sites and make a sites only vcf
    # while extracting and calculating the an/ac/af & sc by subpopulation into a tsv
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        bash ~{monitoring_script} > monitoring.log &

        # custom function to prepend the current datetime to an echo statement
        echo_date () { echo "`date "+%Y/%m/%d %H:%M:%S"` $1"; }

        echo_date "VAT: Convert input to BCF format"
        bcftools convert --threads 4 -O b -o sites_only.bcf ~{sites_only_vcf}

        echo_date "VAT: Calculating number of sites with Ns"

        ## track the dropped variants with N's in the reference (Since Nirvana cant handle N as a base, drop them for now)
        bcftools view --threads 4 -i 'REF~"N"' -O u sites_only.bcf | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' > track_dropped.tsv

        echo_date "VAT: filter out sites with N's in the reference AND sites with AC=0"
        ## NOTE: Sites that were filtered out because of AC=0 are not recorded in the 'track_dropped.tsv' file, but can be
        ##       determined by examining the sites-only VCF provided to this WDL.
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

    runtime {
        docker: variants_docker
        maxRetries: 3
        memory: "16 GB"
        preemptible: 3
        cpu: "8"
        disks: "local-disk " + disk_size + " HDD"
    }

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
        String cromwell_root

        # Mentioning this path in the inputs section of the task combined with checking the 'Use reference disks' option
        # in Terra UI tells Cromwell to arrange for the Nirvana reference disk to be attached to this VM.
        File summon_reference_disk =
            "gs://gcp-public-data--broad-references/hg38/v0/Nirvana/3.18.1_2024-03-06/SupplementaryAnnotation/GRCh38/MITOMAP_20200819.nsa.idx"

        String variants_nirvana_docker

        File omim_annotations = "gs://gcp-public-data--broad-references/hg38/v0/Nirvana/3.18.1_2024-03-06/SupplementaryAnnotation/GRCh38/OMIM_20220516.nga"
        File cosmic_gene_fusion_annotations = "gs://gcp-public-data--broad-references/hg38/v0/Nirvana/3.18.1_2024-03-06/SupplementaryAnnotation/GRCh38/COSMIC_GeneFusions_94.gfj"
        File primate_ai_annotations = "gs://gcp-public-data--broad-references/hg38/v0/Nirvana/3.18.1_2024-03-06/SupplementaryAnnotation/GRCh38/PrimateAI_0.2.nsa"
        File primate_ai_annotations_idx = "gs://gcp-public-data--broad-references/hg38/v0/Nirvana/3.18.1_2024-03-06/SupplementaryAnnotation/GRCh38/PrimateAI_0.2.nsa.idx"
        File splice_ai_annotations = "gs://gcp-public-data--broad-references/hg38/v0/Nirvana/3.18.1_2024-03-06/SupplementaryAnnotation/GRCh38/SpliceAi_1.3.nsa"
        File splice_ai_annotations_idx = "gs://gcp-public-data--broad-references/hg38/v0/Nirvana/3.18.1_2024-03-06/SupplementaryAnnotation/GRCh38/SpliceAi_1.3.nsa.idx"
        Boolean use_reference_disk
    }

    File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"

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
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        bash ~{monitoring_script} > monitoring.log &

        if [[ "~{use_reference_disk}" == "true" ]]
        then
            # There's an issue with how the projects/broad-dsde-cromwell-dev/global/images/nirvana-3-18-1-references-2023-01-03
            # disk image was built: while all the reference files do exist on the image they are not at the expected
            # locations. The following code works around this issue and should continue to work even after a corrected
            # version of the Nirvana reference image is deployed into Terra.

            # Find where the reference disk should have been mounted on this VM.  Note this is referred to as a "candidate
            # mount point" because we do not actually confirm this is a reference disk until the following code block.
            if [[ -e /cromwell_root/gcs_delocalization.sh ]]
            then
              # PAPI mount points
              CANDIDATE_MOUNT_POINT=$(lsblk | grep -v cromwell_root | sed -E -n 's!.*(/mnt/[a-f0-9]+).*!\1!p')
            elif [[ -e /mnt/disks/cromwell_root/gcs_delocalization.sh ]]
            then
              # GCP Batch mount points
              CANDIDATE_MOUNT_POINT=$(lsblk | grep -v cromwell_root | sed -E -n 's!.*(/mnt/disks/[a-f0-9]+).*!\1!p')
            else
              >&2 echo "Could not find a mounted volume that looks like a reference disk, exiting."
              exit 1
            fi

            # Find one particular reference under the mount path. Note this is not the same reference as was specified in the
            # `inputs` section, so this would only be present if the volume we're looking at is in fact a reference disk.
            REFERENCE_FILE="Homo_sapiens.GRCh38.Nirvana.dat"
            REFERENCE_PATH=$(find ${CANDIDATE_MOUNT_POINT} -name "${REFERENCE_FILE}")
            if [[ -z ${REFERENCE_PATH} ]]; then
                >&2 echo "Could not find reference file '${REFERENCE_FILE}' under candidate reference disk mount point '${CANDIDATE_MOUNT_POINT}', exiting."
                exit 1
            fi

            # Take the parent of the parent directory of this file as root of the locally mounted references:
            DATA_SOURCES_FOLDER="$(dirname $(dirname ${REFERENCE_PATH}))"
        else
            DATA_SOURCES_FOLDER=~{cromwell_root}/nirvana_references
            mkdir ${DATA_SOURCES_FOLDER}

            # Download the references
            dotnet /Nirvana/Downloader.dll --ga GRCh38 --out ${DATA_SOURCES_FOLDER}

            # As of 2024-01-24 OMIM is no longer included among the bundle of annotation resources pulled down by the
            # Nirvana downloader. As this annotation set is currently central for our VAT logic, special-case link in
            # the OMIM .nsa bundle we downloaded back when we made the Delta reference disk:
            ln ~{omim_annotations} ${DATA_SOURCES_FOLDER}/SupplementaryAnnotation/GRCh38/
            # Similarly, the following annotations were removed from the latest Nirvana annotations (3.18.1), but we
            # re-add them as desired by Lee
            ln ~{cosmic_gene_fusion_annotations} ${DATA_SOURCES_FOLDER}/SupplementaryAnnotation/GRCh38/
            ln ~{primate_ai_annotations} ${DATA_SOURCES_FOLDER}/SupplementaryAnnotation/GRCh38/
            ln ~{primate_ai_annotations_idx} ${DATA_SOURCES_FOLDER}/SupplementaryAnnotation/GRCh38/
            ln ~{splice_ai_annotations} ${DATA_SOURCES_FOLDER}/SupplementaryAnnotation/GRCh38/
            ln ~{splice_ai_annotations_idx} ${DATA_SOURCES_FOLDER}/SupplementaryAnnotation/GRCh38/
        fi

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

    runtime {
        docker: variants_nirvana_docker
        memory: "128 GB"
        cpu: "4"
        preemptible: 1
        maxRetries: 1
        disks: "local-disk 2000 HDD"
    }

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
        String variants_docker
    }

    File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"

    String output_vt_json = "vat_vt_bq_load" + output_file_suffix
    String output_vt_gcp_path = output_path + 'vt/'

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
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

    runtime {
        docker: variants_docker
        memory: "16 GB"
        preemptible: 2
        cpu: "1"
        disks: "local-disk 500 HDD"
    }

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
        String variants_docker
    }

    File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"

    String output_genes_json = "vat_genes_bq_load" + output_file_suffix
    String output_genes_gcp_path = output_path + 'genes/'

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
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

    runtime {
        docker: variants_docker
        memory: "7 GB"
        preemptible: 3
        cpu: "1"
        disks: "local-disk 500 HDD"
    }

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
        String base_vat_table_name
        File? nirvana_schema
        File? vt_schema
        File? genes_schema
        String project_id
        String dataset_name
        String output_path
        Array[Boolean] prep_vt_json_done
        Array[Boolean] prep_genes_json_done
        String cloud_sdk_docker
    }

    # This is the name of the vat table. Due to sharding (VS-1191) there may be some duplicated entries.
    # So we create it here, and then deduplicate it in a later step
    String vat_table_name = base_vat_table_name + "_w_dups"

    String variant_transcript_table = base_vat_table_name + "_variants"
    String genes_table = base_vat_table_name + "_genes"

    String vt_path = output_path + 'vt/*'
    String genes_path = output_path + 'genes/*'

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        echo "project_id = ~{project_id}" > ~/.bigqueryrc

        DATE=86400 ## 24 hours in seconds

        set +o errexit
        bq --apilog=false show --project_id=~{project_id} ~{dataset_name}.~{variant_transcript_table} > /dev/null
        BQ_SHOW_RC=$?
        set -o errexit

        if [ $BQ_SHOW_RC -ne 0 ]; then
            echo "Creating a pre-vat table ~{dataset_name}.~{variant_transcript_table}"
            bq --apilog=false mk --expiration=$DATE --project_id=~{project_id}  ~{dataset_name}.~{variant_transcript_table} ~{vt_schema}
        fi

        echo "Loading data into a pre-vat table ~{dataset_name}.~{variant_transcript_table}"
        echo ~{vt_path}
        echo ~{genes_path}
        bq --apilog=false load --project_id=~{project_id} --source_format=NEWLINE_DELIMITED_JSON ~{dataset_name}.~{variant_transcript_table} ~{vt_path}

        set +o errexit
        bq --apilog=false show --project_id=~{project_id} ~{dataset_name}.~{genes_table} > /dev/null
        BQ_SHOW_RC=$?
        set -o errexit

        if [ $BQ_SHOW_RC -ne 0 ]; then
            echo "Creating a pre-vat table ~{dataset_name}.~{genes_table}"
            bq --apilog=false mk --expiration=$DATE --project_id=~{project_id}  ~{dataset_name}.~{genes_table} ~{genes_schema}
        fi

        echo "Loading data into a pre-vat table ~{dataset_name}.~{genes_table}"
        bq --apilog=false load  --project_id=~{project_id} --source_format=NEWLINE_DELIMITED_JSON  ~{dataset_name}.~{genes_table} ~{genes_path}

        set +e
        bq --apilog=false show --project_id=~{project_id} ~{dataset_name}.~{vat_table_name} > /dev/null
        BQ_SHOW_RC=$?
        set -e

        if [ $BQ_SHOW_RC -ne 0 ]; then
            echo "Creating the vat table ~{dataset_name}.~{vat_table_name}"
        else
            echo "Dropping and recreating the vat table ~{dataset_name}.~{vat_table_name}"
            bq --apilog=false rm -t -f --project_id=~{project_id} ~{dataset_name}.~{vat_table_name}
        fi

        CLUSTERING_STRING="--clustering_fields=contig"
        bq --apilog=false mk ${CLUSTERING_STRING} --expiration=$DATE --project_id=~{project_id} ~{dataset_name}.~{vat_table_name} ~{nirvana_schema}
        echo "Loading data into it"


        # Now we run a giant query in BQ to get this all in the right table and join the genes properly
        # Note the genes table join includes the group by to avoid the duplicates that get created from genes that span shards
        # Commented out columns in the query are to be added in the next release
        # We want the vat creation query to overwrite the destination table because if new data has been put into the pre-vat tables
        # and this workflow has been run an additional time, we dont want duplicates being appended from the original run

        # bq query --max_rows check: ok selecting into a table
        bq --apilog=false query --nouse_legacy_sql --destination_table=~{dataset_name}.~{vat_table_name} --replace --project_id=~{project_id} \
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
            v.gnomad_nfe_ac,
            v.gnomad_nfe_an,
            v.gnomad_nfe_af,
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

    runtime {
        docker: cloud_sdk_docker
        memory: "3 GB"
        preemptible: 3
        cpu: "1"
        disks: "local-disk 100 HDD"
    }

    output {
        String vat_table = vat_table_name
        Boolean done = true
    }
}

task DeduplicateVatInBigQuery {
    meta {
        # since the WDL will not see the updated data (it's getting put in a gcp bucket)
        volatile: true
    }

    input {
        String input_vat_table_name
        String output_vat_table_name
        File? nirvana_schema

        String project_id
        String dataset_name
        String cloud_sdk_docker
    }


    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        echo "project_id = ~{project_id}" > ~/.bigqueryrc

        DATE=86400 ## 24 hours in seconds

        set +e
        bq --apilog=false show --project_id=~{project_id} ~{dataset_name}.~{output_vat_table_name} > /dev/null
        BQ_SHOW_RC=$?
        set -e

        if [ $BQ_SHOW_RC -ne 0 ]; then
            echo "Creating the final vat table ~{dataset_name}.~{output_vat_table_name}"
        else
            bq --apilog=false rm -t -f --project_id=~{project_id} ~{dataset_name}.~{output_vat_table_name}
        fi
        bq --apilog=false mk --project_id=~{project_id} ~{dataset_name}.~{output_vat_table_name} ~{nirvana_schema}
        echo "Loading data into it"

        # Now we query the original VAT table and recreate it, but remove any rows that appear twice.

        # bq query --max_rows check: ok selecting into a table
        bq --apilog=false query --nouse_legacy_sql --destination_table=~{dataset_name}.~{output_vat_table_name} --replace --project_id=~{project_id} \
        ' SELECT * EXCEPT(row_number) FROM (
            SELECT
                *,
                row_number()
                    over (partition by vid, transcript)
                    row_number
            FROM
                `~{dataset_name}.~{input_vat_table_name}`
            )
            where row_number = 1'
    >>>

    runtime {
        docker: cloud_sdk_docker
        memory: "3 GB"
        preemptible: 3
        cpu: "1"
        disks: "local-disk 100 HDD"
    }

    output {
        String vat_table = output_vat_table_name
        Boolean done = true
    }
}
