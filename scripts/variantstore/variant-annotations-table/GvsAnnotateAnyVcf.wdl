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

        call AnnotateVCF {
            input:
                cromwell_root = GetToolVersions.cromwell_root,
                input_vcf = CopySitesOnlyVcf.output_file_path,
                output_annotated_file_name = "${sites_only_vcf_basename}_annotated",
                variants_nirvana_docker = effective_variants_nirvana_docker,
                use_reference_disk = use_reference_disk,
        }
    }

    output {
        File? genes_annotation_json = AnnotateVCF.genes_annotation_json
        File? positions_annotation_json = AnnotateVCF.positions_annotation_json
        File? monitoring_log = AnnotateVCF.monitoring_log
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
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        bash ~{monitoring_script} > monitoring.log &

        echo_date () { echo "`date "+%Y/%m/%d %H:%M:%S"` $1"; }

        echo_date "VAT: Convert input to BCF format"
        bcftools convert --threads 4 -O b -o sites_only.bcf ~{sites_only_vcf}

        echo_date "VAT: Calculating number of sites with Ns"
        bcftools view --threads 4 -i 'REF~"N"' -O u sites_only.bcf | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' > track_dropped.tsv

        echo_date "VAT: filter out sites with N's in the reference AND sites with AC=0"
        bcftools view --threads 4 -e 'REF~"N" || AC=0' -O b sites_only.bcf -o filtered_sites_only.bcf
        rm sites_only.bcf

        echo_date "VAT: normalize, left align and split multi allelic sites to new lines, remove duplicate lines"
        bcftools norm --threads 4 -m- --check-ref w -f ~{ref} filtered_sites_only.bcf -O b -o normalized.bcf
        rm filtered_sites_only.bcf

        echo_date "VAT: detecting and removing duplicate rows from sites-only VCF"
        bcftools query normalized.bcf -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' | sort | uniq -d > duplicates.tsv

        echo_date "VAT: done with duplicate detection"
        wc -l duplicates.tsv
        echo_date "VAT: Duplicates may have been found"

        if [ -s duplicates.tsv ]; then
        echo_date "VAT: Removing those rows"
        bcftools view --threads 4 normalized.bcf | grep -v -wFf duplicates.tsv > deduplicated.vcf
        else
        echo_date "VAT: No duplicates found"
        bcftools view --threads 4 normalized.bcf -o deduplicated.vcf
        fi
        rm normalized.bcf

        cat duplicates.tsv >> track_dropped.tsv
        rm duplicates.tsv

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