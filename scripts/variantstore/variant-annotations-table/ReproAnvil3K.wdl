version 1.0


workflow repro {
    File input_vcf = "gs://fc-34a42a2b-29aa-43f8-a6dd-2a50f7f3408b/submissions/bef81325-9321-4b4f-a741-b5d151122ceb/GvsCreateVATfromVDS/35a55d99-84e3-4cd4-ad57-1929d99f350c/call-StripCustomAnnotationsFromSitesOnlyVCF/shard-44/0000000044-scattered.sites-only-vcf-c05be276-20b4.unannotated.sites_only.vcf"
    String output_annotated_file_name = "0000000044-scattered.sites-only-vcf-c05be276-20b4_annotated"
    File custom_annotations_file = "gs://fc-34a42a2b-29aa-43f8-a6dd-2a50f7f3408b/submissions/bef81325-9321-4b4f-a741-b5d151122ceb/GvsCreateVATfromVDS/35a55d99-84e3-4cd4-ad57-1929d99f350c/call-StripCustomAnnotationsFromSitesOnlyVCF/shard-44/0000000044-scattered.sites-only-vcf-c05be276-20b4.custom_annotations.tsv"
    String cromwell_root = "/mnt/disks/cromwell_root"
    String variants_nirvana_docker = "us.gcr.io/broad-dsde-methods/variantstore:nirvana_2022_10_19"

    call AnnotateVCF {
        input:
            input_vcf = input_vcf,
            output_annotated_file_name = output_annotated_file_name,
            custom_annotations_file = custom_annotations_file,
            cromwell_root = cromwell_root,
            variants_nirvana_docker = variants_nirvana_docker,
            use_reference_disk = true
    }

    output {
        Boolean done = true
    }
}

task AnnotateVCF {
    input {
        File input_vcf
        String output_annotated_file_name
        File custom_annotations_file
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

        set +o errexit
        cat ~{custom_annotations_file} | grep -v '^#' > content_check_file.txt
        set -o errexit

        if [ ! -s content_check_file.txt ]; then
        echo "Found NO custom annotations in ~{custom_annotations_file} skipping annotation of input VCF"
        echo "Creating empty ennotation jsons for subsequent tasks"
        touch ~{gene_annotation_json_name}
        touch ~{positions_annotation_json_name}
        exit 0
        fi

        if [[ "~{use_reference_disk}" == "true" ]]
        then
        # There's an issue with how the projects/broad-dsde-cromwell-dev/global/images/nirvana-3-18-1-references-2023-01-03
        # disk image was built: while all the reference files do exist on the image they are not at the expected
        # locations. The following code works around this issue and should continue to work even after a corrected
        # version of the Nirvana reference image is deployed into Terra.

        # Find where the reference disk should have been mounted on this VM.  Note this is referred to as a "candidate
        # mount point" because we do not actually confirm this is a reference disk until the following code block.
        if [[ -e /mnt/disks/cromwell_root/gcs_delocalization.sh ]]
        then
        # GCP Batch mount points
        CANDIDATE_MOUNT_POINT=$(lsblk | grep -v cromwell_root | sed -E -n 's!.*(/mnt/disks/[a-f0-9]+).*!\1!p')
        elif [[ -e /cromwell_root/gcs_delocalization.sh ]]
        then
        # PAPI mount points
        CANDIDATE_MOUNT_POINT=$(lsblk | grep -v cromwell_root | sed -E -n 's!.*(/mnt/[a-f0-9]+).*!\1!p')
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

    runtime {
        docker: variants_nirvana_docker
        memory: "128 GB"
        cpu: "4"
        preemptible: 1
        maxRetries: 1
        disks: "local-disk 2000 HDD"
    }

    output {
        File? output_annotated_file_name_out = "~{output_annotated_file_name}"
        File? annotation_json_name_out = "~{annotation_json_name}"
        File genes_annotation_json = "~{gene_annotation_json_name}"
        File positions_annotation_json = "~{positions_annotation_json_name}"
        File monitoring_log = "monitoring.log"
    }
}



