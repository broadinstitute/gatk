version 1.0

workflow GvsTestReference {
  call TestReference {
    input:
#      input_vcf = StripCustomAnnotationsFromSitesOnlyVCF.output_vcf,
#      output_annotated_file_name = "${vcf_filename}_annotated",
#      custom_annotations_file = StripCustomAnnotationsFromSitesOnlyVCF.output_custom_annotations_file,
      variants_nirvana_docker = "us.gcr.io/broad-dsde-methods/variantstore:nirvana_2022_10_19",
      use_reference_disk = true,
  }

}


task TestReference {
  input {
#    File input_vcf
#    String output_annotated_file_name
#    File custom_annotations_file

    # Mentioning this path in the inputs section of the task combined with checking the 'Use reference disks' option
    # in Terra UI tells Cromwell to arrange for the Nirvana reference disk to be attached to this VM.
#    File summon_reference_disk =
#                               "gs://gcp-public-data--broad-references/hg38/v0/Nirvana/3.18.1_2024-03-06/SupplementaryAnnotation/GRCh38/ClinGen_Dosage_Sensitivity_Map_20231105.nsi"
    File summon_reference_disk =
                               "gs://broad-public-datasets/gvs/vat-annotations/Nirvana/3.18.1/SupplementaryAnnotation/GRCh38/ClinGen_Dosage_Sensitivity_Map_20220512.nsi"

    String variants_nirvana_docker

    File omim_annotations = "gs://gcp-public-data--broad-references/hg38/v0/Nirvana/3.18.1_2024-03-06/SupplementaryAnnotation/GRCh38/OMIM_20220516.nga"
    Boolean use_reference_disk
  }

  File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"

#  String annotation_json_name = output_annotated_file_name + ".json.gz"
#  String gene_annotation_json_name = output_annotated_file_name + ".genes.json.gz"
#  String positions_annotation_json_name = output_annotated_file_name + ".positions.json.gz"
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
      echo "This is the /mnt:"
      ls -altrh /mnt
      echo "Done"
      lsblk
      echo "That was lsblk"
      CANDIDATE_MOUNT_POINT=$(lsblk | sed -E -n 's!.*(/mnt/[a-f0-9]+).*!\1!p')
      if [[ -z ${CANDIDATE_MOUNT_POINT} ]]; then
        >&2 echo "Could not find a mounted volume that looks like a reference disk, exiting."
        exit 1
      fi
      echo "C: ${CANDIDATE_MOUNT_POINT}"
      ls -altrh ${CANDIDATE_MOUNT_POINT}
      echo "That was the contents of that mount point"

      # Find one particular reference under the mount path. Note this is not the same reference as was specified in the
      # `inputs` section, so this would only be present if the volume we're looking at is in fact a reference disk.
      REFERENCE_FILE="Homo_sapiens.GRCh38.Nirvana.dat"
      REFERENCE_PATH=$(find ${CANDIDATE_MOUNT_POINT} -name "${REFERENCE_FILE}")
      if [[ -z ${REFERENCE_PATH} ]]; then
        >&2 echo "Could not find reference file '${REFERENCE_FILE}' under candidate reference disk mount point '${CANDIDATE_MOUNT_POINT}', exiting."
        exit 1
      fi
      echo "D: ${REFERENCE_PATH}"

      # Take the parent of the parent directory of this file as root of the locally mounted references:
      DATA_SOURCES_FOLDER="$(dirname $(dirname ${REFERENCE_PATH}))"
    else
      DATA_SOURCES_FOLDER=/cromwell_root/nirvana_references
      mkdir ${DATA_SOURCES_FOLDER}

      # Download the references
      dotnet /Nirvana/Downloader.dll --ga GRCh38 --out ${DATA_SOURCES_FOLDER}

      # As of 2024-01-24 OMIM is no longer included among the bundle of annotation resources pulled down by the
      # Nirvana downloader. As this annotation set is currently central for our VAT logic, special-case link in
      # the OMIM .nsa bundle we downloaded back when we made the Delta reference disk:
      ln ~{omim_annotations} ${DATA_SOURCES_FOLDER}/SupplementaryAnnotation/GRCh38/
      # TODO - should link the other four.
    fi

    echo "Looking in ${DATA_SOURCES_FOLDER}"
    ls -l ${DATA_SOURCES_FOLDER}
    echo "Find in it"
    find ${DATA_SOURCES_FOLDER} -name '*' -exec ls -l {} \;
    echo "That's it"

    # =======================================
    echo "Creating custom annotations"
    mkdir customannotations_dir
    CUSTOM_ANNOTATIONS_FOLDER="$PWD/customannotations_dir"

  >>>

  runtime {
    docker: variants_nirvana_docker
    memory: "64 GB"
    cpu: "4"
    preemptible: 3
    maxRetries: 2
    disks: "local-disk 2000 HDD"
  }

  output {
#    File genes_annotation_json = "~{gene_annotation_json_name}"
#    File positions_annotation_json = "~{positions_annotation_json_name}"
    File monitoring_log = "monitoring.log"
  }
}