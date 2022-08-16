version 1.0

# This WDL is meant to reblock GVCF files that exist in the DRC buckets.
# The localization and delocalization of DRC files has been combined into one task to save money.
# The reblocked Gvcf will be put in an AoU DRC research bucket
workflow GvsAoUReblockGvcf {

  input {
    String gvcf
    String? gvcf_index
    File ref_dict
    File ref_fasta
    File ref_fasta_index

    String? requester_pays_project
    String? site_id
    String docker_image = "us.gcr.io/broad-gatk/gatk:4.2.2.0"

  }

  call ReblockAndCopy {
    input:
     gvcf = gvcf,
     gvcf_index = select_first([gvcf_index, gvcf + ".tbi"]),
     ref_fasta = ref_fasta,
     ref_fasta_index = ref_fasta_index,
     ref_dict = ref_dict,
     output_gvcf_filename = basename(gvcf, ".g.vcf.gz") + ".reblocked.g.vcf.gz",
     site_id = site_id,
     requester_pays_project = requester_pays_project,
     docker_image = docker_image,
     path = "ss_vcf_research/"
  }

  output {
    String reblocked_gvcf = ReblockAndCopy.output_gvcf
    String reblocked_gvcf_index = ReblockAndCopy.output_gvcf_index
  }
}


task ReblockAndCopy {

  input {
    String gvcf
    String gvcf_index
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    String output_gvcf_filename
    String docker_image

    String? site_id
    String? requester_pays_project
    String path
  }

  Int disk_size = (ceil(60 + size(ref_fasta, "GiB") + size(ref_dict, "GiB")) * 2) + 20

  # this is WDL for .toLowerCase() for the three possible site codes
  String site_id_lower = if defined(site_id) then (
    sub(sub(sub(select_first([site_id]), "BI", "bi"), "BCM", "bcm"), "UW", "uw")
  ) else "null"

  String dir =
      if site_id_lower == "bi" then "gs://prod-genomics-data-broad/" else
      (if site_id_lower == "bcm" then "gs://prod-genomics-data-baylor/" else
      (if site_id_lower == "uw" then "gs://prod-genomics-data-northwest/" else "null" ))

  String destination = dir + path

  command {
    set -euo pipefail
    set -x

    if [ ~{site_id} -a ~{dir} == "null" ]; then
      echo "The additional output directory for the reblocked gVCFs cannot be set to a valid value because the site_id is '~{site_id}' and the only valid values are ['bi', 'bcm', 'uw']."
      exit 1
    fi

    gatk --java-options "-Xms3g -Xmx3g" \
      ReblockGVCF \
      -V ~{gvcf} \
      -do-qual-approx \
      --floor-blocks -GQB 20 -GQB 30 -GQB 40 \
      -O ~{output_gvcf_filename} \
      -R ~{ref_fasta} \
      ~{"--gcs-project-for-requester-pays " + requester_pays_project}

    if [ ~{dir} ]; then
      gsutil -m cp ~{output_gvcf_filename} ~{output_gvcf_filename}.tbi ~{destination}
    fi
  }

  runtime {
    memory: "5 GB"
    bootDiskSizeGb: "15"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 3
    docker: docker_image
  }

  output {
    String output_gvcf = destination + output_gvcf_filename
    String output_gvcf_index = destination + output_gvcf_filename + ".tbi"
  }
}
