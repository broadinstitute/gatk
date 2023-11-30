version 1.0

import "GvsUtils.wdl" as Utils

workflow GetToolVersionsAndGenerateProvenanceReport {
  input {
    String? git_branch_or_tag
    String? git_hash
    String? cloud_sdk_docker
    String? variants_docker
  }

  if (!defined(git_hash) || !defined(cloud_sdk_docker) || !defined(variants_docker)) {
    call Utils.GetToolVersions {
      input:
        git_branch_or_tag = git_branch_or_tag,
    }
  }

  String defined_cloud_sdk_docker = select_first([cloud_sdk_docker, GetToolVersions.cloud_sdk_docker])
  String defined_variants_docker = select_first([variants_docker, GetToolVersions.variants_docker])
  String defined_git_hash = select_first([git_hash, GetToolVersions.git_hash])

  call GenerateProvenanceReport {
    input:
      git_branch_or_tag = git_branch_or_tag,
      git_hash = defined_git_hash,
      cloud_sdk_docker = defined_cloud_sdk_docker,
      variants_docker = defined_variants_docker
  }

  output {
    String effective_git_hash = defined_git_hash
    String effective_cloud_sdk_docker = defined_cloud_sdk_docker
    String effective_variants_docker = defined_variants_docker
    File provenance_report = GenerateProvenanceReport.report
  }
}

task GenerateProvenanceReport {
  input {
    String? git_branch_or_tag
    String git_hash
    String cloud_sdk_docker
    String variants_docker
  }

  command <<<
    set -e

    echo "Provenance Report\n" > report_file.txt
    echo "Git Hash: ~{git_hash}" >> report_file.txt

  >>>

  # ------------------------------------------------
  # Runtime settings:
  runtime {
    docker: variants_docker
    memory: "1 GB"
    preemptible: 3
    cpu: "1"
    disks: "local-disk 100 HDD"
  }
  output {
    File report = "report_file.txt"
  }
}

