version 1.0

import "GvsUnified.wdl" as GvsUnified

workflow GvsQuickstartIntegration {

    input {
        String branch_name
        String expected_output_prefix = "gs://broad-dsp-spec-ops/quickstart_integration/2022-04-22/"
    }

    Array[String] external_sample_names = [
        "ERS4367795",
        "ERS4367796",
        "ERS4367797",
        "ERS4367798",
        "ERS4367799",
        "ERS4367800",
        "ERS4367801",
        "ERS4367803",
        "ERS4367804",
        "ERS4367805"
    ]

    Array[File] input_vcfs = [
        "gs://fc-2b4456d7-974b-4b67-90f8-63c2fd2c03d4/gvcfs/HG00405.haplotypeCalls.er.raw.vcf.gz.vcf.gz",
        "gs://fc-2b4456d7-974b-4b67-90f8-63c2fd2c03d4/gvcfs/HG00408.haplotypeCalls.er.raw.vcf.gz.vcf.gz",
        "gs://fc-2b4456d7-974b-4b67-90f8-63c2fd2c03d4/gvcfs/HG00418.haplotypeCalls.er.raw.vcf.gz.vcf.gz",
        "gs://fc-2b4456d7-974b-4b67-90f8-63c2fd2c03d4/gvcfs/HG00420.haplotypeCalls.er.raw.vcf.gz.vcf.gz",
        "gs://fc-2b4456d7-974b-4b67-90f8-63c2fd2c03d4/gvcfs/HG00423.haplotypeCalls.er.raw.vcf.gz.vcf.gz",
        "gs://fc-2b4456d7-974b-4b67-90f8-63c2fd2c03d4/gvcfs/HG00427.haplotypeCalls.er.raw.vcf.gz.vcf.gz",
        "gs://fc-2b4456d7-974b-4b67-90f8-63c2fd2c03d4/gvcfs/HG00429.haplotypeCalls.er.raw.vcf.gz.vcf.gz",
        "gs://fc-2b4456d7-974b-4b67-90f8-63c2fd2c03d4/gvcfs/HG00444.haplotypeCalls.er.raw.vcf.gz.vcf.gz",
        "gs://fc-2b4456d7-974b-4b67-90f8-63c2fd2c03d4/gvcfs/HG00447.haplotypeCalls.er.raw.vcf.gz.vcf.gz",
        "gs://fc-2b4456d7-974b-4b67-90f8-63c2fd2c03d4/gvcfs/HG00450.haplotypeCalls.er.raw.vcf.gz.vcf.gz"
    ]

    Array[File] input_vcf_indexes = [
        "gs://fc-2b4456d7-974b-4b67-90f8-63c2fd2c03d4/gvcfs/HG00405.haplotypeCalls.er.raw.vcf.gz.vcf.gz.tbi",
        "gs://fc-2b4456d7-974b-4b67-90f8-63c2fd2c03d4/gvcfs/HG00408.haplotypeCalls.er.raw.vcf.gz.vcf.gz.tbi",
        "gs://fc-2b4456d7-974b-4b67-90f8-63c2fd2c03d4/gvcfs/HG00418.haplotypeCalls.er.raw.vcf.gz.vcf.gz.tbi",
        "gs://fc-2b4456d7-974b-4b67-90f8-63c2fd2c03d4/gvcfs/HG00420.haplotypeCalls.er.raw.vcf.gz.vcf.gz.tbi",
        "gs://fc-2b4456d7-974b-4b67-90f8-63c2fd2c03d4/gvcfs/HG00423.haplotypeCalls.er.raw.vcf.gz.vcf.gz.tbi",
        "gs://fc-2b4456d7-974b-4b67-90f8-63c2fd2c03d4/gvcfs/HG00427.haplotypeCalls.er.raw.vcf.gz.vcf.gz.tbi",
        "gs://fc-2b4456d7-974b-4b67-90f8-63c2fd2c03d4/gvcfs/HG00429.haplotypeCalls.er.raw.vcf.gz.vcf.gz.tbi",
        "gs://fc-2b4456d7-974b-4b67-90f8-63c2fd2c03d4/gvcfs/HG00444.haplotypeCalls.er.raw.vcf.gz.vcf.gz.tbi",
        "gs://fc-2b4456d7-974b-4b67-90f8-63c2fd2c03d4/gvcfs/HG00447.haplotypeCalls.er.raw.vcf.gz.vcf.gz.tbi",
        "gs://fc-2b4456d7-974b-4b67-90f8-63c2fd2c03d4/gvcfs/HG00450.haplotypeCalls.er.raw.vcf.gz.vcf.gz.tbi",
    ]

    call Prepare {
        input:
            branch_name = branch_name
    }

    call GvsUnified.GvsUnified {
        input:
            dataset_name = Prepare.dataset_name,
            project_id = "spec-ops-aou",
            external_sample_names = external_sample_names,
            gatk_override = Prepare.jar,
            input_vcfs = input_vcfs,
            input_vcf_indexes = input_vcf_indexes,
            filter_set_name = "quickit",
            create_filter_set_scatter_count = 20,
            extract_table_prefix = "quickit",
            extract_scatter_count = 100
    }

    call AssertIdenticalOutputs {
        input:
            expected_output_prefix = expected_output_prefix,
            actual_vcfs = GvsUnified.output_vcfs,
            actual_vcf_indexes = GvsUnified.output_vcf_indexes
    }

    output {
        Array[File] output_vcfs = GvsUnified.output_vcfs
        Array[File] output_vcf_indexes = GvsUnified.output_vcf_indexes
        Float total_vcfs_size_mb = GvsUnified.total_vcfs_size_mb
        File manifest = GvsUnified.manifest
        Boolean done = true
    }
}

task Prepare {
    input {
        String branch_name
    }

    command <<<
        set -o errexit -o nounset -o pipefail

        # git and git-lfs
        apt-get -qq update
        apt-get -qq install git git-lfs

        # Java
        apt-get -qq install wget apt-transport-https gnupg
        wget -O - https://packages.adoptium.net/artifactory/api/gpg/key/public | apt-key add -
        echo "deb https://packages.adoptium.net/artifactory/deb $(awk -F= '/^VERSION_CODENAME/{print$2}' /etc/os-release) main" | tee /etc/apt/sources.list.d/adoptium.list
        apt-get -qq update
        apt -qq install -y temurin-11-jdk

        # GATK
        git clone https://github.com/broadinstitute/gatk.git --depth 1 --branch ~{branch_name} --single-branch
        cd gatk
        ./gradlew shadowJar

        branch=$(git symbolic-ref HEAD 2>/dev/null)
        branch=${branch#refs/heads/}

        hash=$(git rev-parse --short HEAD)

        mv build/libs/gatk-package-unspecified-SNAPSHOT-local.jar "build/libs/gatk-${branch}-${hash}-SNAPSHOT-local.jar"

        # Dataset names must be alphanumeric and underscores only. Convert any dashes to underscores, then delete
        # any remaining characters that are not alphanumeric or underscores.
        dataset="$(echo quickit_${branch}_${hash} | tr '-' '_' | tr -c -d '[:alnum:]_')"

        bq mk --project_id="spec-ops-aou" "$dataset"

        echo -n "$dataset" > dataset.txt
    >>>

    output {
        File jar = glob("gatk/build/libs/*-SNAPSHOT-local.jar")[0]
        String dataset_name = read_string("gatk/dataset.txt")
    }

    runtime {
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:latest"
        disks: "local-disk 500 HDD"
    }
}

task AssertIdenticalOutputs {
    input {
        String expected_output_prefix
        Array[File] actual_vcfs
        Array[File] actual_vcf_indexes
    }
    parameter_meta {
        actual_vcfs:
        {
            localization_optional: true
        }
        actual_vcf_indexes:
        {
            localization_optional: true
        }
    }

    command <<<
        fails=()

        for file in ~{sep=' ' actual_vcfs} ~{sep=' ' actual_vcf_indexes}; do
          cmp <(gsutil hash -c $file | tail -n +2) \
              <(gsutil hash -c "~{expected_output_prefix}/$(basename $file)" | tail -n +2)
          if [[ $? -ne 0 ]]; then
            fails+=( $file )
          fi
        done

        if [[ ${#fails[@]} -ne 0 ]]; then
          echo "Error: the checksums for the following files do not match:"
          for fail in ${fails[@]}; do
            echo $fail
          done
          exit 1
        fi

        echo "All vcf and index files compared and matched!"

    >>>

    runtime {
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:latest"
        disks: "local-disk 500 HDD"
    }
}

