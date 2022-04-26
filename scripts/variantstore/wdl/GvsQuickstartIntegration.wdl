version 1.0

import "GvsUnified.wdl" as GvsUnified

workflow GvsQuickstartIntegration {

    input {
        String branch_name
        String expected_output_prefix = "gs://broad-dsp-spec-ops/quickstart_integration/2022-04-25/"
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
            extract_scatter_count = 100,
            # Force filtering off as it is not deterministic and the initial version of this integration test does not
            # allow for inexact matching of actual and expected results.
            extract_do_not_filter_override = true
    }

    call AssertIdenticalOutputs {
        input:
            expected_output_prefix = expected_output_prefix,
            actual_vcfs = GvsUnified.output_vcfs
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
        # Much of this could/should be put into a Docker image which would be useful not only for integration test runs
        # but also for building nightly GATK jars. This wouldn't be that valuable for this current increment of work as
        # using a Docker image would save maybe 10 minutes from what is currently a ~4 hour workflow, but a Docker image
        # could become more compelling if a scaled down version of this test that could run more frequently was created,
        # in addition to the nightly build use case mentioned above.
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

        # Rename the GATK jar to embed the branch and hash of the most recent commit on the branch.
        mv build/libs/gatk-package-unspecified-SNAPSHOT-local.jar "build/libs/gatk-${branch}-${hash}-SNAPSHOT-local.jar"

        # Build a dataset name based on the branch name and the git hash of the most recent commit on this branch.
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
    }

    command <<<
        failures=()

        # Where the current set of expected results lives in the cloud
        expected_prefix="~{expected_output_prefix}"
        # Remove a trailing slash if there is one
        expected_prefix=${expected_prefix%/}

        # Download all the expected data
        mkdir expected
        cd expected
        # Make a FOFN for more efficient downloading
        for file in ~{sep= ' ' actual_vcfs}; do
            echo "$expected_prefix/$(basename $file)" >> expected_fofn.txt
        done
        # Download and unzip all the expected data
        cat expected_fofn.txt | gsutil -m cp -I .
        gzip -d *
        cd ..

        # Also unzip actual result data
        gzip -d ~{sep= ' ' actual_vcfs}

        # Headers first, these can yield useful diagnostics when there are mismatches.
        for file in ~{sep=' ' actual_vcfs}; do
          unzipped=${file%.gz}
          expected="expected/$(basename $unzipped)"
          cmp <(grep '^#' $unzipped) <(grep '^#' $expected)
          if [[ $? -ne 0 ]]; then
            # If there is a mismatch add it to a list of failures but keep on looking for mismatches.
            failures+=( $unzipped )
          fi
        done

        if [[ ${#failures[@]} -ne 0 ]]; then
          echo "Error: headers for the following files do not match:"
          for failure in ${failures[@]}; do
            echo $(basename $failure)
            expected="expected/$(basename $failure)"
            diff $failure $expected
          done
          exit 1
        fi

        # If the headers all matched look for any mismatches in overall file content.
        fail=0
        for file in ~{sep=' ' actual_vcfs}; do
          unzipped=${file%.gz}
          expected="expected/$(basename $unzipped)"
          cmp $unzipped $expected
          if [[ $? -ne 0 ]]; then
            echo "Error: file contents of expected and actual do not match: $(basename $unzipped)"
            fail=1
          fi
        done

        if [[ $fail -ne 0 ]]; then
          exit 1
        fi

        echo "All vcfs compared and matched!"
    >>>

    runtime {
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:latest"
        disks: "local-disk 500 HDD"
    }

    output {
        File fofn = "expected/expected_fofn.txt"
        Boolean done = true
    }
}
