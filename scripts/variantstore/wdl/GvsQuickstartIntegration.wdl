version 1.0

import "GvsUnified.wdl" as Unified
import "GvsUtils.wdl" as Utils

workflow GvsQuickstartIntegration {

    input {
        String branch_name
        String expected_output_prefix = "gs://gvs-internal-quickstart/integration/2022-06-28/"

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
                                 "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00405.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz",
                                 "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00408.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz",
                                 "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00418.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz",
                                 "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00420.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz",
                                 "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00423.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz",
                                 "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00427.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz",
                                 "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00429.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz",
                                 "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00444.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz",
                                 "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00447.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz",
                                 "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00450.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz"
                                 ]

        Array[File] input_vcf_indexes = [
                                        "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00405.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz.tbi",
                                        "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00408.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz.tbi",
                                        "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00418.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz.tbi",
                                        "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00420.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz.tbi",
                                        "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00423.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz.tbi",
                                        "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00427.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz.tbi",
                                        "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00429.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz.tbi",
                                        "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00444.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz.tbi",
                                        "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00447.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz.tbi",
                                        "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00450.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz.tbi"
                                        ]

        Int? extract_scatter_count

        Int load_data_batch_size = 1
    }

    call Utils.BuildGATKJarAndCreateDataset {
        input:
            branch_name = branch_name,
            dataset_prefix = "quickit"
    }

    call Unified.GvsUnified {
        input:
            dataset_name = BuildGATKJarAndCreateDataset.dataset_name,
            project_id = "gvs-internal",
            external_sample_names = external_sample_names,
            gatk_override = BuildGATKJarAndCreateDataset.jar,
            input_vcfs = input_vcfs,
            input_vcf_indexes = input_vcf_indexes,
            filter_set_name = "quickit",
            extract_table_prefix = "quickit",
            extract_scatter_count = extract_scatter_count,
            # Force filtering off as it is not deterministic and the initial version of this integration test does not
            # allow for inexact matching of actual and expected results.
            extract_do_not_filter_override = true,
            load_data_batch_size = load_data_batch_size,
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


task AssertIdenticalOutputs {
    input {
        String expected_output_prefix
        Array[File] actual_vcfs
    }

    command <<<
        set -o errexit
        set -o nounset
        set -o pipefail
        set -o xtrace

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
        gzip -d *.gz
        cd ..

        # Also unzip actual result data
        gzip -d ~{sep= ' ' actual_vcfs}

        # Headers first, these can yield useful diagnostics when there are mismatches.
        for file in ~{sep=' ' actual_vcfs}; do
          unzipped=${file%.gz}
          expected="expected/$(basename $unzipped)"
          set +o errexit
          cmp <(grep '^#' $unzipped) <(grep '^#' $expected)
          rc=$?
          set -o errexit
          if [[ $rc -ne 0 ]]; then
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
          set +o errexit
          cmp $unzipped $expected
          rc=$?
          set -o errexit
          if [[ $rc -ne 0 ]]; then
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
