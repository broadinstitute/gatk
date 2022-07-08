version 1.0

import "GvsUnified.wdl" as Unified
import "GvsUtils.wdl" as Utils

workflow GvsQuickstartIntegration {

    input {
        String branch_name
        String expected_output_prefix = "gs://gvs-internal-quickstart/integration/2022-07-05/"

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
    }
    String project_id = "gvs-internal"

    call Utils.BuildGATKJarAndCreateDataset {
        input:
            branch_name = branch_name,
            dataset_prefix = "quickit"
    }

    call Unified.GvsUnified {
        input:
            call_set_identifier = branch_name,
            dataset_name = BuildGATKJarAndCreateDataset.dataset_name,
            project_id = project_id,
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
    }

    call AssertIdenticalOutputs {
        input:
            expected_output_prefix = expected_output_prefix,
            actual_vcfs = GvsUnified.output_vcfs
    }

    call AssertCostIsTrackedAndExpected {
        input:
            go = GvsUnified.done,
            dataset_name = BuildGATKJarAndCreateDataset.dataset_name,
            project_id = project_id,
            expected_output_csv = expected_output_prefix + "cost_observability_expected.csv"
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

task AssertCostIsTrackedAndExpected {
    meta {
        # we want to check the database each time this runs
        volatile: true
    }

    input {
        Boolean go = true
        String dataset_name
        String project_id
        File expected_output_csv
    }

    command <<<
        set -o errexit
        set -o nounset
        set -o pipefail
        set -o xtrace

        echo "project_id = ~{project_id}" > ~/.bigqueryrc
        bq query --location=US --project_id=~{project_id} --format=csv --use_legacy_sql=false \
            "SELECT call, step, event_key, sum(event_bytes) \
              FROM \`~{dataset_name}.cost_observability\` \
              GROUP BY call, step, event_key \
              ORDER BY call, step, event_key" > cost_observability_output.csv

        # Put the exit code in a file because we are using a subshell (while) later and changes to the variable *in* the subshell are lost
        echo "0" > ret_val.txt

        paste cost_observability_output.csv ~{expected_output_csv} | while  IFS=$'\t' read observed expected; do
        IFS=, read -ra OBS <<< "$observed"
        IFS=, read -ra EXP <<< "$expected"
        if [[ "${#OBS[@]}" -ne 4  || "${#EXP[@]}" -ne 4 ]]; then
          echo "Unexpected number of rows found in the input files"
          exit 1
        fi

        OBS_KEY=${OBS[0]}.${OBS[1]}.${OBS[2]}
        EXP_KEY=${EXP[0]}.${EXP[1]}.${EXP[2]}
        if [[ "$OBS_KEY" != "$EXP_KEY" ]]; then
          echo "Mismatched keys in results files - were these sorted properly?"
          exit 1
        fi

        if [[ "$OBS_KEY" == "call.step.event_key" ]]; then
          # Skip the header
          continue;
        fi

        OBS_BYTES=${OBS[3]}
        EXP_BYTES=${EXP[3]}

        TOLERANCE=0

        # For these two costs, there is non-determinism in the pipeline - we allow a % difference
        if [[ $OBS_KEY == "ExtractFilterTask.GvsCreateFilterSet.BigQuery Query Scanned" ]]; then
          TOLERANCE=0.015   # 1.5% tolerance  (Note - have seen as high as: 0.0109429)
        elif [[ $OBS_KEY == "ExtractTask.GvsCreateCallset.Storage API Scanned" ]]; then
          TOLERANCE=0.01   # 1% tolerance  (Note - have seen as high as: 0.00608656)
        elif [[ $OBS_KEY == "ExtractFilterTask.GvsCreateFilterSet.Storage API Scanned" ]]; then
          TOLERANCE=0.05  # 5% tolerance (Note - have seen as high as: 0.0281223)
        fi

        if [[ $OBS_BYTES -ne $EXP_BYTES ]]; then
          echo "The bytes observed ($OBS_BYTES) for '$OBS_KEY' differ from those expected ($EXP_BYTES)"

          if [[ $OBS_BYTES -ge $EXP_BYTES ]]; then
            DIFF_FOUND=$(echo $OBS_BYTES $EXP_BYTES | awk '{print ($1-$2)/$1}')
          else
            DIFF_FOUND=$(echo $EXP_BYTES $OBS_BYTES | awk '{print ($1-$2)/$1}')
          fi

          if ! awk "BEGIN{ exit ($DIFF_FOUND > $TOLERANCE) }"
          then
            echo "FAIL!!! The relative difference between these is $DIFF_FOUND, which is greater than the allowed tolerance ($TOLERANCE)"
            echo "1" > ret_val.txt
          else
            echo "However, the relative difference between these is $DIFF_FOUND, which is below the allowed tolerance ($TOLERANCE)"
          fi
        fi
        done

        RET_VAL=`cat ret_val.txt`
        exit $RET_VAL

    >>>

    runtime {
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:latest"
        disks: "local-disk 10 HDD"
    }

    output {
      File cost_observability_output_csv = "cost_observability_output.csv"
    }
}
