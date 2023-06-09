version 1.0

import "GvsUnified.wdl" as Unified
import "GvsUtils.wdl" as Utils

workflow GvsQuickstartVcfIntegration {

    input {
        String branch_name
        File interval_list
        String expected_output_prefix
        Boolean use_VQSR_lite = true
        Boolean extract_do_not_filter_override = true

        Array[String] external_sample_names = [
                                              "ERS4367795",
                                              "ERS4367796",
                                              "ERS4367797",
                                              ]

        Array[File] input_vcfs = [
                                 "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00405.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz",
                                 "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00408.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz",
                                 "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00418.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz",
                                 ]

        Array[File] input_vcf_indexes = [
                                        "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00405.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz.tbi",
                                        "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00408.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz.tbi",
                                        "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00418.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz.tbi",
                                        ]

        Int? extract_scatter_count
        String drop_state = "FORTY"
        String dataset_suffix
        File? gatk_override
    }
    String project_id = "gvs-internal"

    if (!defined(gatk_override)) {
      call Utils.BuildGATKJar {
        input:
          branch_name = branch_name,
      }
    }

    call Utils.CreateDataset {
        input:
            branch_name = branch_name,
            dataset_prefix = "quickit",
            dataset_suffix = dataset_suffix,
    }

    call Unified.GvsUnified {
        input:
            call_set_identifier = branch_name,
            dataset_name = CreateDataset.dataset_name,
            project_id = project_id,
            external_sample_names = external_sample_names,
            gatk_override = select_first([gatk_override, BuildGATKJar.jar]),
            input_vcfs = input_vcfs,
            input_vcf_indexes = input_vcf_indexes,
            filter_set_name = "quickit",
            use_VQSR_lite = use_VQSR_lite,
            extract_table_prefix = "quickit",
            extract_scatter_count = extract_scatter_count,
            # optionally turn off filtering (VQSR Classic is not deterministic)
            # (and the initial version of this integration test does not allow for inexact matching of actual and expected results.)
            extract_do_not_filter_override = extract_do_not_filter_override,
            drop_state = drop_state,
            interval_list = interval_list,
    }

    # Only assert identical outputs if we did not filter (filtering is not deterministic) OR if we are using VQSR Lite (which is deterministic)
    if (extract_do_not_filter_override || use_VQSR_lite) {
        String expected_prefix = expected_output_prefix + dataset_suffix + "/"
        call AssertIdenticalOutputs {
            input:
                expected_output_prefix = expected_prefix,
                actual_vcfs = GvsUnified.output_vcfs,
        }

        call AssertCostIsTrackedAndExpected {
            input:
                go = GvsUnified.done,
                dataset_name = CreateDataset.dataset_name,
                project_id = project_id,
                expected_output_csv = expected_prefix + "cost_observability_expected.csv",
        }

        call AssertTableSizesAreExpected {
            input:
                go = GvsUnified.done,
                dataset_name = CreateDataset.dataset_name,
                project_id = project_id,
                expected_output_csv = expected_prefix + "table_sizes_expected.csv",
        }
    }

    output {
        Array[File] output_vcfs = GvsUnified.output_vcfs
        Array[File] output_vcf_indexes = GvsUnified.output_vcf_indexes
        Float total_vcfs_size_mb = GvsUnified.total_vcfs_size_mb
        File manifest = GvsUnified.manifest
        String dataset_name = CreateDataset.dataset_name
        String filter_set_name = "quickit"
        Boolean done = true
    }
}


task AssertIdenticalOutputs {
    input {
        String expected_output_prefix
        Array[File] actual_vcfs
    }
    parameter_meta {
        actual_vcfs: {
            localization_optional: true
        }
    }
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        failures=()

        # Where the current set of expected results lives in the cloud
        expected_prefix="~{expected_output_prefix}"
        # Remove a trailing slash if there is one
        expected_prefix=${expected_prefix%/}

        # Download all the expected data
        mkdir expected
        cd expected
        gcloud storage cp -r "${expected_prefix}"'/*.vcf.gz' .
        gzip -d *.gz
        cd ..

        mkdir actual
        cd actual
        touch actual_manifest.txt
        # Making the manifest is pretty uninteresting and very noisy so turn off xtrace temporarily.
        set +o xtrace
        for actual in ~{sep=' ' actual_vcfs}
        do
            echo $actual >> actual_manifest.txt
        done
        set -o xtrace

        cat actual_manifest.txt | gcloud storage cp -I .
        # Unzip actual result data.
        ls -1 | grep -E '\.vcf\.gz$' | xargs gzip -d
        cd ..

        echo "Header Check"
        # Headers first, these can yield useful diagnostics when there are mismatches.
        for vcf in $(ls -1 actual | grep -E '\.vcf$')
        do
          actual="actual/$vcf"
          expected="expected/$vcf"
          set +o errexit
          cmp <(grep '^#' $actual) <(grep '^#' $expected)
          rc=$?
          set -o errexit
          if [[ $rc -ne 0 ]]; then
            # If there is a mismatch add it to a list of failures but keep on looking for mismatches.
            failures+=( $vcf )
          fi
        done

        echo "Header Failure Check"
        if [[ ${#failures[@]} -ne 0 ]]; then
          echo "Error: headers for the following files do not match:"
          for failure in ${failures[@]}; do
            echo $failure
            expected="expected/$failure"
            actual="actual/$failure"
            diff <(grep '^#' $actual) <(grep '^#' $expected)
          done
          exit 1
        fi

        echo "Overall Check"
        # If the headers all matched look for any mismatches in overall file content.
        fail=0
        for vcf in $(ls -1 actual | grep -E '\.vcf$')
        do
          expected="expected/$vcf"
          actual="actual/$vcf"
          set +o errexit
          cmp $actual $expected
          rc=$?
          set -o errexit
          if [[ $rc -ne 0 ]]; then
            echo "Error: file contents of expected and actual do not match: $vcf"
            fail=1
          fi
        done

        if [[ $fail -ne 0 ]]; then
          exit 1
        fi

        echo "All vcfs compared and matched!"
    >>>

    runtime {
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:426.0.0-alpine"
        disks: "local-disk 500 HDD"
    }

    output {
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

        echo "project_id = ~{project_id}" > ~/.bigqueryrc
        bq --apilog=false query --project_id=~{project_id} --format=csv --use_legacy_sql=false \
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
          TOLERANCE=0.05   # 5% tolerance  (Note - have seen as high as: 0.0371646)
        elif [[ $OBS_KEY == "ExtractFilterTask.GvsCreateFilterSet.Storage API Scanned" ]]; then
          TOLERANCE=0.05  # 5% tolerance (Note - have seen as high as: 0.0281223)
        elif [[ $OBS_KEY == "ExtractTask.GvsExtractCallset.BigQuery Query Scanned" ]]; then
          TOLERANCE=0.6   # 60% tolerance (Note - have seen as high as: 0.5) - but it's 210 vs 420.
        elif [[ $OBS_KEY == "ExtractTask.GvsExtractCallset.Storage API Scanned" ]]; then
          TOLERANCE=0.05   # 5% tolerance  (Note - have seen as high as: 0.012165)
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
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:426.0.0-alpine"
        disks: "local-disk 10 HDD"
    }

    output {
        File cost_observability_output_csv = "cost_observability_output.csv"
    }
}

task AssertTableSizesAreExpected {
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

        echo "project_id = ~{project_id}" > ~/.bigqueryrc
        bq --apilog=false query --project_id=~{project_id} --format=csv --use_legacy_sql=false \
            "SELECT 'vet_total' AS total_name, sum(total_billable_bytes) AS total_bytes FROM \
            \`~{dataset_name}.INFORMATION_SCHEMA.PARTITIONS\` WHERE table_name LIKE 'vet_%' \
            UNION ALL \
            SELECT 'ref_ranges_total' AS total_name, sum(total_billable_bytes) AS total_bytes \
            FROM \`~{dataset_name}.INFORMATION_SCHEMA.PARTITIONS\` \
            WHERE table_name LIKE 'ref_ranges_%' ORDER BY total_name" > table_size_output.csv

        set +o errexit
        diff -w table_size_output.csv ~{expected_output_csv} > differences.txt
        set -o errexit

        if [[ -s differences.txt ]]; then
            echo "Differences found:"
            cat differences.txt
            exit 1
        fi
    >>>

    runtime {
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:426.0.0-alpine"
        disks: "local-disk 10 HDD"
    }

    output {
        File table_size_output_csv = "table_size_output.csv"
        File differences = "differences.txt"
    }
}
