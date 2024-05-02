version 1.0

import "GvsUtils.wdl" as Utils
import "GvsJointVariantCalling.wdl" as JointVariantCalling

workflow GvsQuickstartVcfIntegration {
    input {
        String git_branch_or_tag
        String? git_hash
        String expected_output_prefix
        Boolean use_VQSR_lite = true
        Boolean extract_do_not_filter_override = true
        Boolean use_compressed_references = false
        Boolean load_vcf_headers = false
        String drop_state = "FORTY"
        Boolean bgzip_output_vcfs = false
        String dataset_suffix
        Boolean is_wgs = true
        File? interval_list
        Boolean use_default_dockers = false
        Boolean check_expected_cost_and_table_size_outputs = true
        String? basic_docker
        String? cloud_sdk_docker
        String? cloud_sdk_slim_docker
        String? variants_docker
        String? gatk_docker
        File? gatk_override
        String? sample_id_column_name ## Note that a column WILL exist that is the <entity>_id from the table name. However, some users will want to specify an alternate column for the sample_name during ingest
        String? vcf_files_column_name
        String? vcf_index_files_column_name
        String? sample_set_name ## NOTE: currently we only allow the loading of one sample set at a time

        String? workspace_bucket
        String? workspace_id
        String? submission_id
    }
    String project_id = "gvs-internal"

    # WDL 1.0 trick to set a variable ('none') to be undefined.
    if (false) {
      File? none = ""
    }

    if (!defined(workspace_bucket) || !defined(workspace_id) || !defined(submission_id) ||
        !defined(git_hash) || !defined(cloud_sdk_docker) || !defined(cloud_sdk_slim_docker) ||
        !defined(variants_docker) || !defined(basic_docker) || !defined(gatk_docker)) {
        call Utils.GetToolVersions {
            input:
                git_branch_or_tag = git_branch_or_tag,
        }
    }

    String effective_basic_docker = select_first([basic_docker, GetToolVersions.basic_docker])
    String effective_cloud_sdk_docker = select_first([cloud_sdk_docker, GetToolVersions.cloud_sdk_docker])
    String effective_cloud_sdk_slim_docker = select_first([cloud_sdk_slim_docker, GetToolVersions.cloud_sdk_slim_docker])
    String effective_variants_docker = select_first([variants_docker, GetToolVersions.variants_docker])
    String effective_gatk_docker = select_first([gatk_docker, GetToolVersions.gatk_docker])
    String effective_git_hash = select_first([git_hash, GetToolVersions.git_hash])

    String effective_workspace_bucket = select_first([workspace_bucket, GetToolVersions.workspace_bucket])
    String effective_workspace_id = select_first([workspace_id, GetToolVersions.workspace_id])
    String effective_submission_id = select_first([submission_id, GetToolVersions.submission_id])

    if (!use_default_dockers && !defined(gatk_override)) {
      call Utils.BuildGATKJar {
        input:
          git_branch_or_tag = git_branch_or_tag,
          cloud_sdk_slim_docker = effective_cloud_sdk_slim_docker,
      }
    }

    call Utils.CreateDatasetForTest {
        input:
            git_branch_or_tag = git_branch_or_tag,
            dataset_prefix = "quickit",
            dataset_suffix = dataset_suffix,
            cloud_sdk_docker = effective_cloud_sdk_docker,
    }

    call JointVariantCalling.GvsJointVariantCalling as JointVariantCalling {
        input:
            call_set_identifier = git_branch_or_tag,
            dataset_name = CreateDatasetForTest.dataset_name,
            project_id = project_id,
            gatk_override = if (use_default_dockers) then none else select_first([gatk_override, BuildGATKJar.jar]),
            use_classic_VQSR = !use_VQSR_lite,
            use_compressed_references = use_compressed_references,
            load_vcf_headers = load_vcf_headers,
            extract_output_file_base_name = "quickit",
            filter_set_name = "quickit",
            extract_table_prefix = "quickit",
            # optionally turn off filtering (VQSR Classic is not deterministic)
            # (and the initial version of this integration test does not allow for inexact matching of actual and expected results.)
            extract_do_not_filter_override = extract_do_not_filter_override,
            drop_state = drop_state,
            bgzip_output_vcfs = bgzip_output_vcfs,
            is_wgs = is_wgs,
            interval_list = interval_list,
            sample_id_column_name = sample_id_column_name,
            vcf_files_column_name = vcf_files_column_name,
            vcf_index_files_column_name = vcf_index_files_column_name,
            sample_set_name = sample_set_name,
            basic_docker = effective_basic_docker,
            cloud_sdk_docker = effective_cloud_sdk_docker,
            variants_docker = effective_variants_docker,
            gatk_docker = effective_gatk_docker,
            workspace_bucket = effective_workspace_bucket,
            workspace_id = effective_workspace_id,
            submission_id = effective_submission_id,
            git_branch_or_tag = git_branch_or_tag,
            git_hash = effective_git_hash,
            tighter_gcp_quotas = false,
    }

    # Only assert identical outputs if we did not filter (filtering is not deterministic) OR if we are using VQSR Lite (which is deterministic)
    if (extract_do_not_filter_override || use_VQSR_lite) {
        String expected_prefix = expected_output_prefix + dataset_suffix + "/"
        call AssertIdenticalOutputs {
            input:
                expected_output_prefix = expected_prefix,
                actual_vcfs = JointVariantCalling.output_vcfs,
                cloud_sdk_docker = effective_cloud_sdk_docker,
        }

        if (check_expected_cost_and_table_size_outputs) {
            call AssertCostIsTrackedAndExpected {
                input:
                    go = JointVariantCalling.done,
                    dataset_name = CreateDatasetForTest.dataset_name,
                    project_id = project_id,
                    expected_output_csv = expected_prefix + "cost_observability.csv",
                    cloud_sdk_docker = effective_cloud_sdk_docker,
            }

            call AssertTableSizesAreExpected {
                input:
                    go = JointVariantCalling.done,
                    dataset_name = CreateDatasetForTest.dataset_name,
                    project_id = project_id,
                    expected_output_csv = expected_prefix + "table_sizes.csv",
                    cloud_sdk_docker = effective_cloud_sdk_docker,
            }
        }
    }

    output {
        Array[File] output_vcfs = JointVariantCalling.output_vcfs
        Array[File] output_vcf_indexes = JointVariantCalling.output_vcf_indexes
        Float total_vcfs_size_mb = JointVariantCalling.total_vcfs_size_mb
        File manifest = JointVariantCalling.manifest
        String dataset_name = CreateDatasetForTest.dataset_name
        String filter_set_name = "quickit"
        String recorded_git_hash = effective_git_hash
        Boolean done = true
        Boolean used_tighter_gcp_quotas = JointVariantCalling.used_tighter_gcp_quotas
    }
}


task AssertIdenticalOutputs {
    input {
        String expected_output_prefix
        Array[File] actual_vcfs
        String cloud_sdk_docker
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
        docker: cloud_sdk_docker
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
        String cloud_sdk_docker
    }

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        mkdir output

        echo "project_id = ~{project_id}" > ~/.bigqueryrc
        bq --apilog=false query --project_id=~{project_id} --format=csv --use_legacy_sql=false \
            'SELECT call, step, event_key, sum(event_bytes)
              FROM `~{dataset_name}.cost_observability`
              GROUP BY call, step, event_key
              ORDER BY call, step, event_key' > output/cost_observability.csv

        # Put the exit code in a file because we are using a subshell (while) later and changes to the variable *in* the subshell are lost
        echo "0" > ret_val.txt

        paste output/cost_observability.csv ~{expected_output_csv} | while  IFS=$'\t' read observed expected;
        do
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
              TOLERANCE=1.0   # 100% tolerance
            elif [[ $OBS_KEY == "ExtractTask.GvsExtractCallset.Storage API Scanned" ]]; then
              TOLERANCE=0.1   # 10% tolerance
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
        docker: cloud_sdk_docker
        disks: "local-disk 10 HDD"
    }

    output {
        File cost_observability_output_csv = "output/cost_observability.csv"
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
        String cloud_sdk_docker
    }

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        mkdir output

        echo "project_id = ~{project_id}" > ~/.bigqueryrc
        bq --apilog=false query --project_id=~{project_id} --format=csv --use_legacy_sql=false \
            "SELECT 'vet_total' AS total_name, sum(total_billable_bytes) AS total_bytes FROM \
            \`~{dataset_name}.INFORMATION_SCHEMA.PARTITIONS\` WHERE table_name LIKE 'vet_%' \
            UNION ALL \
            SELECT 'ref_ranges_total' AS total_name, sum(total_billable_bytes) AS total_bytes \
            FROM \`~{dataset_name}.INFORMATION_SCHEMA.PARTITIONS\` \
            WHERE table_name LIKE 'ref_ranges_%' ORDER BY total_name" > output/table_sizes.csv

        set +o errexit
        diff -w output/table_sizes.csv ~{expected_output_csv} > differences.txt
        set -o errexit

        if [[ -s differences.txt ]]; then
            echo "Differences found:"
            cat differences.txt
            exit 1
        fi
    >>>

    runtime {
        docker: cloud_sdk_docker
        disks: "local-disk 10 HDD"
    }

    output {
        File table_sizes_output_csv = "output/table_sizes.csv"
        File differences = "differences.txt"
    }
}
