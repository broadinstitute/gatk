version 1.0

import "../GvsUtils.wdl" as Utils
import "../../variant-annotations-table/GvsCreateVATfromVDS.wdl" as CreateVATFromVDS

workflow GvsQuickstartVcfIntegration {
    input {
        String git_branch_or_tag
        String? git_hash
        String vds_path
        String ancestry_path
        String output_path
        String split_intervals_scatter_count = 10
        String expected_output_prefix
        String dataset_suffix
        Boolean use_default_dockers = false
        String? basic_docker
        String? cloud_sdk_docker
        String? cloud_sdk_slim_docker
        String? variants_docker
        String? variants_nirvana_docker
        String? gatk_docker
        File? gatk_override

        String? workspace_bucket
        String? workspace_id
        String? submission_id
    }
    String project_id = "gvs-internal"

    # WDL 1.0 trick to set a variable ('none') to be undefined.
    if (false) {
      File? none = ""
    }

    call Utils.GetToolVersions {
        input:
            git_branch_or_tag = git_branch_or_tag,
    }

    String effective_basic_docker = select_first([basic_docker, GetToolVersions.basic_docker])
    String effective_cloud_sdk_docker = select_first([cloud_sdk_docker, GetToolVersions.cloud_sdk_docker])
    String effective_cloud_sdk_slim_docker = select_first([cloud_sdk_slim_docker, GetToolVersions.cloud_sdk_slim_docker])
    String effective_variants_docker = select_first([variants_docker, GetToolVersions.variants_docker])
    String effective_variants_nirvana_docker = select_first([variants_nirvana_docker, GetToolVersions.variants_nirvana_docker])
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

    String extract_output_gcs_dir = "~{effective_workspace_bucket}/output_vcfs/by_submission_id/~{effective_submission_id}/~{dataset_suffix}"

    call CreateVATFromVDS.GvsCreateVATfromVDS as CreateVATfromVDS {
        input:
            project_id = project_id,
            dataset_name = CreateDatasetForTest.dataset_name,
            ancestry_file = ancestry_path,
            filter_set_name = "quickit",
            vds_path = vds_path,
            output_path = output_path,
            split_intervals_scatter_count = split_intervals_scatter_count,

            git_branch_or_tag = git_branch_or_tag,
            basic_docker = effective_basic_docker,
            cloud_sdk_docker = effective_cloud_sdk_docker,
            gatk_docker = effective_gatk_docker,
            variants_docker = effective_variants_docker,
            variants_nirvana_docker = effective_variants_nirvana_docker,
    }

    String expected_prefix = expected_output_prefix + dataset_suffix + "/"
#    call AssertIdenticalOutputs {
#        input:
#            expected_output_prefix = expected_prefix,
#            expected_output_suffix = if (bgzip_output_vcfs) then ".bgz" else ".gz",
#            actual_vcfs = JointVariantCalling.output_vcfs,
#            gatk_docker = effective_gatk_docker
#    }
#
#    if (check_expected_cost_and_table_size_outputs) {
#        call AssertCostIsTrackedAndExpected {
#            input:
#                go = JointVariantCalling.done,
#                dataset_name = CreateDatasetForTest.dataset_name,
#                project_id = project_id,
#                expected_output_csv = expected_prefix + "cost_observability.csv",
#                cloud_sdk_docker = effective_cloud_sdk_docker,
#        }
#
#        call AssertTableSizesAreExpected {
#            input:
#                go = JointVariantCalling.done,
#                dataset_name = CreateDatasetForTest.dataset_name,
#                project_id = project_id,
#                expected_output_csv = expected_prefix + "table_sizes.csv",
#                cloud_sdk_docker = effective_cloud_sdk_docker,
#        }
#    }


    output {
#        Array[File] output_vcfs = JointVariantCalling.output_vcfs
#        Array[File] output_vcf_indexes = JointVariantCalling.output_vcf_indexes
#        Float total_vcfs_size_mb = JointVariantCalling.total_vcfs_size_mb
#        File manifest = JointVariantCalling.manifest
        String dataset_name = CreateDatasetForTest.dataset_name
        String filter_set_name = "quickit"
        String recorded_git_hash = effective_git_hash
        Boolean done = true
#        Boolean used_tighter_gcp_quotas = JointVariantCalling.used_tighter_gcp_quotas
    }
}

task AssertIdenticalOutputs {
    input {
        String expected_output_prefix
        String expected_output_suffix
        Array[File] actual_vcfs
        String gatk_docker
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
        gcloud storage cp -r "${expected_prefix}"'/*.vcf~{expected_output_suffix}' .
        gzip -S ~{expected_output_suffix} -d *~{expected_output_suffix}
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
        ls -1 | grep -E '\.vcf\~{expected_output_suffix}$' | xargs gzip -S ~{expected_output_suffix} -d
        cd ..

        echo "Header Check"
        # Headers first, these can yield useful diagnostics when there are mismatches.
        for vcf in $(ls -1 actual | grep -E '\.vcf$')
        do
          actual="actual/$vcf"
          expected="expected/$vcf"
          set +o errexit
          cmp <(grep '^#' $actual | grep -E -v '^##GATKCommandLine=') <(grep '^#' $expected | grep -E -v '^##GATKCommandLine=')
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
          cmp <(grep -E -v '^##GATKCommandLine=' $actual) <(grep -E -v '^##GATKCommandLine=' $expected)
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
        docker: gatk_docker
        disks: "local-disk 500 HDD"
    }

    output {
        Boolean done = true
    }
}
