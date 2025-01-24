version 1.0

import "../GvsUtils.wdl" as Utils
import "../../variant-annotations-table/GvsCreateVATfromVDS.wdl" as CreateVATFromVDS
import "../../variant-annotations-table/GvsValidateVAT.wdl" as ValidateVAT

# A comment for debugging.

workflow GvsQuickstartVATIntegration {
    input {
        String git_branch_or_tag
        String? git_hash
        Boolean use_vds = true      # If true, use a VDS, otherwise use a sites only VCF.
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
        String? submission_id
    }
    String project_id = "gvs-internal"

    File input_data_prefix = "gs://gvs-internal-quickstart/integration/test_data/2025-01-17/"
    File ancestry_path =  input_data_prefix + "quickstart_ancestry.tsv"
    File? vds_path = if (use_vds) then input_data_prefix + "gvs_export.vds" else none
    File? sites_only_vcf = if (!use_vds) then input_data_prefix + "quickstart_sites_only.vcf.bgz" else none

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

    call CreateVATFromVDS.GvsCreateVATfromVDS as CreateVATfromVDS {
        input:
            project_id = project_id,
            dataset_name = CreateDatasetForTest.dataset_name,
            ancestry_file = ancestry_path,
            filter_set_name = "quickit",
            vds_path = vds_path,
            sites_only_vcf = sites_only_vcf,
            output_path = output_path,
            split_intervals_scatter_count = split_intervals_scatter_count,

            git_branch_or_tag = git_branch_or_tag,
            basic_docker = effective_basic_docker,
            cloud_sdk_docker = effective_cloud_sdk_docker,
            gatk_docker = effective_gatk_docker,
            variants_docker = effective_variants_docker,
            variants_nirvana_docker = effective_variants_nirvana_docker,
    }

    call ValidateVAT.GvsValidateVat {
        input:
            project_id = project_id,
            dataset_name = CreateDatasetForTest.dataset_name,
            vat_table_name = CreateVATfromVDS.vat_table_name,
            is_small_callset = true,
            cloud_sdk_docker = effective_cloud_sdk_docker,
            variants_docker = effective_variants_docker,
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
    call AssertTableSizeIsAsExpected {
        input:
            dataset_name = CreateDatasetForTest.dataset_name,
            project_id = project_id,
            vat_table_name = CreateVATfromVDS.vat_table_name,
            expected_output_csv = expected_prefix + "table_sizes.csv",
            cloud_sdk_docker = effective_cloud_sdk_docker,
    }

    output {

        String dataset_name = CreateDatasetForTest.dataset_name
        String filter_set_name = "quickit"
        String recorded_git_hash = effective_git_hash
        Boolean done = true
    }
}

#task AssertIdenticalOutputs {
#    input {
#        String expected_output_prefix
#        String expected_output_suffix
#        Array[File] actual_vcfs
#        String gatk_docker
#    }
#    parameter_meta {
#        actual_vcfs: {
#            localization_optional: true
#        }
#    }
#    command <<<
#        # Prepend date, time and pwd to xtrace log entries.
#        PS4='\D{+%F %T} \w $ '
#        set -o errexit -o nounset -o pipefail -o xtrace
#
#        failures=()
#
#        # Where the current set of expected results lives in the cloud
#        expected_prefix="~{expected_output_prefix}"
#        # Remove a trailing slash if there is one
#        expected_prefix=${expected_prefix%/}
#
#        # Download all the expected data
#        mkdir expected
#        cd expected
#        gcloud storage cp -r "${expected_prefix}"'/*.vcf~{expected_output_suffix}' .
#        gzip -S ~{expected_output_suffix} -d *~{expected_output_suffix}
#        cd ..
#
#        mkdir actual
#        cd actual
#        touch actual_manifest.txt
#        # Making the manifest is pretty uninteresting and very noisy so turn off xtrace temporarily.
#        set +o xtrace
#        for actual in ~{sep=' ' actual_vcfs}
#        do
#            echo $actual >> actual_manifest.txt
#        done
#        set -o xtrace
#
#        cat actual_manifest.txt | gcloud storage cp -I .
#        # Unzip actual result data.
#        ls -1 | grep -E '\.vcf\~{expected_output_suffix}$' | xargs gzip -S ~{expected_output_suffix} -d
#        cd ..
#
#        echo "Header Check"
#        # Headers first, these can yield useful diagnostics when there are mismatches.
#        for vcf in $(ls -1 actual | grep -E '\.vcf$')
#        do
#          actual="actual/$vcf"
#          expected="expected/$vcf"
#          set +o errexit
#          cmp <(grep '^#' $actual | grep -E -v '^##GATKCommandLine=') <(grep '^#' $expected | grep -E -v '^##GATKCommandLine=')
#          rc=$?
#          set -o errexit
#          if [[ $rc -ne 0 ]]; then
#            # If there is a mismatch add it to a list of failures but keep on looking for mismatches.
#            failures+=( $vcf )
#          fi
#        done
#
#        echo "Header Failure Check"
#        if [[ ${#failures[@]} -ne 0 ]]; then
#          echo "Error: headers for the following files do not match:"
#          for failure in ${failures[@]}; do
#            echo $failure
#            expected="expected/$failure"
#            actual="actual/$failure"
#            diff <(grep '^#' $actual) <(grep '^#' $expected)
#          done
#          exit 1
#        fi
#
#        echo "Overall Check"
#        # If the headers all matched look for any mismatches in overall file content.
#        fail=0
#        for vcf in $(ls -1 actual | grep -E '\.vcf$')
#        do
#          expected="expected/$vcf"
#          actual="actual/$vcf"
#          set +o errexit
#          cmp <(grep -E -v '^##GATKCommandLine=' $actual) <(grep -E -v '^##GATKCommandLine=' $expected)
#          rc=$?
#          set -o errexit
#          if [[ $rc -ne 0 ]]; then
#            echo "Error: file contents of expected and actual do not match: $vcf"
#            fail=1
#          fi
#        done
#
#        if [[ $fail -ne 0 ]]; then
#          exit 1
#        fi
#
#        echo "All vcfs compared and matched!"
#    >>>
#
#    runtime {
#        docker: gatk_docker
#        disks: "local-disk 500 HDD"
#    }
#
#    output {
#        Boolean done = true
#    }
#}

task AssertTableSizeIsAsExpected {
    meta {
        # we want to check the database each time this runs
        volatile: true
    }

    input {
        String dataset_name
        String project_id
        String vat_table_name
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
            "SELECT 'vat_total' AS total_name, sum(total_billable_bytes) AS total_bytes \
            FROM \`~{dataset_name}.INFORMATION_SCHEMA.PARTITIONS\` \
            WHERE table_name = '~{vat_table_name}'" > output/table_sizes.csv

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

