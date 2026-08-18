version 1.0

# VS-1966: Validate the VCF headers ingested during a headers-only GVS ingest, as an early
# fail-fast sanity check of the input gVCFs *before* the expensive vet/ref data ingest.
#
# Operational model (see scripts/variantstore/docs/parquet/header_loading_design.md and
# scripts/variantstore/docs/aou/AOU_DELIVERABLES.md): run GvsBulkIngestGenomes / GvsImportGenomes
# once with load_vcf_headers=true / load_vet_and_ref_ranges=false to load only headers, then run
# this workflow to validate them, and only if it passes proceed with the data ingest. This
# automates the manual DRAGEN-version query that AoU previously ran by hand.
import "GvsUtils.wdl" as Utils

workflow GvsValidateVcfHeaders {
    input {
        # Intentionally unused: this input exists solely to enforce workflow ordering. When this
        # workflow is called as a subworkflow, a parent can pass the headers-only ingest's `done`
        # output here so validation cannot start until that ingest has completed. Defaults to true
        # for standalone runs.
        #@ except: UnusedInput
        Boolean go = true
        String dataset_name
        String project_id
        String? expected_dragen_version
        Boolean require_reblocking = true
        Boolean fail_on_validation_errors = true

        String? git_branch_or_tag
        String? variants_docker
        String? basic_docker
    }

    parameter_meta {
        dataset_name: {
            help: "BigQuery dataset holding the ingested header tables (vcf_header_lines, sample_vcf_header, sample_info)."
        }
        project_id: {
            help: "Google project for the GVS dataset."
        }
        expected_dragen_version: {
            help: "Optional expected DRAGEN version. A single triplet, e.g. '3.7.8', requires every sample to match it exactly (AoU); a range requires every sample's triplet to fall within it. Ranges accept interval notation: '3.4.12-3.7.8' (both inclusive), '[3.7.8-3.8)' (inclusive-exclusive), '(3.7-3.8)' (both exclusive). Bounds compare positionally, so an inclusive two-component upper bound is not 'all of that minor' -- '3.7.8-3.8' excludes 3.8.1; use the exclusive form '3.7.8-3.8)' for 'everything below 3.8'. If unset, only cross-sample consistency is enforced."
        }
        require_reblocking: {
            help: "If true (default), a sample without a ReblockGVCF command line fails validation (AoU fails fast on non-reblocked input). Set to false to report non-reblocked samples without failing."
        }
        fail_on_validation_errors: {
            help: "If true (default), the workflow fails when any fatal header check fails, stopping the pipeline before the expensive data ingest. If false, the results are reported without failing the workflow."
        }
    }

    if (!defined(variants_docker) || !defined(basic_docker)) {
        call Utils.GetToolVersions {
            input:
                git_branch_or_tag = git_branch_or_tag,
        }
    }

    String effective_variants_docker = select_first([variants_docker, GetToolVersions.variants_docker])
    String effective_basic_docker = select_first([basic_docker, GetToolVersions.basic_docker])

    call ValidateVcfHeaders {
        input:
            dataset_name = dataset_name,
            project_id = project_id,
            expected_dragen_version = expected_dragen_version,
            require_reblocking = require_reblocking,
            variants_docker = effective_variants_docker,
    }

    # Fail fast: if any fatal check failed and the caller wants errors to stop the pipeline, abort
    # here (before the expensive vet/ref ingest) with the report as the failure message.
    if (!ValidateVcfHeaders.validation_passed && fail_on_validation_errors) {
        call Utils.TerminateWorkflow {
            input:
                message = "VCF header validation failed -- do NOT proceed with data ingest. Report follows:\n" + ValidateVcfHeaders.report_contents,
                basic_docker = effective_basic_docker,
        }
    }

    output {
        Boolean validation_passed = ValidateVcfHeaders.validation_passed
        File validation_report = ValidateVcfHeaders.report
        String validation_report_contents = ValidateVcfHeaders.report_contents
        Boolean done = true
    }
}

task ValidateVcfHeaders {
    input {
        String dataset_name
        String project_id
        String? expected_dragen_version
        Boolean require_reblocking = true
        String variants_docker
    }

    meta {
        # Should always run against the current tables; do not call-cache.
        volatile: true
    }

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        # The script writes 'true'/'false' to pass.txt and a human-readable report to report.txt.
        # It exits 0 even on a failed validation (the pass/fail is surfaced via pass.txt), so the
        # workflow -- not this task -- decides whether a failure should abort the pipeline.
        # expected_dragen_version is single-quoted because interval-notation values contain shell
        # metacharacters -- e.g. unquoted '(3.7-3.8)' is a bash syntax error. Cromwell does not quote
        # interpolated values, so the quotes must be inside the placeholder. The whole placeholder
        # still drops to an empty string when this optional input is undefined.
        python3 /app/check_vcf_headers.py \
            --project_id ~{project_id} \
            --dataset_name ~{dataset_name} \
            ~{"--expected_dragen_version '" + expected_dragen_version + "'"} \
            ~{if require_reblocking then "" else "--allow_non_reblocked"} \
            --pass_file_output pass.txt \
            --report_file_output report.txt
    >>>

    runtime {
        docker: variants_docker
        memory: "3 GB"
        disks: "local-disk 50 HDD"
        cpu: 1
        preemptible: 1
    }

    output {
        Boolean validation_passed = read_boolean("pass.txt")
        File report = "report.txt"
        String report_contents = read_string("report.txt")
    }
}
