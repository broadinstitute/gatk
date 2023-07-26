version 1.0

import "GvsQuickstartVcfIntegration.wdl" as QuickstartVcfIntegration
import "GvsQuickstartHailIntegration.wdl" as QuickstartHailIntegration
import "GvsUtils.wdl" as Utils

# .....x

workflow GvsQuickstartIntegration {
    input {
        String branch_name
        Boolean use_default_dockers = false
        Boolean run_vcf_integration = true
        Boolean run_hail_integration = true
        Boolean run_exome_integration = true
    }

    File full_wgs_interval_list = "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list"
    File full_exome_interval_list = "gs://gcp-public-data--broad-references/hg38/v0/bge_exome_calling_regions.v1.1.interval_list"
    File expected_output_prefix = "gs://gvs-internal-quickstart/integration/2023-07-25-quicker/"

    # WDL 1.0 trick to set a variable ('none') to be undefined.
    if (false) {
        File? none = ""
    }

    call FilterIntervalListChromosomes {
        input:
            full_interval_list = full_wgs_interval_list,
            chromosomes = ["chrX", "chrY", "chr20"],
    }

    if (!use_default_dockers) {
        call Utils.BuildGATKJar {
            input:
                branch_name = branch_name,
        }
    }

    if (run_hail_integration) {
        call QuickstartHailIntegration.GvsQuickstartHailIntegration as GvsQuickstartHailVQSRLiteIntegration {
            input:
                branch_name = branch_name,
                is_wgs = true,
                use_VQSR_lite = true,
                extract_do_not_filter_override = false,
                dataset_suffix = "lite_hail",
                use_default_dockers = use_default_dockers,
                gatk_override = if (use_default_dockers) then none else BuildGATKJar.jar,
                interval_list = FilterIntervalListChromosomes.out,
                expected_output_prefix = expected_output_prefix,
        }
        call QuickstartHailIntegration.GvsQuickstartHailIntegration as GvsQuickstartHailVQSRClassicIntegration {
            input:
                branch_name = branch_name,
                is_wgs = true,
                use_VQSR_lite = false,
                extract_do_not_filter_override = false,
                dataset_suffix = "classic_hail",
                use_default_dockers = use_default_dockers,
                gatk_override = if (use_default_dockers) then none else BuildGATKJar.jar,
                interval_list = FilterIntervalListChromosomes.out,
                expected_output_prefix = expected_output_prefix,
        }
    }

    if (run_vcf_integration) {
        call QuickstartVcfIntegration.GvsQuickstartVcfIntegration as QuickstartVcfVQSRLiteIntegration {
            input:
                branch_name = branch_name,
                is_wgs = true,
                use_VQSR_lite = true,
                extract_do_not_filter_override = false,
                dataset_suffix = "lite_vcf",
                use_default_dockers = use_default_dockers,
                gatk_override = if (use_default_dockers) then none else BuildGATKJar.jar,
                interval_list = FilterIntervalListChromosomes.out,
                expected_output_prefix = expected_output_prefix,
        }
        call QuickstartVcfIntegration.GvsQuickstartVcfIntegration as QuickstartVcfVQSRClassicIntegration {
            input:
                branch_name = branch_name,
                is_wgs = true,
                use_VQSR_lite = false,
                extract_do_not_filter_override = true,
                dataset_suffix = "classic_vcf",
                use_default_dockers = use_default_dockers,
                gatk_override = if (use_default_dockers) then none else BuildGATKJar.jar,
                interval_list = FilterIntervalListChromosomes.out,
                expected_output_prefix = expected_output_prefix,
        }
    }

    if (run_exome_integration) {
        call QuickstartVcfIntegration.GvsQuickstartVcfIntegration as QuickstartVcfExomeIntegration {
            input:
                branch_name = branch_name,
                is_wgs = false,
                use_VQSR_lite = true,
                extract_do_not_filter_override = false,
                dataset_suffix = "exome",
                use_default_dockers = use_default_dockers,
                gatk_override = if (use_default_dockers) then none else BuildGATKJar.jar,
                interval_list = full_exome_interval_list,
                expected_output_prefix = expected_output_prefix,
        }
    }
}

task FilterIntervalListChromosomes {
    input {
        File full_interval_list
        Array[String]+ chromosomes
    }
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        python3 /app/filter_interval_list_chromosomes.py --input-interval-list ~{full_interval_list} \
        --output-interval-list "filtered.interval_list" --chromosome ~{sep=' --chromosome ' chromosomes}
    >>>
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:2023-06-23-alpine"
    }
    output {
        File out = "filtered.interval_list"
    }
}
