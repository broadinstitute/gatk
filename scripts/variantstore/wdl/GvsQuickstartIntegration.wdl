version 1.0

import "GvsQuickstartVcfIntegration.wdl" as QuickstartVcfIntegration
import "GvsQuickstartHailIntegration.wdl" as QuickstartHailIntegration
import "GvsUtils.wdl" as Utils

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
        docker: "us.gcr.io/broad-dsde-methods/variantstore:2023-07-24-alpine-bd6b9d62e"
    }
    output {
        File out = "filtered.interval_list"
    }
}

workflow GvsQuickstartIntegration {
    input {
        String branch_name
        Boolean run_vcf_integration = true
        Boolean run_hail_integration = true
        String? sample_id_column_name ## Note that a column WILL exist that is the <entity>_id from the table name. However, some users will want to specify an alternate column for the sample_name during ingest
        String? vcf_files_column_name
        String? vcf_index_files_column_name
        String? sample_set_name ## NOTE: currently we only allow the loading of one sample set at a time
    }

    File full_interval_list = "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list"
    File expected_output_prefix = "gs://gvs-internal-quickstart/integration/2023-07-19-quicker/"

    call FilterIntervalListChromosomes {
        input:
            full_interval_list = full_interval_list,
            chromosomes = ["chrX", "chrY", "chr20"],
    }

    call Utils.BuildGATKJar {
        input:
            branch_name = branch_name,
    }

    if (run_hail_integration) {
        call QuickstartHailIntegration.GvsQuickstartHailIntegration as GvsQuickstartHailVQSRLiteIntegration {
            input:
                branch_name = branch_name,
                use_VQSR_lite = true,
                extract_do_not_filter_override = false,
                dataset_suffix = "lite_hail",
                gatk_override = BuildGATKJar.jar,
                interval_list = FilterIntervalListChromosomes.out,
                expected_output_prefix = expected_output_prefix,
                sample_id_column_name = sample_id_column_name,
                vcf_files_column_name = vcf_files_column_name,
                vcf_index_files_column_name = vcf_index_files_column_name,
                sample_set_name = sample_set_name,
                use_classic_VQSR = false,
        }
        call QuickstartHailIntegration.GvsQuickstartHailIntegration as GvsQuickstartHailVQSRClassicIntegration {
            input:
                branch_name = branch_name,
                use_VQSR_lite = false,
                extract_do_not_filter_override = false,
                dataset_suffix = "classic_hail",
                gatk_override = BuildGATKJar.jar,
                gatk_override = BuildGATKJar.jar,
                interval_list = FilterIntervalListChromosomes.out,
                expected_output_prefix = expected_output_prefix,
                sample_id_column_name = sample_id_column_name,
                vcf_files_column_name = vcf_files_column_name,
                vcf_index_files_column_name = vcf_index_files_column_name,
                sample_set_name = sample_set_name,
                use_classic_VQSR = true,
        }
    }

    if (run_vcf_integration) {
        call QuickstartVcfIntegration.GvsQuickstartVcfIntegration as QuickstartVcfVQSRLiteIntegration {
            input:
                branch_name = branch_name,
                use_VQSR_lite = true,
                extract_do_not_filter_override = false,
                dataset_suffix = "lite_vcf",
                gatk_override = BuildGATKJar.jar,
                interval_list = FilterIntervalListChromosomes.out,
                expected_output_prefix = expected_output_prefix,
                sample_id_column_name = sample_id_column_name,
                vcf_files_column_name = vcf_files_column_name,
                vcf_index_files_column_name = vcf_index_files_column_name,
                sample_set_name = sample_set_name,
                drop_state = "FORTY",
        }
        call QuickstartVcfIntegration.GvsQuickstartVcfIntegration as QuickstartVcfVQSRClassicIntegration {
            input:
                branch_name = branch_name,
                use_VQSR_lite = false,
                extract_do_not_filter_override = true,
                dataset_suffix = "classic_vcf",
                gatk_override = BuildGATKJar.jar,
                interval_list = FilterIntervalListChromosomes.out,
                expected_output_prefix = expected_output_prefix,
                sample_id_column_name = sample_id_column_name,
                vcf_files_column_name = vcf_files_column_name,
                vcf_index_files_column_name = vcf_index_files_column_name,
                sample_set_name = sample_set_name,
                drop_state = "FORTY",
        }
    }
}
