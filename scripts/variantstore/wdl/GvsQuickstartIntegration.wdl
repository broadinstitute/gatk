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
        docker: "us.gcr.io/broad-dsde-methods/variantstore:2023-06-02-alpine-0772bbec9"
    }
    output {
        File out = "filtered.interval_list"
    }
}

workflow GvsQuickstartIntegration {
    input {
        String branch_name
        Boolean run_vcf_integration = true
        Boolean run_hail_integration = false
    }

    File full_interval_list = "gs://gvs-internal-scratch/mcovarr/scratch/quicker/chr20_chrX_chrY_calling_regions.hg38.noCentromeres.noTelomeres.interval_list"

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
                use_classic_VQSR = false,
                extract_do_not_filter_override = false,
                dataset_suffix = "lite_hail",
                gatk_override = BuildGATKJar.jar,
                interval_list = FilterIntervalListChromosomes.out,
        }
        call QuickstartHailIntegration.GvsQuickstartHailIntegration as GvsQuickstartHailVQSRClassicIntegration{
            input:
                branch_name = branch_name,
                use_classic_VQSR = true,
                extract_do_not_filter_override = false,
                dataset_suffix = "classic_hail",
                gatk_override = BuildGATKJar.jar,
                interval_list = FilterIntervalListChromosomes.out,
        }
    }

    if (run_vcf_integration) {
        call QuickstartVcfIntegration.GvsQuickstartVcfIntegration as QuickstartVcfVQSRLiteIntegration {
            input:
                branch_name = branch_name,
                use_classic_VQSR = false,
                extract_do_not_filter_override = false,
                dataset_suffix = "lite_vcf",
                gatk_override = BuildGATKJar.jar,
                interval_list = FilterIntervalListChromosomes.out,
        }
        call QuickstartVcfIntegration.GvsQuickstartVcfIntegration as QuickstartVcfVQSRClassicIntegration {
            input:
                branch_name = branch_name,
                use_classic_VQSR = true,
                extract_do_not_filter_override = true,
                dataset_suffix = "classic_vcf",
                gatk_override = BuildGATKJar.jar,
                interval_list = FilterIntervalListChromosomes.out,
        }
    }
}
