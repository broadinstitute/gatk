version 1.0

import "GvsQuickstartVcfIntegration.wdl" as QuickstartVcfIntegration
import "GvsQuickstartHailIntegration.wdl" as QuickstartHailIntegration
import "GvsUtils.wdl" as Utils

workflow GvsQuickstartIntegration {
    input {
        String branch_name
        Boolean run_vcf_integration = true
        Boolean run_hail_integration = true
    }

    String project_id = "gvs-internal"
    File interval_list = "gs://gvs-internal-scratch/mcovarr/scratch/quicker/chr20_chrX_chrY_calling_regions.hg38.noCentromeres.noTelomeres.interval_list"

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
                interval_list = interval_list,
        }
        call QuickstartHailIntegration.GvsQuickstartHailIntegration as GvsQuickstartHailVQSRClassicIntegration{
            input:
                branch_name = branch_name,
                use_classic_VQSR = true,
                extract_do_not_filter_override = false,
                dataset_suffix = "classic_hail",
                gatk_override = BuildGATKJar.jar,
                interval_list = interval_list,
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
                interval_list = interval_list,
        }
        call QuickstartVcfIntegration.GvsQuickstartVcfIntegration as QuickstartVcfVQSRClassicIntegration {
            input:
                branch_name = branch_name,
                use_classic_VQSR = true,
                extract_do_not_filter_override = true,
                dataset_suffix = "classic_vcf",
                gatk_override = BuildGATKJar.jar,
                interval_list = interval_list,
        }
    }
}
