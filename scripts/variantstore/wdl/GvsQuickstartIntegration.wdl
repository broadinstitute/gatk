version 1.0

import "GvsQuickstartVcfIntegration.wdl" as QuickstartVcfIntegration
import "GvsQuickstartHailIntegration.wdl" as QuickstartHailIntegration

workflow GvsQuickstartIntegration {
    input {
        String branch_name
        Boolean run_vcf_integration = true
        Boolean run_hail_integration = true
    }

    String project_id = "gvs-internal"

    if (run_hail_integration) {
        call QuickstartHailIntegration.GvsQuickstartHailIntegration as GvsQuickstartHailVQSRLiteIntegration {
            input:
                branch_name = branch_name,
                use_classic_VQSR = false,
                dataset_suffix = "lite.hail"
        }
        call QuickstartHailIntegration.GvsQuickstartHailIntegration as GvsQuickstartHailVQSRClassicIntegration{
            input:
                branch_name = branch_name,
                use_classic_VQSR = true,
                dataset_suffix = "classic.hail"
        }
    }

    if (run_vcf_integration) {
        call QuickstartVcfIntegration.GvsQuickstartVcfIntegration as QuickstartVcfVQSRLiteIntegration {
            input:
                branch_name = branch_name,
                use_classic_VQSR = false,
                dataset_suffix = "lite.vcf"
        }
        call QuickstartVcfIntegration.GvsQuickstartVcfIntegration as QuickstartVcfVQSRClassicIntegration {
            input:
                branch_name = branch_name,
                use_classic_VQSR = true,
                dataset_suffix = "classic.vcf"
        }
    }
}
