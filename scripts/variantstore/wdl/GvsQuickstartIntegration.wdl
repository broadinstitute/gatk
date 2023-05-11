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
        call QuickstartHailIntegration.GvsQuickstartHailIntegration {
            input:
                branch_name = branch_name
        }
    }

    if (run_vcf_integration) {
        call QuickstartVcfIntegration.GvsQuickstartVcfIntegration {
            input:
                branch_name = branch_name,
        }
    }
}
