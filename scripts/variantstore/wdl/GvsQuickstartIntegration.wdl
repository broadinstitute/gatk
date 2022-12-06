version 1.0

import "GvsQuickstartVcfIntegration.wdl" as QuickstartVcfIntegration
import "GvsQuickstartHailIntegration.wdl" as QuickstartHailIntegration

workflow GvsQuickstartIntegration {
    input {
        String branch_name
        Boolean run_vcf_integration = true
        Boolean run_hail_integration = true
        String hail_wheel = "gs://gvs-internal-scratch/hail-wheels/2022-10-18/0.2.102-964bee061eb0/hail-0.2.102-py3-none-any.whl"
    }

    String project_id = "gvs-internal"

    if (run_hail_integration) {
        call QuickstartHailIntegration.GvsQuickstartHailIntegration {
            input:
                branch_name = branch_name,
                hail_wheel = hail_wheel,
        }
    }

    if (run_vcf_integration) {
        call QuickstartVcfIntegration.GvsQuickstartVcfIntegration {
            input:
                branch_name = branch_name,
        }
    }
}
