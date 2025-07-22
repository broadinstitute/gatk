version 1.0

import "../GvsUtils.wdl" as GvsUtils
import "GvsQuickstartVCFIntegration.wdl" as Integration

workflow AssertIdenticalMergedOutputs {
    input {
        File actual
        String expected_output_prefix
    }

    meta {
        description: "Asserts that two merged VCFs are identical."
    }

    call GvsUtils.GetToolVersions {}

    call Integration.AssertIdenticalMergedOutputs as DoIt {
        input:
            merged_output_vcf = actual,
            expected_output_prefix = expected_output_prefix,
            variants_docker = GetToolVersions.variants_docker,
    }
}