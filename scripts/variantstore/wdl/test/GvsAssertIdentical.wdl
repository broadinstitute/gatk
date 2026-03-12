version 1.0

workflow GvsAssertIdentical {
    input {
        File merged_output_vcf
        String expected_output_prefix
        String variants_docker
    }

    call AssertIdenticalMergedOutputs {
        input:
            merged_output_vcf = merged_output_vcf,
            expected_output_prefix = expected_output_prefix,
            variants_docker = variants_docker,
    }
}


task AssertIdenticalMergedOutputs {
    input {
        File merged_output_vcf
        String expected_output_prefix
        String variants_docker
    }

    parameter_meta {
        merged_output_vcf: {
            localization_optional: true
        }
    }
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        # Where the current set of expected results lives in the cloud
        expected_prefix="~{expected_output_prefix}"
        # Remove a trailing slash if there is one
        expected_prefix=${expected_prefix%/}

        # Specify a GCS auth token for htslib/bcftools. Note this token only works for a limited time, but that's fine
        # for this use case. If a persistent solution is required, see
        # https://broadinstitute.slack.com/archives/C0544AAC70D/p1696360070640409

        # Let's not log the token.
        set +o xtrace
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        set -o xtrace

        bcftools view --header-only ~{merged_output_vcf} | grep -E -v '^##GATKCommandLine=|^##bcftools' > actual_header.txt
        bcftools view --header-only ${expected_prefix}/merged.vcf.gz | grep -E -v '^##GATKCommandLine=|^##bcftools' > expected_header.txt

        echo "Header Failure Check"
        diff actual_header.txt expected_header.txt

        bcftools view --no-header ~{merged_output_vcf} > actual_content.txt
        bcftools view --no-header ${expected_prefix}/merged.vcf.gz > expected_content.txt

        echo "Content Failure Check"
        diff actual_content.txt expected_content.txt

        echo "VCFs compared and matched!"
    >>>

    runtime {
        docker: variants_docker
        disks: "local-disk 500 HDD"
    }

    output {
        Boolean done = true
    }
}