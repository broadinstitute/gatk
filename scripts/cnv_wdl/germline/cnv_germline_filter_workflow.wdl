version 1.0

workflow FilterVCFs {
    input {
        Array[File] vcfs  # Input array of VCF files
        String filter_expression = "QUAL < 50"  # Example filter criteria
        String ref_fasta
        String gatk_docker
    }

    scatter (vcf in vcfs) {
        call FilterVCF {
            input:
                vcf_file = vcf,
                filter_expression = filter_expression,
                ref_fasta = ref_fasta,
                gatk_docker = gatk_docker
        }
    }

    output {
        Array[File] filtered_vcfs = FilterVCF.filtered_vcf
    }
}

task FilterVCF {
    input {
        File vcf_file
        String filter_expression
        String ref_fasta
        String gatk_docker
    }

    command {
        gatk --java-options "-Xmx4g" VariantFiltration \
            -R ${ref_fasta} \
            -V ${vcf_file} \
            --filter-expression "${filter_expression}" \
            --filter-name "LowQual" \
            -O filtered_${basename(vcf_file)}
    }

    output {
        File filtered_vcf = "filtered_${basename(vcf_file)}"
    }

    runtime {
        docker: gatk_docker #"biocontainers/bcftools:v1.10.2-1-deb_cv1"
        memory: "8G"
        cpu: 2
    }
}