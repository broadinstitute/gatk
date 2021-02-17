version 1.0

workflow ReblockGVCF {
    input {
        File gvcf
        File ref_dict
        File ref_fasta
        File ref_fasta_index
        Int disk_size
    }

    String sub_strip_path = "gs://.*/"
    String sub_strip_gvcf = ".g.vcf.gz" + "$"
    String sub_sub = sub(sub(gvcf, sub_strip_path, ""), sub_strip_gvcf, "")
    String gatk_docker = "us.gcr.io/broad-dsde-methods/update_reblocking@sha256:3f7fe2a45656260c5a678f3f0e77174d6ec0305056e5fc05e62a3ca38b1f828d"

    call reblock {
        input:
            gvcf = gvcf,
            gvcf_index = gvcf + ".tbi",
            ref_dict = ref_dict,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            output_vcf_filename = sub_sub + ".vcf.gz",
            disk_size = disk_size,
            gatk_docker = gatk_docker
    }

    output {
        File reblocked_vcf = reblock.output_vcf
        File reblocked_vcf_index = reblock.output_vcf_index
    }
}

task reblock {
    input {
        File gvcf
        File gvcf_index

        File ref_dict
        File ref_fasta
        File ref_fasta_index

        String output_vcf_filename

        Int disk_size
        String gatk_docker
    }

    parameter_meta {
        gvcf: {localization_optional: true},
        gvcf_index: {localization_optional: true},
        ref_fasta: {localization_optional: true},
        ref_fasta_index: {localization_optional: true}
    }

    command <<<

        gatk --java-options "-Xms3g -Xmx3g" \
        ReblockGVCF \
        -V ${gvcf} \
        -R ${ref_fasta} \
        -do-qual-approx \
        --floor-blocks -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 \
        -O ${output_vcf_filename}

    >>>
    runtime {
        memory: "3 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 3
        docker: gatk_docker
    }
    output {
        File output_vcf = "${output_vcf_filename}"
        File output_vcf_index = "${output_vcf_filename}.tbi"
    }
}
