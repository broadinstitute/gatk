version 1.0
workflow GvsSitesOnlyVCF {
   input {
        File gvs_extract_cohort_filtered_vcf
        String output_sites_only_file_name
        File? gatk_override
    }

    call SitesOnlyVcf {
                input:
                    vcf_bgz_gts                 = gvs_extract_cohort_filtered_vcf,
                    output_filename             = "${output_sites_only_file_name}.sites_only.vcf.gz",
            }

    }

################################################################################
task SitesOnlyVcf {
    input {
        File vcf_bgz_gts
        String output_filename
    }
    String output_vcf_idx = basename(output_filename) + ".tbi" # or will this be .idx if from .vcf.gz
    command <<<
        set -e
        gatk --java-options "-Xmx2048m" \
            SelectVariants \
                -V ~{vcf_bgz_gts} \
                --exclude-filtered \
                --sites-only-vcf-output \
                -O ~{output_filename}
     >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "broadinstitute/gatk:4.2.0.0"
        memory: "3 GB"
        cpu: "1"
        disks: "local-disk 100 HDD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        File output_vcf="~{output_filename}"
        File output_vcf_idx="~{output_vcf_idx}"
    }
}


