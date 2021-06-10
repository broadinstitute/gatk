version 1.0

workflow GvsSitesOnlyVCF {
   input {
        Array[File] gvs_extract_cohort_filtered_vcfs
        String output_sites_only_file_name
        String output_merged_file_name
        File? gatk_override
    }

    scatter(i in range(length(gvs_extract_cohort_filtered_vcfs)) ) {
        call SitesOnlyVcf {
            input:
                vcf_bgz_gts                 = gvs_extract_cohort_filtered_vcfs[i],
                output_filename             = "${output_sites_only_file_name}_${i}.sites_only.vcf.gz",
        }
    }

    call MergeVCFs {
      input:
          input_vcfs = SitesOnlyVcf.output_vcf,
          input_vcf_indices = SitesOnlyVcf.output_vcf_index,
          output_merged_file_name = output_merged_file_name,
    }
}

################################################################################
task SitesOnlyVcf {
    input {
        File vcf_bgz_gts
        String output_filename
    }
    # this needs to be array-ized
    parameter_meta {
        vcf_bgz_gts: {localization_optional: true}
    }
    command <<<
        set -e
        gatk --java-options "-Xmx2048m" \
            SelectVariants \
                -V ~{vcf_bgz_gts} \
                --exclude-filtered \
                --sites-only-vcf \
                -O ~{output_filename}

        tabix ~{output_filename}

     >>>
     # add an indexing bit here?

    # ------------------------------------------------
    # Runtime settings:
    runtime {
        #docker:"broadinstitute/gatk:4.2.0.0"
        docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:varstore_d8a72b825eab2d979c8877448c0ca948fd9b34c7_change_to_hwe"
        memory: "3 GB"
        cpu: "1"
        disks: "local-disk 100 HDD"
    }

    # ------------------------------------------------
    # Outputs:

    output {
        File output_vcf="~{output_filename}"
        File output_vcf_index = "~{output_filename}.tbi"
    }
}

task MergeVCFs {
    input {
        Array[File] input_vcfs
        Array[File] input_vcf_indices
        String output_merged_file_name
    }
    String output_vcf = basename(output_merged_file_name) + ".vcf.gz"
    String output_vcf_idx = basename(output_vcf) + ".tbi"
    command <<<
        set -e
        gatk --java-options "-Xmx2048m" MergeVcfs -I ~{sep=' -I ' input_vcfs} -O ~{output_vcf}
    >>>
    runtime {
        docker: "broadinstitute/gatk:4.2.0.0"
        memory: "3 GB"
        cpu: "1"
        disks: "local-disk 100 HDD"
    }
    output {
        File merged_vcf = "~{output_vcf}"
        File merged_vcf_idx = "~{output_vcf_idx}"
    }
}






