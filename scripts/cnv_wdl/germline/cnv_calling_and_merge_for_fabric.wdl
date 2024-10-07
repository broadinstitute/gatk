version 1.0

import "cnv_germline_case_filter_workflow.wdl" as cnv_case_and_filter
workflow CNVCallingAndMergeForFabric {
    input {
        File normal_bam
        File normal_bai
        File short_variant_vcf
        File short_variant_vcf_index
        File short_variant_vcf_md5sum

        File contig_ploidy_model_tar
        File filtered_intervals
        Array[File] gcnv_model_tars
        Array[File] gcnv_panel_genotyped_segments

        String gatk_docker
        File intervals

        Int maximum_number_events_per_sample
        Int maximum_number_pass_events_per_sample
        Int num_intervals_per_scatter
        Int ref_copy_number_autosomal_contigs

        File ref_fasta
        File ref_fasta_fai
        File ref_fasta_dict

        Array[String] allosomal_contigs
        Int padding
    }

    call cnv_case_and_filter.SingleSampleGCNVAndFilterVCFs {
        input:
            normal_bam = normal_bam,
            normal_bai = normal_bai,
            contig_ploidy_model_tar = contig_ploidy_model_tar,
            filtered_intervals = filtered_intervals,
            gcnv_model_tars = gcnv_model_tars,
            pon_genotyped_segments_vcfs = gcnv_panel_genotyped_segments,
            gatk_docker = gatk_docker,
            intervals = intervals,
            maximum_number_events_per_sample = maximum_number_events_per_sample,
            maximum_number_pass_events_per_sample = maximum_number_pass_events_per_sample,
            num_intervals_per_scatter = num_intervals_per_scatter,
            ref_copy_number_autosomal_contigs = ref_copy_number_autosomal_contigs,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            ref_fasta_dict = ref_fasta_dict,
            allosomal_contigs = allosomal_contigs,
            padding = padding
    }

    if (SingleSampleGCNVAndFilterVCFs.qc_passed) {

        call MergeVcfs {
            input:
                cnv_vcf = SingleSampleGCNVAndFilterVCFs.filtered_vcf,
                short_variant_vcf = short_variant_vcf,
                gatk_docker = gatk_docker
        }
    }

    output {
        File filtered_cnv_genotyped_segments_vcf = SingleSampleGCNVAndFilterVCFs.filtered_vcf
        File filtered_cnv_genotyped_segments_vcf_index = SingleSampleGCNVAndFilterVCFs.filtered_vcf_index
        File filtered_cnv_genotyped_segments_vcf_md5sum = SingleSampleGCNVAndFilterVCFs.filtered_vcf_md5sum

        File merged_vcf = select_first([MergeVcfs.merged_vcf, short_variant_vcf])
        File merged_vcf_index = select_first([MergeVcfs.merged_vcf_index, short_variant_vcf_index])
        File merged_vcf_md5sum = select_first([MergeVcfs.merged_vcf_md5sum, short_variant_vcf])

        Boolean qc_passed = SingleSampleGCNVAndFilterVCFs.qc_passed

    }
}

task MergeVcfs {
    input {
        File cnv_vcf
        File short_variant_vcf

        String gatk_docker
        Int mem_gb=4
        Int disk_size_gb = 100
    }

    String output_basename = basename(short_variant_vcf, ".hard-filtered.vcf.gz")
    String output_cnv_basename = basename(cnv_vcf, ".filtered.genotyped-segments.vcf.gz")
    command <<<
        set -euo pipefail

        if [ "~{output_basename}" != "~{output_cnv_basename}" ]; then
            echo "input vcf names do not agree"
            exit 1
        fi

        gatk --java-options "-Dsamjdk.create_md5=true" MergeVcfs -I ~{short_variant_vcf} -I ~{cnv_vcf} -O ~{output_basename}.merged.vcf.gz

        mv ~{output_basename}.merged.vcf.gz.md5 ~{output_basename}.merged.vcf.gz.md5sum
    >>>

    output {
        File merged_vcf = "~{output_basename}.merged.vcf.gz"
        File merged_vcf_index = "~{output_basename}.merged.vcf.gz.tbi"
        File merged_vcf_md5sum = "~{output_basename}.merged.vcf.gz.md5sum"
    }

     runtime {
        docker: gatk_docker
        memory: mem_gb + " GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 5
    }
}