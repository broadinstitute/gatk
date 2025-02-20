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

        Float overlap_thresh = 0.5

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
            padding = padding,
            overlap_thresh = overlap_thresh
    }

    if (SingleSampleGCNVAndFilterVCFs.qc_passed) {
        call ReformatGCNVForFabric {
            input:
                cnv_vcf = SingleSampleGCNVAndFilterVCFs.filtered_vcf
        }

        call MergeVcfs {
            input:
                cnv_vcf = ReformatGCNVForFabric.reformatted_vcf,
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
        File merged_vcf_md5sum = select_first([MergeVcfs.merged_vcf_md5sum, short_variant_vcf_md5sum])

        Boolean qc_passed = SingleSampleGCNVAndFilterVCFs.qc_passed
        File cnv_metrics = SingleSampleGCNVAndFilterVCFs.cnv_metrics

    }
}

#Fabric doesn't seem to like ./. genotypes on
task ReformatGCNVForFabric {
    input {
        File cnv_vcf
        Int disk_size_gb = 20
        Int mem_gb = 4
    }

    String output_basename = basename(cnv_vcf, ".filtered.genotyped-segments.vcf.gz")

    command <<<
        set -euo pipefail

        python << CODE
        from pysam import VariantFile

        with VariantFile("~{cnv_vcf}") as cnv_vcf_in:
            header_out = cnv_vcf_in.header
            header_out.info.add("CN", "A", "Integer", "Copy number associated with <CNV> alleles")
            with VariantFile("~{output_basename}.reformatted_for_fabric.vcf.gz",'w', header = header_out) as cnv_vcf_out:
                for rec in cnv_vcf_in.fetch():
                    if 'PASS' in rec.filter:
                        if rec.alts and rec.alts[0] == "<DUP>":
                            for rec_sample in rec.samples.values():
                                ploidy = len(rec_sample.alleles)
                                rec_sample.allele_indices = (None,)*(ploidy - 1) + (1,)
                        cnv_vcf_out.write(rec)
        CODE
    >>>

    runtime {
            docker: "us.gcr.io/broad-dsde-methods/pysam:v1.1"
            preemptible: 3
            cpu: 2
            disks: "local-disk " + disk_size_gb + " HDD"
            memory: mem_gb + " GB"
        }

    output {
        File reformatted_vcf = "~{output_basename}.reformatted_for_fabric.vcf.gz"
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
    String output_cnv_basename = basename(cnv_vcf, ".reformatted_for_fabric.vcf.gz")
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