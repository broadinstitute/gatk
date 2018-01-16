#Inputs are same as mutect2_multi_sample.wdl with the addtion of a truth_list input, which is a tsv of the form
# truth1.vcf    truth1.vcf.idx
# truth2.vcf    truth2.vcf.idx
# . . .
# The rows of this file correspond to the rows of the pair_list input

import "mutect2_multi_sample.wdl" as m2_multi

workflow Mutect2_Multi_Concordance {
    # inputs
	File? intervals
    File ref_fasta
    File ref_fai
    File ref_dict
	Int scatter_count
	File pair_list
	File? pon
	File? pon_index
	File? gnomad
	File? gnomad_index
	File? variants_for_contamination
    File? variants_for_contamination_index
	Boolean is_run_orientation_bias_filter
	Array[String] artifact_modes
    String? m2_extra_args
    String? m2_extra_filtering_args
    File truth_list
    Array[Array[String]] truth = read_tsv(truth_list)

    File? gatk_override

    # runtime
    String gatk_docker
    Int? preemptible_attempts

    call m2_multi.Mutect2_Multi {
        input:
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
        	scatter_count = scatter_count,
        	pair_list = pair_list,
        	intervals = intervals,
        	pon = pon,
        	pon_index = pon_index,
        	gnomad = gnomad,
        	gnomad_index = gnomad_index,
        	variants_for_contamination = variants_for_contamination,
            variants_for_contamination_index = variants_for_contamination_index,
        	is_run_orientation_bias_filter = is_run_orientation_bias_filter,
        	is_run_oncotator = false,
            preemptible_attempts = preemptible_attempts,
            artifact_modes = artifact_modes,
            m2_extra_args = m2_extra_args,
            m2_extra_filtering_args = m2_extra_filtering_args,
            gatk_override = gatk_override,
            gatk_docker = gatk_docker,
            oncotator_docker = "NO_ONCOTATOR"
    }

    scatter (n in range(length(truth))) {
        call Concordance {
            input:
                intervals = intervals,
                truth_vcf = truth[n][0],
                truth_vcf_idx = truth[n][1],
                eval_vcf = Mutect2_Multi.filtered_vcf[n],
                eval_vcf_idx = Mutect2_Multi.filtered_vcf_idx[n],
                preemptible_attempts = preemptible_attempts,
                gatk_override = gatk_override,
                gatk_docker = gatk_docker
        }
    }

    output {
        Array[File] tpfn = Concordance.tpfn
        Array[File] tpfn_idx = Concordance.tpfn_idx
        Array[File] tpfp = Concordance.tpfp
        Array[File] tpfp_idx = Concordance.tpfp_idx
        Array[File] ftnfn = Concordance.ftnfn
        Array[File] ftnfn_idx = Concordance.ftnfn_idx
        Array[File] summary = Concordance.summary
    }
}

task Concordance {
    # inputs
    File? intervals
    File truth_vcf
    File truth_vcf_idx
    File eval_vcf
    File eval_vcf_idx

    File? gatk_override

    # runtime
    String gatk_docker
    Int? preemptible_attempts

    command {
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        gatk --java-options "-Xmx2g" Concordance \
            ${"-L " + intervals} \
            -truth ${truth_vcf} -eval ${eval_vcf} \
            -tpfn "tpfn.vcf" \
            -tpfp "tpfp.vcf" \
            -ftnfn "ftnfn.vcf" \
            -summary summary.tsv
    }

    runtime {
        memory: "5 GB"
        docker: "${gatk_docker}"
        disks: "local-disk " + 100 + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        File tpfn = "tpfn.vcf"
        File tpfn_idx = "tpfn.vcf.idx"
        File tpfp = "tpfp.vcf"
        File tpfp_idx = "tpfp.vcf.idx"
        File ftnfn = "ftnfn.vcf"
        File ftnfn_idx = "ftnfn.vcf.idx"
        File summary = "summary.tsv"
    }
}