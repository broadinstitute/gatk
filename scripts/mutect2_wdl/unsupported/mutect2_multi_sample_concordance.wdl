#Inputs are same as mutect2_multi_sample.wdl with the addtion of a truth_list input, which is a tsv of the form
# truth1.vcf    truth1.vcf.idx
# truth2.vcf    truth2.vcf.idx
# . . .
# The rows of this file correspond to the rows of the pair_list input

import "mutect2_multi_sample.wdl" as m2_multi

workflow Mutect2_Multi_Concordance {
	File? gatk_override
	File picard
	String gatk_docker
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
    Int? preemptible_attempts
    Array[String] artifact_modes

    String? m2_args
    String? m2_filtering_args

    File truth_list
    Array[Array[String]] truth = read_tsv(truth_list)

    call m2_multi.Mutect2_Multi {
        input:
            gatk_override = gatk_override,
            picard = picard,
            gatk_docker = gatk_docker,
            oncotator_docker = "NO_ONCOTATOR",
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
            m2_args = m2_args,
            m2_filtering_args = m2_filtering_args
    }

    scatter (n in range(length(truth))) {
        call Concordance {
            input:
                gatk_override = gatk_override,
                gatk_docker = gatk_docker,
                intervals = intervals,
                truth_vcf = truth[n][0],
                truth_vcf_idx = truth[n][1],
                eval_vcf = Mutect2_Multi.filtered_vcf_files[n],
                eval_vcf_idx = Mutect2_Multi.filtered_vcf_index_files[n],
                preemptible_attempts = preemptible_attempts
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
    File? gatk_override
    String gatk_docker
    File? intervals
    File truth_vcf
    File truth_vcf_idx
    File eval_vcf
    File eval_vcf_idx
    Int preemptible_attempts

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
        disks: "local-disk " + 400 + " HDD"
        preemptible: "${preemptible_attempts}"
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