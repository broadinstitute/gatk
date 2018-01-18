# Given trios consisting of: a normal sample, a relatively reliable tumor sample, and a questionable tumor sample
# run Mutect in tumor-normal mode on the "good" sample, then run it on the "bad" sample, and compare results
# The idea is to see how well Mutect performs under difficult conditions.  For example, the "good" sample might
# be from a solid tumor while the "bad" sample might be low-allele-fraction cfDNA

# Exept for the different format of the input bams table, the inputs for this script are identical to those of
# mutect2_multi_sample.wdl
import "mutect2.wdl" as m2

workflow Mutect2Trio {
	File? intervals
	File ref_fasta
	File ref_fai
	File ref_dict
	Int scatter_count
	# trio_list file is a tsv file with the following nine columns in this order.
	# normal_bam, normal_bai, good_tumor_bam, good_tumor_bai, bad_tumor_bam, bad_tumor_bai
	File trio_list
	Array[Array[String]] trios = read_tsv(trio_list)
	File? pon
	File? pon_index
	File? gnomad
	File? gnomad_index
	Boolean is_run_orientation_bias_filter
	Array[String] artifact_modes

    File? gatk_override

    # runtime
    String gatk_docker
    Int? preemptible_attempts

	scatter(trio in trios) {
		call m2.Mutect2 as GoodTumor {
			input:
				intervals = intervals,
				ref_fasta = ref_fasta,
				ref_fai = ref_fai,
				ref_dict = ref_dict,
				tumor_bam = trio[2],
				tumor_bai = trio[3],
				normal_bam = trio[0],
				normal_bai = trio[1],
				pon = pon,
				pon_index = pon_index,
				scatter_count = scatter_count,
				gnomad = gnomad,
				gnomad_index = gnomad_index,
                is_run_orientation_bias_filter = is_run_orientation_bias_filter,
                is_run_oncotator = false,
                artifact_modes = artifact_modes,
                gatk_override = gatk_override,
                gatk_docker = gatk_docker,
                oncotator_docker = "NO_ONCOTATOR",
                preemptible_attempts = preemptible_attempts
		}

		call m2.Mutect2 as BadTumor {
            input:
        		intervals = intervals,
        		ref_fasta = ref_fasta,
       			ref_fai = ref_fai,
        		ref_dict = ref_dict,
        		tumor_bam = trio[4],
        		tumor_bai = trio[5],
        		normal_bam = trio[0],
        		normal_bai=trio[1],
        		pon = pon,
        		pon_index = pon_index,
        		scatter_count = scatter_count,
        		gnomad = gnomad,
        		gnomad_index = gnomad_index,
                is_run_orientation_bias_filter = is_run_orientation_bias_filter,
                is_run_oncotator = false,
                artifact_modes = artifact_modes,
                gatk_override = gatk_override,
        	    gatk_docker = gatk_docker,
        	    oncotator_docker = "NO_ONCOTATOR",
                preemptible_attempts = preemptible_attempts
        }

		call Concordance {
		    input:
                intervals = intervals,
                truth_vcf = GoodTumor.filtered_vcf, #note, no orientation bias since it's optional output
                truth_vcf_idx = GoodTumor.filtered_vcf_index,
                eval_vcf = BadTumor.filtered_vcf,
                eval_vcf_idx = BadTumor.filtered_vcf_index,
                gatk_override = gatk_override
		}
	}

	output {
		Array[File] tpfn = Concordance.tpfn
        Array[File] tpfn_idx = Concordance.tpfn_idx
        Array[File] tpfp = Concordance.tpfp
        Array[File] tpf_idx = Concordance.tpfp_idx
        Array[File] summary = Concordance.summary
	}
}

task Concordance {
    File? intervals
    File truth_vcf
    File truth_vcf_idx
    File eval_vcf
    File eval_vcf_idx

    File? gatk_override

    String gatk_docker

    command {
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        gatk --java-options "-Xmx2g" Concordance \
            ${"-L " + intervals} \
            -truth ${truth_vcf} -eval ${eval_vcf} -tpfn "true_positives_and_false_negatives.vcf" \
            -tpfp "true_positives_and_false_positives.vcf" \
            -summary summary.tsv
    }

    runtime {
        docker: gatk_docker
        memory: "5 GB"
    }

    output {
        File tpfn = "true_positives_and_false_negatives.vcf"
        File tpfn_idx = "true_positives_and_false_negatives.vcf.idx"
        File tpfp = "true_positives_and_false_positives.vcf"
        File tpfp_idx = "true_positives_and_false_positives.vcf.idx"
        File summary = "summary.tsv"
    }
}