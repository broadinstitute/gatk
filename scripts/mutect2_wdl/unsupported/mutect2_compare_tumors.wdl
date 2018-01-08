# Given trios consisting of: a normal sample, a relatively reliable tumor sample, and a questionable tumor sample
# run Mutect in tumor-normal mode on the "good" sample, then run it on the "bad" sample, and compare results
# The idea is to see how well Mutect performs under difficult conditions.  For example, the "good" sample might
# be from a solid tumor while the "bad" sample might be low-allele-fraction cfDNA

# Exept for the different format of the input bams table, the inputs for this script are identical to those of
# mutect2_multi_sample.wdl
import "mutect2.wdl" as m2


task Concordance {
    String gatk
    File? gatk_override
    File? intervals
    File truth_vcf
    File truth_vcf_idx
    File eval_vcf
    File eval_vcf_idx

    command {
        GATK_JAR=${gatk}
        if [[ "${gatk_override}" == *.jar ]]; then
            GATK_JAR=${gatk_override}
        fi

        java -jar $GATK_JAR Concordance ${"-L " + intervals} \
            -truth ${truth_vcf} -eval ${eval_vcf} -tpfn "true_positives_and_false_negatives.vcf" \
            -tpfp "true_positives_and_false_positives.vcf" \
            -summary summary.tsv
    }

    runtime {
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

workflow Mutect2Trio {
	String gatk
	Int scatter_count
	# trio_list file is a tsv file with the following nine columns in this order.
	# normal_bam, normal_bai, normal_sample_name, good_tumor_bam, good_tumor_bai, good_tumor_sample_name, bad_tumor_bam, bad_tumor_bai, bad_tumor_sample_name,
	File trio_list
	Array[Array[String]] trios = read_tsv(trio_list)
	File? intervals
	File ref_fasta
	File ref_fai
	File ref_dict
	File? pon
	File? pon_index
	File? gnomad
	File? gnomad_index
	Boolean is_run_orientation_bias_filter
	Boolean is_run_oncotator
	String oncotator_docker
	String m2_docker
	File? gatk_override
	Int preemptible_attempts
	Array[String] artifact_modes

	scatter(trio in trios) {
		call m2.Mutect2 as GoodTumor {
			input:
				gatk=gatk,
				intervals=intervals,
				ref_fasta=ref_fasta,
				ref_fai=ref_fai,
				ref_dict=ref_dict,
				tumor_bam=trio[3],
				tumor_bai=trio[4],
				tumor_sample_name=trio[5],
				normal_bam=trio[0],
				normal_bai=trio[1],
				normal_sample_name=trio[2],
				pon=pon,
				pon_index=pon_index,
				scatter_count=scatter_count,
				gnomad=gnomad,
				gnomad_index=gnomad_index,
				picard = picard,
                is_run_orientation_bias_filter = is_run_orientation_bias_filter,
                is_run_oncotator=is_run_oncotator,
                oncotator_docker=oncotator_docker,
                m2_docker = m2_docker,
                gatk_override = gatk_override,
                preemptible_attempts = preemptible_attempts,
                artifact_modes = artifact_modes
		}

		call m2.Mutect2 as BadTumor {
            input:
        	    gatk=gatk,
        		intervals=intervals,
        		ref_fasta=ref_fasta,
       			ref_fai=ref_fai,
        		ref_dict=ref_dict,
        		tumor_bam=trio[6],
        		tumor_bai=trio[7],
        		tumor_sample_name=trio[8],
        		normal_bam=trio[0],
        		normal_bai=trio[1],
        		normal_sample_name=trio[2],
        		pon=pon,
        		pon_index=pon_index,
        		scatter_count=scatter_count,
        		gnomad=gnomad,
        		gnomad_index=gnomad_index,
        		picard = picard,
                is_run_orientation_bias_filter = is_run_orientation_bias_filter,
                is_run_oncotator=is_run_oncotator,
                oncotator_docker=oncotator_docker,
                m2_docker = m2_docker,
                gatk_override = gatk_override,
                preemptible_attempts = preemptible_attempts,
                artifact_modes = artifact_modes
        }

		call Concordance {
		    input:
		        gatk = gatk,
                gatk_override = gatk_override,
                intervals = intervals,
                truth_vcf = GoodTumor.filtered_vcf, #note, no orientation bias since it's optional output
                truth_vcf_idx = GoodTumor.filtered_vcf_index,
                eval_vcf = BadTumor.filtered_vcf,
                eval_vcf_idx = BadTumor.filtered_vcf_index,
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