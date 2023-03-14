version 1.0

# This is a workflow for filtering a joint callset VCF using INFO level annotations (so filtering is at the site level).
# Note that the input VCFs here may be sharded by genomic position which may be helpful for large cohorts. The script
# will output the same number of shards that are input.
# This portion of the filtering pipeline will assign a SCORE INFO field annotation to each site, but does not yet apply
# the filtering threshold to the final VCF.

workflow JointVcfFiltering {
	input {
		Array[File] vcf
		Array[File] vcf_index
		File sites_only_vcf
		File sites_only_vcf_index
		String basename

		String? model_backend
		File? training_python_script
		File? scoring_python_script
		File? hyperparameters_json

		String gatk_docker
		File? extract_interval_list
		File? score_interval_list

		String snp_annotations
		String indel_annotations
		File? gatk_override

		Boolean use_allele_specific_annotations

		String snp_resource_args = "--resource:hapmap,training=true,calibration=true gs://gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz --resource:omni,training=true,calibration=true gs://gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz --resource:1000G,training=true,calibration=false gs://gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
		String indel_resource_args = "--resource:mills,training=true,calibration=true gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --resource:axiom,training=true,calibration=true gs://gcp-public-data--broad-references/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz"
	}

	parameter_meta {
		vcf: "An array of input VCFs that are one callset sharded by genomic region."
		sites_only_vcf: "The full VCF callset without any genotype or sample level information."
		basename: "Desired output file basename."
	}

	call ExtractVariantAnnotations as ExtractVariantAnnotationsSNPs {
		input:
			input_vcf = sites_only_vcf,
			input_vcf_index = sites_only_vcf_index,
			mode = "SNP",
			annotations = snp_annotations,
			resource_args = snp_resource_args,
			basename = basename,
			interval_list = extract_interval_list,
			use_allele_specific_annotations = use_allele_specific_annotations,
			gatk_override = gatk_override,
			gatk_docker = gatk_docker
	}

	call ExtractVariantAnnotations as ExtractVariantAnnotationsINDELs {
		input:
			input_vcf = sites_only_vcf,
			input_vcf_index = sites_only_vcf_index,
			mode = "INDEL",
			annotations = indel_annotations,
			resource_args = indel_resource_args,
			basename = basename,
			interval_list = extract_interval_list,
			use_allele_specific_annotations = use_allele_specific_annotations,
			gatk_override = gatk_override,
			gatk_docker = gatk_docker
	}

	call TrainVariantAnnotationModel as TrainVariantAnnotationModelSNPs {
		input:
			annots = ExtractVariantAnnotationsSNPs.annots,
			basename = basename,
			mode = "snp",
			model_backend = model_backend,
			python_script = training_python_script,
			hyperparameters_json = hyperparameters_json,
			gatk_override = gatk_override,
			gatk_docker = gatk_docker
	}

	call TrainVariantAnnotationModel as TrainVariantAnnotationModelINDELs {
		input:
			annots = ExtractVariantAnnotationsINDELs.annots,
			basename = basename,
			mode = "indel",
			model_backend = model_backend,
			python_script = training_python_script,
			hyperparameters_json = hyperparameters_json,
			gatk_override = gatk_override,
			gatk_docker = gatk_docker
	}

	scatter(idx in range(length(vcf))) {
		call ScoreVariantAnnotations as ScoreVariantAnnotationsSNPs {
			input:
				vcf = vcf[idx],
				vcf_index = vcf_index[idx],
				basename = basename,
				mode = "SNP",
				model_backend = model_backend,
				python_script = scoring_python_script,
				annotations = snp_annotations,
				extracted_training_vcf = ExtractVariantAnnotationsSNPs.extracted_training_vcf,
				extracted_training_vcf_index = ExtractVariantAnnotationsSNPs.extracted_training_vcf_index,
				interval_list = score_interval_list,
				model_files = TrainVariantAnnotationModelSNPs.outputs,
				resource_args = snp_resource_args,
				use_allele_specific_annotations = use_allele_specific_annotations,
				gatk_override = gatk_override,
				gatk_docker = gatk_docker
		}

		call ScoreVariantAnnotations as ScoreVariantAnnotationsINDELs {
			input:
				vcf = vcf[idx],
				vcf_index = vcf_index[idx],
				basename = basename,
				mode = "INDEL",
				model_backend = model_backend,
				python_script = scoring_python_script,
				annotations = indel_annotations,
				extracted_training_vcf = ExtractVariantAnnotationsINDELs.extracted_training_vcf,
				extracted_training_vcf_index = ExtractVariantAnnotationsINDELs.extracted_training_vcf_index,
				interval_list = score_interval_list,
				model_files = TrainVariantAnnotationModelINDELs.outputs,
				resource_args = indel_resource_args,
				use_allele_specific_annotations = use_allele_specific_annotations,
				gatk_override = gatk_override,
				gatk_docker = gatk_docker
		}

	}

	output {
		Array[File] indels_variant_scored_vcf = ScoreVariantAnnotationsINDELs.output_vcf
		Array[File] indels_variant_scored_vcf_index = ScoreVariantAnnotationsINDELs.output_vcf_index
		Array[File] snps_variant_scored_vcf = ScoreVariantAnnotationsSNPs.output_vcf
		Array[File] snps_variant_scored_vcf_index = ScoreVariantAnnotationsSNPs.output_vcf_index
	}

}

task ExtractVariantAnnotations {
	input {
		String gatk_docker
		File? gatk_override
		File input_vcf
		File input_vcf_index
		String basename
		String mode
		String annotations
		String resource_args
		File? interval_list
		Boolean use_allele_specific_annotations

		Int memory_mb = 28000
		Int command_mem = memory_mb - 1000
	}
	Int disk_size = ceil(size(input_vcf, "GB") + size(input_vcf_index, "GB") + 100)

	File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"

	command {
		set -e

		bash ~{monitoring_script} > monitoring.log &

		export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

		gatk --java-options "-Xmx~{command_mem}m" \
			ExtractVariantAnnotations \
			-V ~{input_vcf} \
			-O ~{basename}.~{mode} \
			~{annotations} \
			~{if use_allele_specific_annotations then "--use-allele-specific-annotations" else ""} \
			~{"-L " + interval_list} \
			--mode ~{mode} \
			~{resource_args}
	}
	output {
		File annots = "~{basename}.~{mode}.annot.hdf5"
		File extracted_training_vcf = "~{basename}.~{mode}.vcf.gz"
		File extracted_training_vcf_index = "~{basename}.~{mode}.vcf.gz.tbi"
		Array[File] outputs = glob("~{basename}.~{mode}.*")
		File monitoring_log = "monitoring.log"
	}
	runtime {
		docker: gatk_docker
		disks: "local-disk " + disk_size + " HDD"
		memory: memory_mb + " MiB"
	}
}

task TrainVariantAnnotationModel {
	input {
		String gatk_docker
		File? gatk_override
		File annots
		String basename
		String mode
		String? model_backend
		File? python_script
		File? hyperparameters_json

		Int memory_mb = 28000
		Int command_mem = memory_mb - 1000
	}
	Int disk_size = ceil(size(annots, "GB") + 100)

	File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"

	command <<<
		set -e

		bash ~{monitoring_script} > monitoring.log &

		export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

		mode=$(echo "~{mode}" | awk '{print toupper($0)}')

		gatk --java-options "-Xmx~{command_mem}m" \
			TrainVariantAnnotationsModel \
			--annotations-hdf5 ~{annots} \
			-O ~{basename} \
			~{"--model-backend " + model_backend} \
			~{"--python-script " + python_script} \
			~{"--hyperparameters-json " + hyperparameters_json} \
			--mode $mode

	>>>
	output {
		Array[File] outputs = glob("~{basename}.~{mode}.*")
		File monitoring_log = "monitoring.log"
	}
	runtime {
		docker: gatk_docker
		disks: "local-disk " + disk_size + " HDD"
		memory: memory_mb + " MiB"
	}
}

task ScoreVariantAnnotations {
	input {
		String gatk_docker
		File? gatk_override
		File vcf
		File vcf_index
		String basename
		String mode
		String? model_backend
		File? python_script
		String annotations
		String resource_args
		File extracted_training_vcf
		File extracted_training_vcf_index
		File? interval_list
		Array[File] model_files
		Boolean use_allele_specific_annotations

		Int memory_mb = 16000
		Int command_mem = memory_mb - 1000
	}
	Int disk_size = ceil(size(vcf, "GB") * 2 + 50)

	File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"

	command {
		zgrep -v '#' ~{vcf} > empty.txt
		set -e

		bash ~{monitoring_script} > monitoring.log &

		if [ -s empty.txt ]; then
			ln -s ~{sep=" . && ln -s " model_files} .

			export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

			gatk --java-options "-Xmx~{command_mem}m" \
				ScoreVariantAnnotations \
				~{"-L " + interval_list} \
				-V ~{vcf} \
				-O ~{basename}.~{mode} \
				~{"--model-backend " + model_backend} \
				~{"--python-script " + python_script} \
				--model-prefix ~{basename} \
				~{annotations} \
				~{if use_allele_specific_annotations then "--use-allele-specific-annotations" else ""} \
				-mode ~{mode} \
				--resource:extracted,extracted=true ~{extracted_training_vcf} \
				~{resource_args}
		else
			echo "Input VCF was empty so we'll return the same VCF that was input."
			echo "Scores and annot hdf5 files will not be produced since the input was empty."
			ln -s ~{vcf} ~{basename}.~{mode}.vcf.gz
			ln -s ~{vcf_index} ~{basename}.~{mode}.vcf.gz.tbi
		fi
	}
	output {
		File? scores = "~{basename}.~{mode}.scores.hdf5"
		File? annots = "~{basename}.~{mode}.annot.hdf5"
		File output_vcf = "~{basename}.~{mode}.vcf.gz"
		File output_vcf_index = "~{basename}.~{mode}.vcf.gz.tbi"
		File monitoring_log = "monitoring.log"
	}
	runtime {
		docker: gatk_docker
		disks: "local-disk " + disk_size + " HDD"
		memory: memory_mb + " MiB"
	}
}

