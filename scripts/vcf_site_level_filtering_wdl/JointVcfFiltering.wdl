version 1.0

workflow JointVcfFiltering {
	input {
		Array[File] vcf
		Array[File] vcf_index
		File sites_only_vcf
		File sites_only_vcf_index
		String basename

		#TODO delete these inputs
		File python_script
		File hyperparameters_json

		String gatk_docker
		String? interval_contig

		Float indel_sensitivity_threshold
		Float snp_sensitivity_threshold
		String snp_annotations
		String indel_annotations
		File? gatk_override
	}

	call ExtractVariantAnnotations as ExtractVariantAnnotationsSNPs {
		input:
			input_vcf = sites_only_vcf,
			input_vcf_index = sites_only_vcf_index,
			mode = "SNP",
			annotations = snp_annotations,
			resources = "-resource:hapmap,training=true,calibration=true gs://gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz -resource:omni,training=true,calibration=true gs://gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz -resource:1000G,training=true,calibration=false gs://gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
			basename = basename,
			interval_contig = interval_contig,
			gatk_override = gatk_override,
			gatk_docker = gatk_docker
	}

	call ExtractVariantAnnotations as ExtractVariantAnnotationsINDELs {
		input:
			input_vcf = sites_only_vcf,
			input_vcf_index = sites_only_vcf_index,
			mode = "INDEL",
			annotations = indel_annotations,
			resources = "--resource:mills,training=true,calibration=true gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
			basename = basename,
			interval_contig = interval_contig,
			gatk_override = gatk_override,
			gatk_docker = gatk_docker
	}

	call TrainVariantAnnotationModel as TrainVariantAnnotationModelSNPs {
		input:
			annots = ExtractVariantAnnotationsSNPs.annots,
			basename = basename,
			mode_uc = "SNP",
			mode_lc = "snp",
			python_script = python_script,
			hyperparameters_json = hyperparameters_json,
			gatk_override = gatk_override,
			gatk_docker = gatk_docker
	}

	call TrainVariantAnnotationModel as TrainVariantAnnotationModelINDELs {
		input:
			annots = ExtractVariantAnnotationsINDELs.annots,
			basename = basename,
			mode_uc = "INDEL",
			mode_lc = "indel",
			python_script = python_script,
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
				scoring_python_script = python_script,
				annotations = snp_annotations,
				extracted_training_vcf = ExtractVariantAnnotationsSNPs.extracted_training_vcf,
				extracted_training_vcf_index = ExtractVariantAnnotationsSNPs.extracted_training_vcf_index,
				interval_contig = interval_contig,
				model = TrainVariantAnnotationModelSNPs.scorer,
				resources = "-resource:hapmap,training=false,calibration=true,prior=15 gs://gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz -resource:omni,training=false,calibration=true,prior=12 gs://gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz",
				gatk_override = gatk_override,
				gatk_docker = gatk_docker
		}

		call ScoreVariantAnnotations as ScoreVariantAnnotationsINDELs {
			input:
				vcf = ScoreVariantAnnotationsSNPs.output_vcf,
				vcf_index = ScoreVariantAnnotationsSNPs.output_vcf_index,
				basename = basename,
				mode = "INDEL",
				scoring_python_script = python_script,
				annotations = indel_annotations,
				extracted_training_vcf = ExtractVariantAnnotationsINDELs.extracted_training_vcf,
				extracted_training_vcf_index = ExtractVariantAnnotationsINDELs.extracted_training_vcf_index,
				model = TrainVariantAnnotationModelINDELs.scorer,
				resources = "--resource:mills,training=false,calibration=true gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
				gatk_override = gatk_override,
				gatk_docker = gatk_docker
		}
	}

	output {
		Array[File] variant_filtered_vcf = ScoreVariantAnnotationsINDELs.output_vcf
		Array[File] variant_filtered_vcf_index = ScoreVariantAnnotationsINDELs.output_vcf_index
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
		String resources
		String? interval_contig

		Int command_mem = 13000
	}
	Int disk_size = ceil(size(input_vcf, "GB") + 50)
	String interval_arg = if defined(interval_contig) then "-L " + interval_contig else ""
	command {
		set -e
		export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

		gatk --java-options "-Xmx${command_mem}m" \
		ExtractVariantAnnotations \
		~{interval_arg} \
		-V ~{input_vcf} \
		-O ~{basename}.~{mode} \
		~{annotations} \
		-mode ~{mode} \
		~{resources}
	}
	output {
		File annots = "~{basename}.~{mode}.annot.hdf5"
		File extracted_training_vcf = "~{basename}.~{mode}.vcf.gz"
		File extracted_training_vcf_index = "~{basename}.~{mode}.vcf.gz.tbi"
	}
	runtime {
		docker: gatk_docker
		disks: "local-disk " + disk_size + " LOCAL"
		memory: "14 GB"
	}
}

task TrainVariantAnnotationModel {
	input {
		String gatk_docker
		File? gatk_override
		File annots
		String basename
		String mode_uc
		String mode_lc
		File python_script
		File hyperparameters_json

		Int command_mem = 13000
	}
	Int disk_size = ceil(size(annots, "GB") + 100)
	command {
		set -e

		conda install -y --name gatk dill

		export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

		gatk --java-options "-Xmx${command_mem}m" \
		TrainVariantAnnotationsModel \
		--annotations-hdf5 ~{annots} \
		-O ~{basename} \
		--python-script ~{python_script} \
		--hyperparameters-json ~{hyperparameters_json} \
		-mode ~{mode_uc}

	}
	output {
		File scorer = "~{basename}.~{mode_lc}.scorer.pkl"
		File training_scores = "~{basename}.~{mode_lc}.trainingScores.hdf5"
		File truth_scores = "~{basename}.~{mode_lc}.calibrationScores.hdf5"
	}
	runtime {
		docker: gatk_docker
		disks: "local-disk " + disk_size + " LOCAL"
		memory: "14 GB"
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
		File scoring_python_script
		String annotations
		String resources
		File extracted_training_vcf
		File extracted_training_vcf_index
		String? interval_contig
		File model

		Int command_mem = 15000
	}
	Int disk_size = ceil(size(vcf, "GB") *2 + 50)
	String interval_arg = if defined(interval_contig) then "-L " + interval_contig else ""
	String model_basename = basename(model)
	command {
		set -e

		ln -s ~{model} ~{model_basename}

		conda install -y --name gatk dill

		export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

		gatk --java-options "-Xmx${command_mem}m" \
		ScoreVariantAnnotations \
		~{interval_arg} \
		-V ~{vcf} \
		-O ~{basename}.~{mode} \
		--python-script ~{scoring_python_script} \
		--model-prefix ~{basename} \
		~{annotations} \
		-mode ~{mode} \
		--resource:extracted-training,training=true,calibration=false ~{extracted_training_vcf} \
		~{resources}
	}
	output {
		File scores = "~{basename}.~{mode}.scores.hdf5"
		File annots = "~{basename}.~{mode}.annot.hdf5"
		File output_vcf = "~{basename}.~{mode}.vcf.gz"
		File output_vcf_index = "~{basename}.~{mode}.vcf.gz.tbi"
	}
	runtime {
		docker: gatk_docker
		disks: "local-disk " + disk_size + " LOCAL"
		memory: "16 GB"
	}
}

task ApplyVariantAnnotationScores {
	input {
		String gatk_docker
		File? gatk_override
		File input_vcf
		File input_vcf_index
		File tranches
		File recal_vcf
		File recal_vcf_idx
		String basename
		String filename_mode
		String mode
		String? interval_contig
		Float sensitivity

		Int command_mem = 15000
	}
	parameter_meta {
		input_vcf: {localization_optional: true}
	}
	Int disk_size = ceil(size(input_vcf, "GB") + 50)
	String interval_arg = if defined(interval_contig) then "-L " + interval_contig else ""
	command {
		set -e
		export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

		gatk --java-options "-Xmx${command_mem}m" \
		ApplyVariantAnnotationScores \
		~{interval_arg} \
		-V ~{input_vcf} \
		-O ~{basename}.~{filename_mode}.vcf.gz \
		--tranches-file ~{tranches} \
		--recal-file ~{recal_vcf} \
		-mode ~{mode} \
		--truth-sensitivity-filter-level ~{sensitivity}
	}
	output {
		File output_vcf = "~{basename}.~{filename_mode}.vcf.gz"
		File output_vcf_index = "~{basename}.~{filename_mode}.vcf.gz.tbi"
	}
	runtime {
		docker: gatk_docker
		disks: "local-disk " + disk_size + " LOCAL"
		memory: "16 GB"
	}
}

