version 1.0

# Workflow for scoring and optionally filtering a VCF based on site-level annotations using the
# ExtractVariationAnnotations-TrainVariantAnnotationsModel-ScoreVariantAnnotations toolchain,
# which supersedes the corresponding VariantRecalibrator-ApplyVQSR toolchain.
# See the parameter_meta section below for descriptions of the workflow inputs.
# Also see the GATK documentation for these tools for descriptions of the corresponding methods and additional details.

struct RuntimeAttributes {
    Int? cpu
    Int? command_mem_gb
    Int? additional_mem_gb
    Int? disk_size_gb
    Int? boot_disk_size_gb
    Boolean? use_ssd
    Int? preemptible
    Int? max_retries
}

workflow JointVcfFiltering {
    input {
        Array[File] input_vcfs
        Array[File] input_vcf_idxs
        File sites_only_vcf
        File sites_only_vcf_idx
        String output_prefix

        Array[String] annotations
        String resource_args

        String? model_backend
        File? python_script
        File? hyperparameters_json

        String? extract_extra_args
        String? train_extra_args
        String? score_extra_args

        String gatk_docker
        File? gatk_override

        RuntimeAttributes? extract_runtime_attributes
        RuntimeAttributes? train_runtime_attributes
        RuntimeAttributes? score_runtime_attributes
    }

    parameter_meta {
        input_vcfs: "Sharded input VCFs to be scored and optionally filtered."
        sites_only_vcf: "A concatenated, sites-only version of the sharded input VCFs; used for extracting training and calibration sets."
        output_prefix: "Base prefix for output files. Sharded output VCFs will be named following the pattern \"{output_prefix}.{zero_based_shard_index}.score.vcf.gz\"."
        annotations: "Annotations to be used for extraction, training, and scoring."
        resource_args: "Resource arguments to be used for extraction and scoring. For example, \"--resource:training_and_calibration_set,training=true,calibration=true gs://path-to-training-and-calibration-set ...\".\n See GATK documentation for the ExtractVariantAnnotations and ScoreVariantAnnotations tools."
        model_backend: "(Optional) Model backend to be used by TrainVariantAnnotationsModel. See GATK documentation for this tool."
        python_script: "(Optional) Python script specifying custom model backend to be used by TrainVariantAnnotationsModel. See GATK documentation for this tool."
        hyperparameters_json: "(Optional) JSON file specifying model hyperparameters to be used by TrainVariantAnnotationsModel. See GATK documentation for this tool."
        extract_extra_args: "(Optional) Catch-all string to provide additional arguments for ExtractVariantAnnotations. This can include intervals (as string arguments or non-localized files), variant-type modes, arguments for enabling positive-negative training, etc. The \"do-not-gzip-vcf-output\" argument is not supported by this workflow. See GATK documentation for this tool."
        train_extra_args: "(Optional) Catch-all string to provide additional arguments for TrainVariantAnnotationsModel. This can include variant-type modes, arguments for enabling positive-negative training, etc. See GATK documentation for this tool."
        score_extra_args: "(Optional) Catch-all string to provide additional arguments for ScoreVariantAnnotations. This can include intervals (as string arguments or non-localized files), variant-type modes, arguments for enabling positive-negative training and hard filtering, etc. The \"do-not-gzip-vcf-output\" argument is not supported by this workflow. See GATK documentation for this tool."
    }

    call ExtractVariantAnnotations {
        input:
            input_vcf = sites_only_vcf,
            input_vcf_idx = sites_only_vcf_idx,
            output_prefix = output_prefix,
            annotations = annotations,
            resource_args = resource_args,
            extra_args = extract_extra_args,
            gatk_docker = gatk_docker,
            gatk_override = gatk_override,
            runtime_attributes = extract_runtime_attributes
    }

    call TrainVariantAnnotationsModel {
        input:
            annotations_hdf5 = ExtractVariantAnnotations.annotations_hdf5,
            unlabeled_annotations_hdf5 = ExtractVariantAnnotations.unlabeled_annotations_hdf5,
            model_backend = model_backend,
            python_script = python_script,
            hyperparameters_json = hyperparameters_json,
            output_prefix = output_prefix,
            extra_args = train_extra_args,
            gatk_docker = gatk_docker,
            gatk_override = gatk_override,
            runtime_attributes = train_runtime_attributes
    }

    scatter (i in range(length(input_vcfs))) {
        call ScoreVariantAnnotations {
            input:
                input_vcf = input_vcfs[i],
                input_vcf_idx = input_vcf_idxs[i],
                output_prefix = "~{output_prefix}.~{i}",
                annotations = annotations,
                resource_args = resource_args,
                extracted_vcf = ExtractVariantAnnotations.extracted_vcf,
                extracted_vcf_idx = ExtractVariantAnnotations.extracted_vcf_idx,
                model_prefix = output_prefix,
                model_files = TrainVariantAnnotationsModel.model_files,
                extra_args = score_extra_args,
                gatk_docker = gatk_docker,
                gatk_override = gatk_override,
                runtime_attributes = score_runtime_attributes
        }
    }

    output {
        File extracted_annotations_hdf5 = ExtractVariantAnnotations.annotations_hdf5
        File? extracted_unlabeled_annotations_hdf5 = ExtractVariantAnnotations.unlabeled_annotations_hdf5
        File extracted_vcf = ExtractVariantAnnotations.extracted_vcf
        File extracted_vcf_idx = ExtractVariantAnnotations.extracted_vcf_idx

        Array[File] model_files = TrainVariantAnnotationsModel.model_files

        Array[File] scored_vcfs = ScoreVariantAnnotations.scored_vcf
        Array[File] scored_vcf_idxs = ScoreVariantAnnotations.scored_vcf_idx
        Array[File?] annotations_hdf5s = ScoreVariantAnnotations.annotations_hdf5
        Array[File?] scores_hdf5s = ScoreVariantAnnotations.scores_hdf5
    }
}

task ExtractVariantAnnotations {
    input {
        File input_vcf
        File input_vcf_idx
        String output_prefix
        Array[String] annotations
        String resource_args
        String? extra_args

        String gatk_docker
        File? gatk_override

        RuntimeAttributes runtime_attributes = {}
    }

    parameter_meta {
        input_vcf: {localization_optional: true}
        input_vcf_idx: {localization_optional: true}
    }

    command {
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

        gatk --java-options "-Xmx~{default=6 runtime_attributes.command_mem_gb}G" \
            ExtractVariantAnnotations \
                -V ~{input_vcf} \
                -O ~{output_prefix}.extract \
                -A ~{sep=" -A " annotations} \
                ~{resource_args} \
                ~{extra_args}
    }

    runtime {
        docker: gatk_docker
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 100]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File annotations_hdf5 = "~{output_prefix}.extract.annot.hdf5"
        File? unlabeled_annotations_hdf5 = "~{output_prefix}.extract.unlabeled.annot.hdf5"
        File extracted_vcf = "~{output_prefix}.extract.vcf.gz"          # this line will break if extra_args includes the do-not-gzip-vcf-output argument
        File extracted_vcf_idx = "~{output_prefix}.extract.vcf.gz.tbi"  # this line will break if extra_args includes the do-not-gzip-vcf-output argument
    }
}

task TrainVariantAnnotationsModel {
    input {
        File annotations_hdf5
        File? unlabeled_annotations_hdf5
        String? model_backend
        File? python_script
        File? hyperparameters_json
        String output_prefix
        String? extra_args

        String gatk_docker
        File? gatk_override

        RuntimeAttributes runtime_attributes = {}
    }

    command {
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

        gatk --java-options "-Xmx~{default=6 runtime_attributes.command_mem_gb}G" \
            TrainVariantAnnotationsModel \
                --annotations-hdf5 ~{annotations_hdf5} \
                ~{"--unlabeled-annotations-hdf5 " + unlabeled_annotations_hdf5} \
                ~{"--model-backend " + model_backend} \
                ~{"--python-script " + python_script} \
                ~{"--hyperparameters-json " + hyperparameters_json} \
                -O ~{output_prefix}.train \
                ~{extra_args}
    }

    runtime {
        docker: gatk_docker
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 100]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        Array[File] model_files = glob("~{output_prefix}.train.*")
    }
}

task ScoreVariantAnnotations {
    input {
        File input_vcf
        File input_vcf_idx
        String output_prefix
        Array[String] annotations
        String resource_args
        File extracted_vcf
        File extracted_vcf_idx
        String model_prefix
        Array[File] model_files
        String? extra_args

        String gatk_docker
        File? gatk_override

        RuntimeAttributes runtime_attributes = {}
    }

    parameter_meta {
        input_vcf: {localization_optional: true}
        input_vcf_idx: {localization_optional: true}
        extracted_vcf: {localization_optional: true}
        extracted_vcf_idx: {localization_optional: true}
    }

    command {
        set -e

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

        mkdir model-files
        ln -s ~{sep=" model-files && ln -s " model_files} model-files

        gatk --java-options "-Xmx~{default=2 runtime_attributes.command_mem_gb}G" \
            ScoreVariantAnnotations \
                -V ~{input_vcf} \
                -O ~{output_prefix}.score \
                -A ~{sep=" -A " annotations} \
                ~{resource_args} \
                --resource:extracted,extracted=true ~{extracted_vcf} \
                --model-prefix model-files/~{model_prefix}.train \
                ~{extra_args}
    }

    runtime {
        docker: gatk_docker
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 2]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 100]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File scored_vcf = "~{output_prefix}.score.vcf.gz"               # this line will break if extra_args includes the do-not-gzip-vcf-output argument
        File scored_vcf_idx = "~{output_prefix}.score.vcf.gz.tbi"       # this line will break if extra_args includes the do-not-gzip-vcf-output argument
        File? annotations_hdf5 = "~{output_prefix}.score.annot.hdf5"    # this file will only be produced if the number of sites scored is nonzero
        File? scores_hdf5 = "~{output_prefix}.score.scores.hdf5"        # this file will only be produced if the number of sites scored is nonzero
    }
}
