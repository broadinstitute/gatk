# CRAM to trained CNNVariant Model

import "cnn_variant_common_tasks.wdl" as CNNTasks

workflow Cram2TrainedModel {
    File input_cram
    File reference_fasta
    File reference_dict
    File reference_fasta_index
    File truth_vcf
    File truth_vcf_index
    File truth_bed
    String output_prefix
    String tensor_type
    Int epochs
    File calling_intervals
    Int scatter_count
    String extra_args

    # Runtime parameters
    File? gatk_override
    String gatk_docker
    Int? mem_gb
    Int? preemptible_attempts
    Int? disk_space_gb
    Int? cpu

    Int? increase_disk_size
    Int additional_disk = select_first([increase_disk_size, 20])
    Float ref_size = size(reference_fasta, "GB") + size(reference_fasta_index, "GB") + size(reference_dict, "GB")

    call CNNTasks.CramToBam {
        input:
          reference_fasta = reference_fasta,
          reference_dict = reference_dict,
          reference_fasta_index = reference_fasta_index,
          cram_file = input_cram,
          output_prefix = output_prefix,
          disk_space_gb = disk_space_gb,
          preemptible_attempts = preemptible_attempts
    }

    call CNNTasks.SplitIntervals {
        input:
            scatter_count = scatter_count,
            intervals = calling_intervals,
            ref_fasta = reference_fasta,
            ref_dict = reference_dict,
            ref_fai = reference_fasta_index,
            gatk_docker = gatk_docker,
            gatk_override = gatk_override,
            preemptible_attempts = preemptible_attempts
    }

    Float bam_size = size(CramToBam.output_bam, "GB")

    scatter (calling_interval in SplitIntervals.interval_files) {
        call CNNTasks.RunHC4 {
            input:
                input_bam = CramToBam.output_bam,
                input_bam_index = CramToBam.output_bam_index,
                reference_fasta = reference_fasta,
                reference_dict = reference_dict,
                reference_fasta_index = reference_fasta_index,
                output_prefix = output_prefix,
                interval_list = calling_interval,
                gatk_docker = gatk_docker,
                gatk_override = gatk_override,
                preemptible_attempts = preemptible_attempts,
                extra_args = extra_args,
                disk_space_gb = round(bam_size + ref_size + additional_disk)                                 
        }

        call WriteTensors {
            input:
                input_vcf = RunHC4.raw_vcf,
                input_vcf_index = RunHC4.raw_vcf_index,
                input_bam = RunHC4.bamout,
                input_bam_index = RunHC4.bamout_index,
                truth_vcf = truth_vcf,
                truth_vcf_index = truth_vcf_index,
                truth_bed = truth_bed,
                tensor_type = tensor_type,
                reference_fasta = reference_fasta,
                reference_dict = reference_dict,
                reference_fasta_index = reference_fasta_index,
                output_prefix = output_prefix,
                interval_list = calling_interval,
                gatk_docker = gatk_docker,
                gatk_override = gatk_override,
                preemptible_attempts = preemptible_attempts,
                disk_space_gb = disk_space_gb
        }
    }

    call CNNTasks.MergeVCFs as MergeVCF_HC4 {
        input:
            input_vcfs = RunHC4.raw_vcf,
            output_prefix = output_prefix,
            gatk_override = gatk_override,
            gatk_docker = gatk_docker,
            preemptible_attempts = preemptible_attempts,
            disk_space_gb = disk_space_gb
    }

    call CNNTasks.SamtoolsMergeBAMs {
        input:
            input_bams = RunHC4.bamout,
            output_prefix = output_prefix,
            disk_space_gb = round(bam_size + ref_size + additional_disk)
    }

    call TrainModel {
        input:
            tar_tensors = WriteTensors.tensors,
            output_prefix = output_prefix,
            tensor_type = tensor_type,
            gatk_docker = gatk_docker,
            gatk_override = gatk_override,
            disk_space_gb = disk_space_gb,
            epochs = epochs
    }

    output {
        MergeVCF_HC4.*
        SamtoolsMergeBAMs.*
        TrainModel.*
    }

}

task WriteTensors {
    File input_bam
    File input_bam_index
    File input_vcf
    File input_vcf_index
    File reference_fasta
    File reference_dict
    File reference_fasta_index
    File truth_vcf
    File truth_vcf_index
    File truth_bed
    String output_prefix
    String tensor_type
    File interval_list

    # Runtime parameters
    String gatk_docker
    File? gatk_override
    Int? mem_gb
    Int? preemptible_attempts
    Int? disk_space_gb
    Int? cpu

    Int default_ram_mb = 8000

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem_gb) then mem_gb *1000 else default_ram_mb
    Int command_mem = machine_mem - 1000
    command {
        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        mkdir "/root/tensors/"

        gatk --java-options "-Xmx${command_mem}m" \
        CNNVariantWriteTensors \
        -R ${reference_fasta} \
        -V ${input_vcf} \
        -truth-vcf ${truth_vcf} \
        -truth-bed ${truth_bed} \
        -tensor-type ${tensor_type} \
        -output-tensor-dir "/root/tensors/" \
        -bam-file ${input_bam}
        
        tar -czf "tensors.tar.gz" "/root/tensors/"
    }

    output {
      File tensors = "tensors.tar.gz"
    }
    runtime {
        docker: "${gatk_docker}"
        memory: machine_mem + " MB"
        disks: "local-disk " + disk_space_gb + " SSD"
    }

}

task TrainModel {
    Array[File] tar_tensors
    String output_prefix
    String tensor_type
    Int epochs

    # Runtime parameters
    String gatk_docker
    File? gatk_override
    Int? mem_gb
    Int? preemptible_attempts
    Int? disk_space_gb
    Int? cpu

    Int default_ram_mb = 8000

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem_gb) then mem_gb *1000 else default_ram_mb
    Int command_mem = machine_mem - 1000
    command {
        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        for tensors in ${sep=' ' tar_tensors}  ; do
            tar -xzf $tensors 
        done

        gatk --java-options "-Xmx${command_mem}m" \
        CNNVariantTrain \
        -input-tensor-dir "./tensors/" \
        -model-name ${output_prefix} \
        -image-dir "./" \
        -tensor-type ${tensor_type} \
        -epochs ${epochs}
    }

    output {
        File model_json = "${output_prefix}.json"
        File model_hd5 = "${output_prefix}.hd5"
        File roc_png = "per_class_roc_${output_prefix}.png"
        File training_png = "metric_history_${output_prefix}.png"
    }

    runtime {
      docker: "${gatk_docker}"
      #gpuType: "nvidia-tesla-k80" # This will require PAPI v2 and CUDA on VM
      #gpuCount: 1
      #zones: ["us-central1-c"]
      memory: machine_mem + " MB"
      disks: "local-disk 400 SSD"
      bootDiskSizeGb: "16"
    }
}