# The CNNScoreVariants tool annotates a VCF with scores from a Neural Net as part of a single-sample workflow.
# The site-level scores are added to the INFO field of the VCF.
# The architecture arguments, info_key and tensor type arguments MUST be in agreement
# (e.g. 2D models must have tensor_type of read_tensor and info_key CNN_2D, 1D models have tensor_type reference and info_key CNN_1D)
# The INFO field key will be "1D_CNN" or "2D_CNN" depending on the neural net architecture used for inference.
# The architecture arguments specify pre-trained networks.
# New networks can be trained by the GATK tools: CNNVariantWriteTensors and CNNVariantTrain
# The bam file and index are only required by 2D CNNs which take a read-level tensor_type such as "read_tensor".
# For 1D CNNs the tensor_type is typically "reference".
# Parallelization over sites is controlled by the scatter_count variable.

import "cnn_variant_common_tasks.wdl" as CNNTasks

workflow CNNScoreVariantsWorkflow {
    File input_vcf                  # The VCF to annotate with scores
    File input_vcf_index
    File reference_fasta
    File reference_dict
    File reference_fasta_index
    File resource_fofn              # File of VCF file names of resources of known SNPs and INDELs, (e.g. Mills, gnomAD)
    File resource_fofn_index        # File of VCF index file names for resources
    File? bam_file                  # Bam (or HaplotypeCaller-generated "bamout") file from which input_vcf was called, required by read-level architectures
    File? bam_file_index
    File? architecture_json         # Neural Net configuration for CNNScoreVariants
    File? architecture_hd5          # Pre-Trained weights and architecture for CNNScoreVariants
    String? tensor_type             # Keyword indicating the shape of the input tensor (e.g. read_tensor, reference)
    String info_key                 # The score key for the INFO field of the vcf (e.g. CNN_1D, CNN_2D)
    String snp_tranches             # Filtering threshold(s) for SNPs in terms of sensitivity to overlapping known variants in resources
    String indel_tranches           # Filtering threshold(s) for INDELs in terms of sensitivity to overlapping known variants in resources
    String output_prefix            # Identifying string for this run which will be used to name output files (the gzipped VCF and, for the 2D CNN, bamout)
    Int? inference_batch_size       # Batch size for python in CNNScoreVariants
    Int? transfer_batch_size        # Batch size for java transfers to python in CNNScoreVariants
    Int? intra_op_threads           # Tensorflow threading within nodes
    Int? inter_op_threads           # Tensorflow threading between nodes
    File? gatk_override
    String gatk_docker
    File calling_intervals
    Int scatter_count 
    Int? preemptible_attempts
    Int? cnn_task_mem_gb
    Int? cnn_task_cpu
    Int? mem_gb

    call CNNTasks.SplitIntervals {
        input:
            gatk_override = gatk_override,
            scatter_count = scatter_count,
            intervals = calling_intervals,
            ref_fasta = reference_fasta,
            ref_dict = reference_dict,
            ref_fai = reference_fasta_index,
            preemptible_attempts = preemptible_attempts,
            gatk_docker = gatk_docker
    }

    scatter (calling_interval in SplitIntervals.interval_files) {

        call CNNTasks.CNNScoreVariants {
            input:
                input_vcf = input_vcf,
                input_vcf_index = input_vcf_index,
                reference_fasta = reference_fasta,
                reference_dict = reference_dict,
                reference_fasta_index = reference_fasta_index,
                bam_file = bam_file,
                bam_file_index = bam_file_index,
                architecture_json = architecture_json,
                architecture_hd5 = architecture_hd5,
                tensor_type = tensor_type,
                inference_batch_size = inference_batch_size,
                transfer_batch_size = transfer_batch_size,
                intra_op_threads = intra_op_threads,
                inter_op_threads = inter_op_threads,
                output_prefix = output_prefix,
                interval_list = calling_interval,
                gatk_override = gatk_override,
                gatk_docker = gatk_docker,
                preemptible_attempts = preemptible_attempts,
                mem_gb = cnn_task_mem_gb,
                cpu = cnn_task_cpu
        }
    }

    call CNNTasks.MergeVCFs as MergeVCF_CNN {
        input: 
            input_vcfs = CNNScoreVariants.cnn_annotated_vcf,
            output_prefix = output_prefix,
            preemptible_attempts = preemptible_attempts,
            gatk_override = gatk_override,
            gatk_docker = gatk_docker
    }

    call CNNTasks.FilterVariantTranches {
        input:
            input_vcf = MergeVCF_CNN.merged_vcf,
            input_vcf_index = MergeVCF_CNN.merged_vcf_index,
            resource_fofn = resource_fofn,
            resource_fofn_index = resource_fofn_index,
            output_prefix = output_prefix,
            snp_tranches = snp_tranches,
            indel_tranches = indel_tranches,
            info_key = info_key,
            gatk_override = gatk_override,
            preemptible_attempts = preemptible_attempts,
            gatk_docker = gatk_docker
    }

    output {
        FilterVariantTranches.*
    }
}
