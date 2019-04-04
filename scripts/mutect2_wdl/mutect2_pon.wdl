#  Create a Mutect2 panel of normals
#
#  Description of inputs
#  intervals: genomic intervals
#  ref_fasta, ref_fai, ref_dict: reference genome, index, and dictionary
#  normal_bams: arrays of normal bams
#  scatter_count: number of parallel jobs when scattering over intervals
#  pon_name: the resulting panel of normals is {pon_name}.vcf
#  m2_extra_args: additional command line parameters for Mutect2.  This should not involve --max-mnp-distance,
#  which the wdl hard-codes to 0 because GenpmicsDBImport can't handle MNPs

import "mutect2_nio.wdl" as m2
workflow Mutect2_Panel {
    # inputs
	File? intervals
	File ref_fasta
	File ref_fai
	File ref_dict
	Int scatter_count
	Array[String] normal_bams
	Array[String] normal_bais
	String gnomad
    String? m2_extra_args
    String? create_pon_extra_args
    Boolean? compress
    String pon_name

    Int? min_contig_size
    Int? num_contigs
    Int contig_size = select_first([min_contig_size, 1000000])

    File? gatk_override

    # runtime
    String gatk_docker
    Int? preemptible_attempts
    Int? max_retries

    scatter (normal_bam in zip(normal_bams, normal_bais)) {
        call m2.Mutect2 {
            input:
                intervals = intervals,
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                ref_dict = ref_dict,
                tumor_bam = normal_bam.left,
                tumor_bai = normal_bam.right,
                scatter_count = scatter_count,
                m2_extra_args = select_first([m2_extra_args, ""]) + "--max-mnp-distance 0",
                gatk_override = gatk_override,
                gatk_docker = gatk_docker,
                preemptible_attempts = preemptible_attempts,
                max_retries = max_retries
        }
    }

    call m2.SplitIntervals {
        input:
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            scatter_count = select_first([num_contigs, 24]),
            split_intervals_extra_args = "--subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION --min-contig-size " + contig_size,
            gatk_override = gatk_override,
            gatk_docker = gatk_docker
    }

    scatter (subintervals in SplitIntervals.interval_files ) {
            call CreatePanel {
                input:
                    input_vcfs = Mutect2.filtered_vcf,
                    intervals = subintervals,
                    ref_fasta = ref_fasta,
                    ref_fai = ref_fai,
                    ref_dict = ref_dict,
                    gnomad = gnomad,
                    output_vcf_name = pon_name,
                    create_pon_extra_args = create_pon_extra_args,
                    gatk_override = gatk_override,
                    preemptible_attempts = preemptible_attempts,
                    max_retries = max_retries,
                    gatk_docker = gatk_docker
            }
    }

    call m2.MergeVCFs {
        input:
            input_vcfs = CreatePanel.output_vcf,
            input_vcf_indices = CreatePanel.output_vcf_index,
            output_name = pon_name,
            compress = select_first([compress, false]),
            gatk_override = gatk_override,
            gatk_docker = gatk_docker,
            preemptible_attempts = preemptible_attempts,
            max_retries = max_retries
    }

    output {
        File pon = MergeVCFs.merged_vcf
        File pon_idx = MergeVCFs.merged_vcf_index
        Array[File] normal_calls = Mutect2.filtered_vcf
        Array[File] normal_calls_idx = Mutect2.filtered_vcf_index
    }
}

task CreatePanel {
    # inputs
    File intervals
    Array[String] input_vcfs
    File ref_fasta
    File ref_fai
    File ref_dict
    String output_vcf_name
    String gnomad
    String? create_pon_extra_args

    File? gatk_override

    # runtime
    String gatk_docker
    Int? mem
    Int? preemptible_attempts
    Int? max_retries
    Int? disk_space

    Int machine_mem = select_first([mem, 8])
    Int command_mem = machine_mem - 1

    command {
        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        gatk GenomicsDBImport --genomicsdb-workspace-path pon_db -R ${ref_fasta} -V ${sep=' -V ' input_vcfs} -L ${intervals}

        gatk --java-options "-Xmx${command_mem}g"  CreateSomaticPanelOfNormals -R ${ref_fasta} --germline-resource ${gnomad} \
            -V gendb://pon_db -O ${output_vcf_name}.vcf ${create_pon_extra_args}
    }

    runtime {
        docker: gatk_docker
        bootDiskSizeGb: 12
        memory: machine_mem + " GB"
        disks: "local-disk " + select_first([disk_space, 100]) + " HDD"
        preemptible: select_first([preemptible_attempts, 3])
        maxRetries: select_first([max_retries, 0])
    }

    output {
        File output_vcf = "${output_vcf_name}.vcf"
        File output_vcf_index = "${output_vcf_name}.vcf.idx"
    }
}