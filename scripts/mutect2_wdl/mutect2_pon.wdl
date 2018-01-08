#  Create a Mutect2 panel of normals
#
#  Description of inputs
#  gatk: java jar file containing gatk 4
#  picard: java jar file containing picard
#  intervals: genomic intervals
#  ref_fasta, ref_fai, ref_dict: reference genome, index, and dictionary
#  normal_bams, normal_bais: arrays of normal bams and bam indices
#  scatter_count: number of parallel jobs when scattering over intervals
#  pon_name: the resulting panel of normals is {pon_name}.vcf
#  duplicate_sample_strategy: THROW_ERROR (default if left empty) to fail if mulitple bams have the same sample name,
#                               CHOOSE_FIRST to use only the first bam with a given sample name, ALLOW_ALL to use all bams
#
# This WDL needs to decide whether to use the ``gatk_jar`` or ``gatk_jar_override`` for the jar location.  As of cromwell-0.24,
#   this logic *must* go into each task.  Therefore, there is a lot of duplicated code.  This allows users to specify a jar file
#   independent of what is in the docker file.  See the README.md for more info.
#
import "mutect2.wdl" as m2

workflow Mutect2_Panel {
    # gatk needs to be a String input to the workflow in order to work in a Docker image
	String gatk
	Int scatter_count
	Array[File] normal_bams
	Array[File] normal_bais
	File? intervals
	File ref_fasta
	File ref_fai
	File ref_dict
    String gatk_docker
    File? gatk_override
    Int? preemptible_attempts
    File picard
    String? m2_extra_args
    String? duplicate_sample_strategy
    String pon_name

    Array[Pair[File,File]] normal_bam_pairs = zip(normal_bams, normal_bais)

    scatter (normal_bam_pair in normal_bam_pairs) {
        File normal_bam = normal_bam_pair.left
        File normal_bai = normal_bam_pair.right

        call m2.Mutect2 {
            input:
                gatk = gatk,
                intervals = intervals,
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                ref_dict = ref_dict,
                tumor_bam = normal_bam,
                tumor_bai = normal_bai,
                scatter_count = scatter_count,
                gatk_docker = gatk_docker,
                gatk_override = gatk_override,
                preemptible_attempts = preemptible_attempts,
                m2_extra_args = m2_extra_args,
                is_run_orientation_bias_filter = false,
                is_run_oncotator = false,
                picard = picard,
                oncotator_docker = gatk_docker,
                artifact_modes = [""]
        }
    }

    call CreatePanel {
        input:
            gatk = gatk,
            input_vcfs = Mutect2.unfiltered_vcf,
            input_vcfs_idx = Mutect2.unfiltered_vcf_index,
            duplicate_sample_strategy = duplicate_sample_strategy,
            output_vcf_name = pon_name,
            gatk_override = gatk_override,
            preemptible_attempts = preemptible_attempts,
            gatk_docker = gatk_docker
    }

    output {
        File pon = CreatePanel.output_vcf
        File pon_idx = CreatePanel.output_vcf_index
        Array[File] normal_calls = Mutect2.unfiltered_vcf
        Array[File] normal_calls_idx = Mutect2.unfiltered_vcf_index
    }
}

task CreatePanel {
    String gatk
    Array[File] input_vcfs
    Array[File] input_vcfs_idx
    String? duplicate_sample_strategy
    String output_vcf_name
    File? gatk_override
    Int? preemptible_attempts
    String gatk_docker

    String gatk_local_jar = select_first([gatk_override, gatk])

    command {
        # Use GATK Jar override if specified
        export GATK_LOCAL_JAR=${gatk_local_jar}

        gatk --java-options "-Xmx2g"  CreateSomaticPanelOfNormals -vcfs ${sep=' -vcfs ' input_vcfs} ${"-duplicate-sample-strategy " + duplicate_sample_strategy} -O ${output_vcf_name}.vcf
    }

    runtime {
        docker: "${gatk_docker}"
        memory: "5 GB"
        disks: "local-disk " + 300 + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        File output_vcf = "${output_vcf_name}.vcf"
        File output_vcf_index = "${output_vcf_name}.vcf.idx"
    }
}