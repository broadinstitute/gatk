version 1.0

import "https://raw.githubusercontent.com/gatk-workflows/gatk4-somatic-snvs-indels/2.6.0/mutect2.wdl" as m2

workflow Mutect3TrainingData {
    input {
        File? intervals
        File? masks
        File ref_fasta
        File ref_fai
        File ref_dict
        Int scatter_count
        File tumor_bam
        File tumor_bai
        File? normal_bam
        File? normal_bai
        File? pon
        File? gnomad
        File? variants_for_contamination
        String ref_downsample
        Boolean? run_orientation_bias_mixture_model_filter
        File? realignment_index_bundle
        String? realignment_extra_args
        String? m2_extra_args
        String? m2_extra_filtering_args
        String? normal_artifact_extra_args
        String? split_intervals_extra_args
        File? truth_vcf
        File? truth_vcf_idx
        Boolean? make_bamout

        # runtime
        String gatk_docker
        File? gatk_override
        Int? preemptible
        Int? max_retries
    }

    String m2_extra_args_with_training_mode = select_first([m2_extra_args, ""]) + " --training-data-mode --training-data-mode-ref-downsample " + ref_downsample

    Runtime small_runtime = {"gatk_docker": gatk_docker, "gatk_override": gatk_override,
                                   "max_retries": 2, "preemptible": 0, "cpu": 2,
                                   "machine_mem": 4000, "command_mem": 3500,
                                   "disk": 100, "boot_disk_size": 12}

    # call on the tumor (with normal if present) to get tumor read data and M2 filtering
    call m2.Mutect2 as Tumor {
        input:
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            scatter_count = scatter_count,
            tumor_reads = tumor_bam,
            tumor_reads_index = tumor_bai,
            normal_reads = normal_bam,
            normal_reads_index = normal_bai,
            intervals = intervals,
            pon = pon,
            gnomad = gnomad,
            variants_for_contamination = variants_for_contamination,
            run_orientation_bias_mixture_model_filter = run_orientation_bias_mixture_model_filter,
            realignment_index_bundle = realignment_index_bundle,
            realignment_extra_args = realignment_extra_args,
            preemptible = preemptible,
            max_retries = max_retries,
            m2_extra_args = m2_extra_args_with_training_mode,
            m2_extra_filtering_args = m2_extra_filtering_args,
            make_bamout = make_bamout,
            gatk_override = gatk_override,
            gatk_docker = gatk_docker
    }

    if(defined(truth_vcf)) {
        call Concordance  {
            input:
                intervals = intervals,
                masks = masks,
                truth_vcf = select_first([truth_vcf]),
                truth_vcf_idx = select_first([truth_vcf_idx]),
                eval_vcf = Tumor.filtered_vcf,
                eval_vcf_idx = Tumor.filtered_vcf_idx,
                preemptible = preemptible,
                gatk_override = gatk_override,
                gatk_docker = gatk_docker
        }

        call MakeTableFromConcordance as TumorConcordanceTable {
            input:
                tpfp = Concordance.tpfp,
                tpfp_idx = Concordance.tpfp_idx,
                ftnfn = Concordance.ftnfn,
                ftnfn_idx = Concordance.ftnfn_idx,
                gatk_override = gatk_override,
                gatk_docker = gatk_docker,
                preemptible = preemptible
        }
    }

    if(!defined(truth_vcf)) {
        call MakeTableFromMutect2 as TumorTable {
            input:
                filtered_vcf = Tumor.filtered_vcf,
                filtered_vcf_idx = Tumor.filtered_vcf_idx,
                gatk_override = gatk_override,
                gatk_docker = gatk_docker,
                preemptible = preemptible
        }
    }

    # call on the normal, with tumor as "matched normal", to get normal read data and M2 filtering
    if(defined(normal_bam)) {
        call m2.Mutect2 as Normal {
            input:
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                ref_dict = ref_dict,
                scatter_count = scatter_count,
                tumor_reads = select_first([normal_bam]),
                tumor_reads_index = select_first([normal_bai]),
                normal_reads = tumor_bam,
                normal_reads_index = tumor_bai,
                intervals = intervals,
                pon = pon,
                gnomad = gnomad,
                variants_for_contamination = variants_for_contamination,
                run_orientation_bias_mixture_model_filter = run_orientation_bias_mixture_model_filter,
                realignment_index_bundle = realignment_index_bundle,
                realignment_extra_args = realignment_extra_args,
                preemptible = preemptible,
                max_retries = max_retries,
                m2_extra_args = m2_extra_args_with_training_mode,
                m2_extra_filtering_args = m2_extra_filtering_args,
                make_bamout = make_bamout,
                gatk_override = gatk_override,
                gatk_docker = gatk_docker
        }

        # there's no reason to call concordance on the normal because the calls will have no relation to the truth VCF

        call MakeTableFromMutect2 as NormalTable {
            input:
                filtered_vcf = Normal.filtered_vcf,
                filtered_vcf_idx = Normal.filtered_vcf_idx,
                gatk_override = gatk_override,
                gatk_docker = gatk_docker,
                preemptible = preemptible
        }

        call m2.SplitIntervals as Split {
            input:
                intervals = intervals,
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                ref_dict = ref_dict,
                scatter_count = scatter_count,
                split_intervals_extra_args = split_intervals_extra_args,
                runtime_params = small_runtime
        }

        scatter (subintervals in Split.interval_files ) {
            call GetNormalArtifactData {
                input:
                    ref_fasta = ref_fasta,
                    ref_fai = ref_fai,
                    ref_dict = ref_dict,
                    tumor_reads = select_first([normal_bam]),
                    tumor_reads_index = select_first([normal_bai]),
                    normal_reads = tumor_bam,
                    normal_reads_index = tumor_bai,
                    intervals = subintervals,
                    preemptible = preemptible,
                    max_retries = max_retries,
                    extra_args = normal_artifact_extra_args,
                    gatk_override = gatk_override,
                    gatk_docker = gatk_docker
            }
        }

        call MergeNormalArtifactData {
            input:
                input_tables = GetNormalArtifactData.table,
                runtime_params = small_runtime
        }
    }

    output {
        File tumor_table = select_first([TumorConcordanceTable.table, TumorTable.table])
        File? normal_table = NormalTable.table
        File? normal_artifact_table = MergeNormalArtifactData.merged_table
    }
}

task Concordance {
    input {
        File? intervals
        File? masks
        File truth_vcf
        File truth_vcf_idx
        File eval_vcf
        File eval_vcf_idx

        File? gatk_override

        # runtime
        String gatk_docker
        Int? preemptible
    }

    command {
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        gatk --java-options "-Xmx2g" Concordance \
        ~{"-L " + intervals} \
        ~{"-XL " + masks} \
        -truth ~{truth_vcf} -eval ~{eval_vcf} \
        -tpfp "tpfp.vcf" \
        -ftnfn "ftnfn.vcf" \
        -summary "summary.txt"
    }

    runtime {
        memory: "5 GB"
        bootDiskSizeGb: 12
        docker: "${gatk_docker}"
        disks: "local-disk " + 100 + " HDD"
        preemptible: select_first([preemptible, 2])
    }

    output {
        File tpfp = "tpfp.vcf"
        File tpfp_idx = "tpfp.vcf.idx"
        File ftnfn = "ftnfn.vcf"
        File ftnfn_idx = "ftnfn.vcf.idx"
        File summary = "summary.txt"
    }
}

task MakeTableFromMutect2 {
    input {
        File filtered_vcf
        File filtered_vcf_idx

        File? gatk_override
        String gatk_docker
        Int? preemptible
    }

    command {
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        gatk --java-options "-Xmx2g" SelectVariants -V ~{filtered_vcf} --restrict-alleles-to BIALLELIC -O biallelic.vcf
        gatk --java-options "-Xmx2g" VariantsToTable -V biallelic.vcf \
          -F CHROM -F POS -F REF -F ALT -F POPAF -F TLOD -F STATUS -F REF_BASES -F HEC -F HAPDOM -F HAPCOMP -GF DP -F FILTER -GF FRS \
          --show-filtered \
          -O output.table
    }

    runtime {
        memory: "5 GB"
        bootDiskSizeGb: 12
        docker: "${gatk_docker}"
        disks: "local-disk " + 100 + " HDD"
        preemptible: select_first([preemptible, 2])
    }

    output {
        File table = "output.table"
    }
}

task MakeTableFromConcordance {
    input {
        File tpfp
        File tpfp_idx
        File ftnfn
        File ftnfn_idx

        File? gatk_override
        String gatk_docker
        Int? preemptible
    }

    command {
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        for file in ~{tpfp} ~{ftnfn}; do
            gatk --java-options "-Xmx2g" SelectVariants -V $file --restrict-alleles-to BIALLELIC -O biallelic.vcf
            gatk --java-options "-Xmx2g" VariantsToTable -V biallelic.vcf \
              -F CHROM -F POS -F REF -F ALT -F POPAF -F TLOD -F STATUS -F REF_BASES -F HEC -F HAPDOM -F HAPCOMP -GF DP -F FILTER -GF FRS \
              --show-filtered \
              -O tmp.table

            # if it's the first table, copy it to the output; otherwise copy all but the header line
            if [ ! -f output.table ]; then
                mv tmp.table output.table
            else
                tail -n +2 tmp.table >> output.table
            fi
        done
    }

    runtime {
        memory: "5 GB"
        bootDiskSizeGb: 12
        docker: "${gatk_docker}"
        disks: "local-disk " + 100 + " HDD"
        preemptible: select_first([preemptible, 2])
    }

    output {
        File table = "output.table"
    }
}

task GetNormalArtifactData {
    input {
        File? intervals
        File ref_fasta
        File ref_fai
        File ref_dict
        File tumor_reads
        File tumor_reads_index
        File? normal_reads
        File? normal_reads_index
        String? extra_args

        File? gatk_override
        String? gcs_project_for_requester_pays

        # runtime
        String gatk_docker
        Int? mem
        Int? preemptible
        Int? max_retries
        Int? disk_space
        Int? cpu
        Boolean use_ssd = false
    }

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1000 else 3500
    Int command_mem = machine_mem - 500

    parameter_meta{
        intervals: {localization_optional: true}
        ref_fasta: {localization_optional: true}
        ref_fai: {localization_optional: true}
        ref_dict: {localization_optional: true}
        tumor_reads: {localization_optional: true}
        tumor_reads_index: {localization_optional: true}
        normal_reads: {localization_optional: true}
        normal_reads_index: {localization_optional: true}
    }

    command <<<
        set -e

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

        if [[ ! -z "~{normal_reads}" ]]; then
            gatk --java-options "-Xmx~{command_mem}m" GetSampleName -R ~{ref_fasta} -I ~{normal_reads} -O normal_name.txt -encode \
            ~{"--gcs-project-for-requester-pays " + gcs_project_for_requester_pays}
            normal_sample="`cat normal_name.txt`"
        fi

        gatk --java-options "-Xmx~{command_mem}m" GetNormalArtifactData \
            -R ~{ref_fasta} ~{"-L " + intervals} -I ~{tumor_reads} -I ~{normal_reads} -O normal_artifact.table \
            -normal $normal_sample \
            ~{extra_args} ~{"--gcs-project-for-requester-pays " + gcs_project_for_requester_pays}
    >>>

    runtime {
        docker: gatk_docker
        bootDiskSizeGb: 12
        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, 100]) + if use_ssd then " SSD" else " HDD"
        preemptible: select_first([preemptible, 10])
        maxRetries: select_first([max_retries, 0])
        cpu: select_first([cpu, 1])
    }

    output {
        File table = "normal_artifact.table"
    }
}

task MergeNormalArtifactData {
    input {
        Array[File] input_tables
        Runtime runtime_params
    }

    command {
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}

        gatk --java-options "-Xmx~{runtime_params.command_mem}m" GatherNormalArtifactData \
            -I ~{sep=' -I ' input_tables} \
            -O normal_artifact.table
    }

    runtime {
        docker: runtime_params.gatk_docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.machine_mem + " MB"
        disks: "local-disk " + runtime_params.disk + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    output {
        File merged_table = "normal_artifact.table"
    }
}