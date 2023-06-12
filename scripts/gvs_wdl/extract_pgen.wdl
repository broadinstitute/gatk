version 1.0

workflow ExtractPgen {

    input {
        File input_vcf
        File input_vcf_index

        File? sample_list
        File? interval_list

        # Reference for splitting interval list
        File reference_fasta
        File reference_dict
        File reference_fasta_index

        String gatk_docker
        String plink_docker

        Int scatter_count
    }

    call SplitIntervals {
        input:
            intervals = interval_list,
            ref_fasta = reference_fasta,
            ref_fai = reference_fasta_index,
            ref_dict = reference_dict,
            scatter_count = scatter_count,
            gatk_docker = gatk_docker
    }

    scatter (calling_interval in SplitIntervals.interval_files) {
        call GatkGvsPgenExtract {
            input:
                input_vcf = input_vcf,
                input_vcf_index = input_vcf_index,
                sample_list = sample_list,
                interval_list = calling_interval,
                gatk_docker = gatk_docker
        }
    }

    call MergePgen {
        input:
            pgen_files = GatkGvsPgenExtract.pgen_file,
            pvar_files = GatkGvsPgenExtract.pvar_file,
            psam_files = GatkGvsPgenExtract.psam_file,
            plink_docker = plink_docker
    }

    output {
        File pgen_file = MergePgen.pgen_file
        File pvar_file = MergePgen.pvar_file
        File psam_file = MergePgen.psam_file
    }
}

task SplitIntervals {
    input {
        File? intervals
        File ref_fasta
        File ref_fai
        File ref_dict
        Int scatter_count

        # runtime
        String gatk_docker
        Int? mem
        Int? preemptible_attempts
        Int? disk_space
        Int? cpu
    }

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1000 else 3500
    Int command_mem = machine_mem - 500

    command {
        set -e
        export GATK_LOCAL_JAR="/root/gatk.jar"

        gatk --java-options "-Xmx${command_mem}m" \
        SplitIntervals \
        -R ${ref_fasta} \
        ${"-L " + intervals} \
        -scatter ${scatter_count} \
        -O ./
    }

    runtime {
        docker: "${gatk_docker}"
        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, 100]) + " HDD"
        preemptible: select_first([preemptible_attempts, 10])
        cpu: select_first([cpu, 1])
        bootDiskSizeGb: "16"
    }

    output {
        Array[File] interval_files = glob("*.interval_list")
    }
}

task GatkGvsPgenExtract {

    input {
        File input_vcf
        File input_vcf_index
        File? sample_list
        File? interval_list

        String gatk_docker
        Int? mem
        Int? preemptible_attempts
        Int? disk_space
        Int? cpu
    }

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1000 else 3500
    Int command_mem = machine_mem - 500

    command {
        set -e
        export GATK_LOCAL_JAR="/root/gatk.jar"

        gatk --java-options "-Xmx${command_mem}m" \
        GvsPgenExtractor \
        -V ${input_vcf} \
        ${"-sn " + sample_list} \
        ${"-L " + interval_list} \
        -O output.pgen
    }

    runtime {
        docker: "${gatk_docker}"
        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, 100]) + " HDD"
        preemptible: select_first([preemptible_attempts, 10])
        cpu: select_first([cpu, 1])
        bootDiskSizeGb: "16"
    }

    output {
        File pgen_file = "output.pgen"
        File pvar_file = "output.pvar"
        File psam_file = "output.psam"
    }
}

task MergePgen {
    input {
        Array[File] pgen_files
        Array[File] pvar_files
        Array[File] psam_files
        String plink_docker
    }

    command <<<
        set -e
        PGEN_ARRAY=(~{sep=" " pgen_files})
        touch mergelist.txt
        for pgen in PGEN_ARRAY
        do
            echo ${pgen%.pgen} >> mergelist.txt
        done

        ./plink2 --pmerge-list mergelist.txt --out merged
    >>>

    output {
        File pgen_file = "merged.pgen"
        File pvar_file = "merged.pvar"
        File psam_file = "merged.psam"
    }

    runtime {
        docker: "${plink_docker}"
    }
}
