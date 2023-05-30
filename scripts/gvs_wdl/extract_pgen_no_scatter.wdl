version 1.0

workflow ExtractPgen {

    input {
        File input_vcf
        File input_vcf_index

        File? sample_list
        File? interval_list

        String gatk_docker
    }

    call GatkGvsPgenExtract {
        input:
            input_vcf = input_vcf,
            input_vcf_index = input_vcf_index,
            sample_list = sample_list,
            interval_list = interval_list,
            gatk_docker = gatk_docker
    }

    output {
        File pgen_file = GatkGvsPgenExtract.pgen_file
        File pvar_file = GatkGvsPgenExtract.pvar_file
        File psam_file = GatkGvsPgenExtract.psam_file
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
