version 1.0

workflow ExtractPgen {

    input {
        Array[Array[File]] input_vcfs_with_indexes

        File? sample_list
        File? interval_list

        String gatk_docker
        String plink_docker
    }

    scatter (vcf_with_index in input_vcfs_with_indexes) {
        call GatkGvsPgenExtract {
            input:
                input_vcf = vcf_with_index[0],
                input_vcf_index = vcf_with_index[1],
                sample_list = sample_list,
                interval_list = interval_list,
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
        export GATK_LOCAL_JAR="/root/gatk.jar"

        gatk --java-options "-Xmx${command_mem}m" \
        GvsPgenExtractor \
        -V ${input_vcf} \
        ${"-sn " + sample_list} \
        ${"-L " + interval_list} \
        -O output.pgen

        if [[ ! -f output.pgen ]]
        then
            touch output.pgen
            touch output.psam
            touch output.pvar
        fi
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
        count=0
        for pgen in "${PGEN_ARRAY[@]}"
        do
            if [ -s ${pgen} ]
            then
                count=$((count+1))
                echo -e "${pgen%.pgen}" >> mergelist.txt
            fi
        done

        case $count in
            0)
                echo "No pgen files so creating empty ones"
                touch merged.pgen
                touch merged.pvar
                touch merged.psam
                ;;
            1)
                echo "Only one pgen file so renaming the files for output"
                pgen_basename=$(cat mergelist.txt | xargs)
                mv ${pgen_basename}.pgen merged.pgen
                mv ${pgen_basename}.psam merged.psam
                mv ${pgen_basename}.pvar merged.pvar
                ;;
            *)
                echo "${count} pgen files, merging"
                plink2 --pmerge-list mergelist.txt --out merged
                ;;
        esac

    >>>

    output {
        File pgen_file = "merged.pgen"
        File pvar_file = "merged.pvar"
        File psam_file = "merged.psam"
        File mergelist = "mergelist.txt"
    }

    runtime {
        docker: "${plink_docker}"
    }
}
