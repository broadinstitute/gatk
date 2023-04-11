version 1.0

import "Compare.wdl" as Compare

workflow ComparePgensAndVcfs {
    input {
        Array[File] pgens
        Array[File] pvars
        Array[File] psams
        Array[File] vcfs
    }

    call Convert {
        input:
            pgens = pgens,
            pvars = pvars,
            psams = psams,
            vcfs_for_size = vcfs
    }

    call Compare.CompareFiles {
        input:
            actual = Convert.output_vcfs,
            expected = vcfs
    }

    output {
        Array[File] diffs = CompareFiles.diffs
    }
}

# Uses plink to convert the supplied pgen and related files to vcfs
task Convert {
    input {
        Array[File] pgens
        Array[File] pvars
        Array[File] psams

        Array[File] vcfs_for_size
    }

    parameter_meta {
        vcfs_for_size: {
            localization_optional: true
        }
    }

    Int disk_in_gb = ceil(10 + 2 * size(vcfs_for_size, "GB"))

    command <<<
        set -e
        PGEN_ARRAY=(~{sep=" " pgens})
        # Get the extension of the first pvar file so we know if they're compressed
        PVAR_ARRAY=(~{sep=" " pvars})
        PVAR_EXT="${PVAR_ARRAY[0]##*.}"
        for i in "${!PGEN_ARRAY[@]}"
        do
            FILE_BASENAME="$(basename ${PGEN_ARRAY[$i]} .pgen)"
            FILENAME_NO_EXT="${PGEN_ARRAY[$i]%.pgen}"
            # If the pvar files are .pvar, use the normal mode, otherwise assume they are .pvar.zst and use vzs mode
            if [ "${PVAR_EXT}" == "pvar" ]
            then
                plink2 --pfile ${FILENAME_NO_EXT} --export vcf bgz --output-chr chrMT --out ${FILE_BASENAME}
            else
                plink2 --pfile ${FILENAME_NO_EXT} vzs --export vcf bgz --output-chr chrMT --out ${FILE_BASENAME}
            fi
        done
    >>>

    output {
        Array[File] output_vcfs = glob("*.vcf.gz")
    }

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/klydon/plink2:test"
        memory: "4 GB"
        disks: "local-disk ${disk_in_gb} HDD"
        preemptible: 3
        cpu: 1
    }
}
