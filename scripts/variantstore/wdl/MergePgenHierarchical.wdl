version 1.0

workflow MergePgenWorkflow {
    input {
        File pgen_file_list
        File pvar_file_list
        File psam_file_list
        String plink_docker
        String output_file_base_name
        Int? threads
        Int merge_disk_size
        Int split_count
        Boolean zero_padded_prefix
    }

    call SortFileLists {
        input:
            pgen_list = pgen_file_list,
            psam_list = psam_file_list,
            pvar_list = pvar_file_list,
            zero_padded_prefix = zero_padded_prefix
    }

    call SplitFileLists {
        input:
            pgen_list = SortFileLists.sorted_pgen_list,
            psam_list = SortFileLists.sorted_psam_list,
            pvar_list = SortFileLists.sorted_pvar_list,
            split_count = split_count
    }

    scatter(i in range(length(SplitFileLists.pgen_lists))) {
        call MergePgen as ScatterMerge {
            input:
                pgen_list = SplitFileLists.pgen_lists[i],
                psam_list = SplitFileLists.psam_lists[i],
                pvar_list = SplitFileLists.pvar_lists[i],
                plink_docker = plink_docker,
                output_file_base_name = "${output_file_base_name}${i}",
                threads = threads,
                disk_in_gb = merge_disk_size
           }
    }

    call MakeFileLists as MakeListsForFinal {
        input:
            pgen_files = ScatterMerge.pgen_file,
            psam_files = ScatterMerge.psam_file,
            pvar_files = ScatterMerge.pvar_file
    }

    call MergePgen as FinalMerge {
        input:
            pgen_list = MakeListsForFinal.pgen_list,
            psam_list = MakeListsForFinal.psam_list,
            pvar_list = MakeListsForFinal.pvar_list,
            plink_docker = plink_docker,
            output_file_base_name = output_file_base_name,
            threads = threads,
            disk_in_gb = merge_disk_size
    }

    output {
        File pgen_file = FinalMerge.pgen_file
        File pvar_file = FinalMerge.pvar_file
        File psam_file = FinalMerge.psam_file
    }
}

task MergePgen {
    input {
        File pgen_list
        File psam_list
        File pvar_list
        String plink_docker
        String output_file_base_name
        Int threads = 1
        Int disk_in_gb
    }

    Int cpu = threads + 1

    command <<<
        set -euxo pipefail

        # Download files using gsutil
        mkdir pgen_dir
        cat ~{pgen_list} | gsutil -m cp -I pgen_dir
        cat ~{psam_list} | gsutil -m cp -I pgen_dir
        cat ~{pvar_list} | gsutil -m cp -I pgen_dir

        # Create a file with a list of all the pgen basenames for merging
        touch mergelist.txt
        count=0
        for pgen in pgen_dir/*.pgen
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
            touch ~{output_file_base_name}.pgen
            touch ~{output_file_base_name}.pvar
            touch ~{output_file_base_name}.psam
            ;;
        1)
            echo "Only one pgen file so renaming the files for output"
            pgen_basename=$(cat mergelist.txt | xargs)
            mv ${pgen_basename}.pgen ~{output_file_base_name}.pgen
            mv ${pgen_basename}.psam ~{output_file_base_name}.psam
            if test -f ${pgen_basename}.pvar
            then
                mv ${pgen_basename}.pvar ~{output_file_base_name}.pvar
            else
                mv ${pgen_basename}.pvar.zst ~{output_file_base_name}.pvar.zst
            fi
            ;;
        *)
            echo "${count} pgen files, merging"
            if head -n 1 "~{pvar_list}" | grep -q "pvar$"
            then
                plink2 --pmerge-list mergelist.txt --memory 10000 --threads ~{threads} --out ~{output_file_base_name} --pmerge-output-vzs
            else
                plink2 --pmerge-list mergelist.txt pfile-vzs --memory 10000 --threads ~{threads} --out ~{output_file_base_name} --pmerge-output-vzs
            fi
            ;;
        esac

    >>>

    output {
        File pgen_file = "${output_file_base_name}.pgen"
        File pvar_file = "${output_file_base_name}.pvar.zst"
        File psam_file = "${output_file_base_name}.psam"
    }

    runtime {
        docker: "${plink_docker}"
        memory: "12 GB"
        disks: "local-disk ${disk_in_gb} HDD"
        bootDiskSizeGb: 15
        cpu: "${cpu}"
    }
}

task MakeFileLists {
    input {
        Array[File] pgen_files
        Array[File] pvar_files
        Array[File] psam_files
    }

    parameter_meta {
        pgen_files: {
            localization_optional: true
        }
        pvar_files: {
            localization_optional: true
        }
        psam_files: {
            localization_optional: true
        }
    }
    
    meta {
        # This causes issues when call cached for some reason, so we don't want to do that
        volatile: true
    }

    command <<<
        set -euxo pipefail

        touch pgen_list.txt
        touch psam_list.txt
        touch pvar_list.txt

        PGEN_ARRAY=(~{sep=" " pgen_files})
        PSAM_ARRAY=(~{sep=" " psam_files})
        PVAR_ARRAY=(~{sep=" " pvar_files})

        for i in "${!PGEN_ARRAY[@]}"
        do
            echo "${PGEN_ARRAY[$i]}" >> pgen_list.txt
            echo "${PSAM_ARRAY[$i]}" >> psam_list.txt
            echo "${PVAR_ARRAY[$i]}" >> pvar_list.txt
        done
    >>>

    output {
        File pgen_list = "pgen_list.txt"
        File psam_list = "psam_list.txt"
        File pvar_list = "pvar_list.txt"
    }

    runtime {
        docker: "ubuntu:22.04"
        memory: "1GB"
        bootDiskSizeGb: 15
    }
}

task SortFileLists {
    input {
        File pgen_list
        File psam_list
        File pvar_list

        Boolean zero_padded_prefix = true
    }

    command <<<
        set -euxo pipefail

        touch sorted_pgen_list.txt
        touch sorted_psam_list.txt
        touch sorted_pvar_list.txt

        readarray -t PGEN_ARRAY < ~{pgen_list}
        readarray -t PSAM_ARRAY < ~{psam_list}
        readarray -t PVAR_ARRAY < ~{pvar_list}

        # Maps of numbers to uris
        declare -A pgen_map
        declare -A psam_map
        declare -A pvar_map

        # Loop through arrays and add them to maps with the number prefix/suffix as the key
        for i in "${!PGEN_ARRAY[@]}"
        do
            # If zero_padded_prefix, that means we have to get the number from the start of the filename, before the first -
            # ex. gs://path/to/file/0001-basename.pgen
            if [ ~{zero_padded_prefix} = true ]
            then
                PGEN_NUM=$(echo "${PGEN_ARRAY[$i]}" | sed 's/.*\///' | sed 's/\-.*//')
                PSAM_NUM=$(echo "${PSAM_ARRAY[$i]}" | sed 's/.*\///' | sed 's/\-.*//')
                PVAR_NUM=$(echo "${PVAR_ARRAY[$i]}" | sed 's/.*\///' | sed 's/\-.*//')
            # If not zero_padded_prefix, that means we have to get the number from the end of the filename, after the last _
            # ex. gs://path/to/file/basename_1.pgen
            else
                PGEN_NUM=$(echo "${PGEN_ARRAY[$i]}" | sed -n 's/.*_\([^_]*\)\..*/\1/p')
                PSAM_NUM=$(echo "${PSAM_ARRAY[$i]}" | sed -n 's/.*_\([^_]*\)\..*/\1/p')
                PVAR_NUM=$(echo "${PVAR_ARRAY[$i]}" | sed -n 's/.*_\([^_]*\)\..*/\1/p')
            fi
            pgen_map[$PGEN_NUM]="${PGEN_ARRAY[$i]}"
            psam_map[$PSAM_NUM]="${PSAM_ARRAY[$i]}"
            pvar_map[$PVAR_NUM]="${PVAR_ARRAY[$i]}"
        done

        # Sort the keys numerically
        SORTED_PGEN_NUM=($(printf "%s\n" "${!pgen_map[@]}" | sort -n))
        SORTED_PSAM_NUM=($(printf "%s\n" "${!psam_map[@]}" | sort -n))
        SORTED_PVAR_NUM=($(printf "%s\n" "${!pvar_map[@]}" | sort -n))

        # Write to files
        for num in "${SORTED_PGEN_NUM[@]}"
        do
            echo "${pgen_map[$num]}" >> sorted_pgen_list.txt
        done
        for num in "${SORTED_PSAM_NUM[@]}"
        do
            echo "${psam_map[$num]}" >> sorted_psam_list.txt
        done
        for num in "${SORTED_PVAR_NUM[@]}"
        do
            echo "${pvar_map[$num]}" >> sorted_pvar_list.txt
        done
    >>>

    output {
        File sorted_pgen_list = "sorted_pgen_list.txt"
        File sorted_psam_list = "sorted_psam_list.txt"
        File sorted_pvar_list = "sorted_pvar_list.txt"
    }

    runtime {
        docker: "ubuntu:22.04"
        memory: "1GB"
        bootDiskSizeGb: 15
    }
}

task SplitFileLists {
    input {
        File pgen_list
        File psam_list
        File pvar_list

        Int split_count
    }

    command <<<
        set -euxo pipefail
        # Get the count of files and divide by split count (rounded up) to get number of files per split list
        FILE_COUNT=$(wc -l < ~{pgen_list})
        SPLIT_LINES=$(((FILE_COUNT+~{split_count}-1)/~{split_count}))
        # Split the lists
        split -l ${SPLIT_LINES} ~{pgen_list} split_pgen_files
        split -l ${SPLIT_LINES} ~{psam_list} split_psam_files
        split -l ${SPLIT_LINES} ~{pvar_list} split_pvar_files
    >>>

    output {
        Array[File] pgen_lists = glob("split_pgen_files*")
        Array[File] psam_lists = glob("split_psam_files*")
        Array[File] pvar_lists = glob("split_pvar_files*")
    }

    runtime {
        docker: "ubuntu:22.04"
        memory: "1GB"
        bootDiskSizeGb: 15
    }
}