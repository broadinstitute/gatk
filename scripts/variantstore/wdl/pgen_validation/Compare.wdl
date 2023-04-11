version 1.0

workflow CompareFilesWorkflow {
    input{
        Array[File] actual
        Array[File] expected
    }

    call CompareFiles {
        input:
            actual = actual,
            expected = expected
    }

    output {
        Array[File] diffs = CompareFiles.diffs
    }
}

task CompareFiles {

    input{
        Array[File] actual
        Array[File] expected
    }

    Int file_count = length(actual)

    Int disk_in_gb = 2 * ceil(10 + size(actual, "GB") + size(expected, "GB"))

    command <<<
        set -euxo pipefail

        sort_by_basename() {
            local unsorted="$1"
            # Make a map of basenames to filenames
            declare -A basename_map
            while read -r filename
            do
                file_basename=$(basename "$filename")
                basename_map["$file_basename"]="$filename"
            done <"$unsorted"
            # Sort the basenames
            for basename in "${!basename_map[@]}"
            do 
                echo "$basename"; 
            done | sort > "sorted.txt"
            # Write the filenames sorted by basename to a file
            local sorted_filename="$2"
            while read -r file_basename
            do
                echo "${basename_map[$file_basename]}" >> "$sorted_filename"
            done <"sorted.txt"
        }

        generate_diff_file() {
            # Generate the diff file
            java -jar -Xmx1g /comparator/pgen_vcf_comparator.jar "$1" "$2" > "$(basename $1).diff"
            # If the diff file is empty, delete it
            if ! [ -s "$(basename $1).diff" ]
            then
                rm "$(basename $1).diff"
            fi
        }

        # Write the vcf lists to files and sort them by vcf basename so the vcfs match up correctly
        unsorted_actual=~{write_lines(actual)}
        sort_by_basename "$unsorted_actual" "sorted_actual.txt"
        unsorted_expected=~{write_lines(expected)}
        sort_by_basename "$unsorted_expected" "sorted_expected.txt"

        # Generate diff files for each pair of files
        for ((i = 1; i <= ~{file_count}; i++))
        do
            actual_file=$(sed -n "${i}p" "sorted_actual.txt")
            expected_file=$(sed -n "${i}p" "sorted_expected.txt")
            # Generate the diff file
            $(generate_diff_file "${actual_file}" "${expected_file}") &
        done

        wait
    >>>

    output {
        Array[File] diffs = glob("*.diff")
    }

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/klydon/pgen_vcf_comparator:test"
        memory: "12 GB"
        disks: "local-disk ${disk_in_gb} HDD"
        cpu: 10
    }
}