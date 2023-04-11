version 1.0

import "ComparePgenAndVcf.wdl" as ComparePgenAndVcf

workflow ComparePgensAndVcfsScattered {
    input {
        File pgen_file_list
        File pvar_file_list
        File psam_file_list
        File vcfs_file_list

        Int split_count
        Int? report_disk_size
    }

    call SplitFileLists {
        input:
            pgen_list = pgen_file_list,
            psam_list = psam_file_list,
            pvar_list = pvar_file_list,
            vcf_list = vcfs_file_list,
            split_count = split_count
    }

    scatter(idx in range(length(SplitFileLists.pgen_lists))) {
        call ComparePgenAndVcf.ComparePgensAndVcfs {
            input:
                pgens = read_lines(SplitFileLists.pgen_lists[idx]),
                pvars = read_lines(SplitFileLists.pvar_lists[idx]),
                psams = read_lines(SplitFileLists.psam_lists[idx]),
                vcfs = read_lines(SplitFileLists.vcf_lists[idx]),
        }
    }

    call Report {
        input:
            diff_file_uris = flatten(ComparePgensAndVcfs.diffs),
            disk_in_gb = report_disk_size
    }

    output {
        Array[File] diffs = flatten(ComparePgensAndVcfs.diffs)
        File report = Report.report
        Int count = Report.count
    }
}

# Generates a report file based on the input diff files that lists the files with differences and the count of those
# files
task Report {
    input {
        Array[String] diff_file_uris

        Int disk_in_gb = 20
    }

    command <<<
        touch report.txt
        DIFF_ARRAY=(~{sep=" " diff_file_uris})
        count=0
        for i in "${!DIFF_ARRAY[@]}"
        do
            count=$((count+1))
            echo -e "${DIFF_ARRAY[$i]}" >> report.txt
        done
        echo -e "${count} files with differences" >> report.txt
        touch count.txt
        echo -e "${count}" > count.txt
    >>>

    output {
        File report = "report.txt"
        Int count = read_int("count.txt")
    }

    runtime {
        docker: "ubuntu:22.04"
        memory: "3 GB"
        disks: "local-disk ${disk_in_gb} HDD"
        preemptible: 3
        cpu: 1
    }
}

task SplitFileLists {
    input {
        File pgen_list
        File psam_list
        File pvar_list
        File vcf_list

        Int split_count
    }

    command <<<
        # Get the count of files and divide by split count (rounded up) to get number of files per split list
        FILE_COUNT=$(wc -l < ~{pgen_list})
        SPLIT_LINES=$(((FILE_COUNT+~{split_count}-1)/~{split_count}))
        # Split the lists
        split -l ${SPLIT_LINES} ~{pgen_list} split_pgen_files
        split -l ${SPLIT_LINES} ~{psam_list} split_psam_files
        split -l ${SPLIT_LINES} ~{pvar_list} split_pvar_files
        split -l ${SPLIT_LINES} ~{vcf_list} split_vcf_files
    >>>

    output {
        Array[File] pgen_lists = glob("split_pgen_files*")
        Array[File] psam_lists = glob("split_psam_files*")
        Array[File] pvar_lists = glob("split_pvar_files*")
        Array[File] vcf_lists = glob("split_vcf_files*")
    }

    runtime {
        docker: "ubuntu:22.04"
        memory: "1GB"
        bootDiskSizeGb: 15
    }
}