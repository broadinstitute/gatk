workflow CNVOncotateCalledSegments {
    File called_file
    Int? mem
    String? oncotator_docker
    String? additional_args

    call OncotateSegments {
        input:
            called_file = called_file,
            mem = mem,
            oncotator_docker = oncotator_docker,
            additional_args = additional_args
    }

    output {
        File oncotated_called_file = OncotateSegments.oncotated_called_file
        File oncotated_called_gene_list_file = OncotateSegments.oncotated_called_gene_list_file
    }
}

task OncotateSegments {
    File called_file
    String? additional_args
    String basename_called_file = basename(called_file)

    # Runtime parameters
    Int? mem
    String? oncotator_docker
    Int? preemptible_attempts
    Int? disk_space_gb

    command <<<
        set -e
        /root/oncotator_venv/bin/oncotator --db-dir /root/onco_dbdir/ -c /root/tx_exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt \
          -u file:///root/onco_cache/ -r -v ${called_file} ${basename_called_file}.per_segment.oncotated.txt hg19 \
          -i SEG_FILE -o SIMPLE_TSV ${default="" additional_args}

        /root/oncotator_venv/bin/oncotator --db-dir /root/onco_dbdir/ -c /root/tx_exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt \
          -u file:///root/onco_cache/ -r -v ${called_file} ${basename_called_file}.gene_list.txt hg19 \
          -i SEG_FILE -o GENE_LIST ${default="" additional_args}
    >>>

    runtime {
        docker: select_first([oncotator_docker, "broadinstitute/oncotator:1.9.5.0-eval-gatk-protected"])
        memory: select_first([mem, 3]) + " GB"
        disks: "local-disk " + select_first([disk_space_gb, 50]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        File oncotated_called_file = "${basename_called_file}.per_segment.oncotated.txt"
        File oncotated_called_gene_list_file = "${basename_called_file}.gene_list.txt"
    }
}