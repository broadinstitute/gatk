version 1.0

workflow CNVOncotatorWorkflow {

    input {
      ##################################
      #### required basic arguments ####
      ##################################
      File called_file

      ##################################
      #### optional basic arguments ####
      ##################################
      String? additional_args
      String? oncotator_docker
      Int? mem_gb_for_oncotator
      Int? boot_disk_space_gb_for_oncotator
      Int? preemptible_attempts
    }

    call OncotateSegments {
        input:
            called_file = called_file,
            additional_args = additional_args,
            oncotator_docker = oncotator_docker,
            mem_gb = mem_gb_for_oncotator,
            boot_disk_space_gb = boot_disk_space_gb_for_oncotator,
            preemptible_attempts = preemptible_attempts
    }

    output {
        File oncotated_called_file = OncotateSegments.oncotated_called_file
        File oncotated_called_gene_list_file = OncotateSegments.oncotated_called_gene_list_file
    }
}

task OncotateSegments {

    input {
      File called_file
      String? additional_args

      # Runtime parameters
      String? oncotator_docker
      Int? mem_gb
      Int? disk_space_gb
      Int? boot_disk_space_gb
      Boolean use_ssd = false
      Int? cpu
      Int? preemptible_attempts
    }

    Int machine_mem_mb = select_first([mem_gb, 3]) * 1000

    String basename_called_file = basename(called_file)

    command <<<
        set -e

        # Get rid of the sequence dictionary at the top of the file
        egrep -v "^\@" ~{called_file} > ~{basename_called_file}.seq_dict_removed.seg

        echo "Starting the simple_tsv..."

        /root/oncotator_venv/bin/oncotator --db-dir /root/onco_dbdir/ -c /root/tx_exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt \
          -u file:///root/onco_cache/ -r -v ~{basename_called_file}.seq_dict_removed.seg ~{basename_called_file}.per_segment.oncotated.txt hg19 \
          -i SEG_FILE -o SIMPLE_TSV ~{default="" additional_args}

        echo "Starting the gene list..."

        /root/oncotator_venv/bin/oncotator --db-dir /root/onco_dbdir/ -c /root/tx_exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt \
          -u file:///root/onco_cache/ -r -v ~{basename_called_file}.seq_dict_removed.seg ~{basename_called_file}.gene_list.txt hg19 \
          -i SEG_FILE -o GENE_LIST ~{default="" additional_args}
    >>>

    runtime {
        docker: select_first([oncotator_docker, "broadinstitute/oncotator:1.9.5.0-eval-gatk-protected"])
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 50]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 2])
        bootDiskSizeGb: select_first([boot_disk_space_gb, 20])
    }

    output {
        File oncotated_called_file = "~{basename_called_file}.per_segment.oncotated.txt"
        File oncotated_called_gene_list_file = "~{basename_called_file}.gene_list.txt"
    }
}