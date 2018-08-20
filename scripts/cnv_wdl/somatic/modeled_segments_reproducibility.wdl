workflow CallModeledSegmentsReproducibilityValidation {
    File called_segs_1
    File called_segs_2
    String group_id
    String eval_docker
    File targets_file
    File? gatk4_jar_override
    File ref_fasta
    File ref_fasta_dict
    File ref_fasta_fai

    call ReproducibilityValidation {
        input: 
            called_segs_1 = called_segs_1,
            called_segs_2 = called_segs_2,
            group_id = group_id,
            eval_docker = eval_docker,
            targets_file = targets_file,
            gatk4_jar_override = gatk4_jar_override,
            ref_fasta = ref_fasta,
            ref_fasta_dict = ref_fasta_dict,
            ref_fasta_fai = ref_fasta_fai
    }
}


task ReproducibilityValidation {
    File called_segs_1
    File called_segs_2
    String group_id
    String eval_docker
    File targets_file
    File? gatk4_jar_override
    File ref_fasta
    File ref_fasta_dict
    File ref_fasta_fai

    # This should be a optional, but cromwell 30 croaks.
    Float ploidy = 3.7

    String base_targets_file = basename(targets_file)
    String sample1_name = basename(called_segs_1)
    String sample2_name = basename(called_segs_2)
    Boolean is_cr = false
    command <<<
        set -e
        # TODO: We need runtime parameters

        # Changing extension to work with CombineSegmentBreakpoints
        cp ${targets_file} targets_file.seg

        java -Xmx4g -jar ${default="/root/gatk.jar" gatk4_jar_override} CombineSegmentBreakpoints \
            --segments ${called_segs_1} --segments ${called_segs_2}  \
			--columns-of-interest CALL --columns-of-interest MEAN_LOG2_COPY_RATIO \
            -O reproducibility.tsv.seg -R ${ref_fasta}

        java -Xmx4g -jar ${default="/root/gatk.jar" gatk4_jar_override} CombineSegmentBreakpoints \
            --segments reproducibility.tsv.seg --segments targets_file.seg  \
			--columns-of-interest CALL_1 --columns-of-interest CALL_2 --columns-of-interest MEAN_LOG2_COPY_RATIO_1 --columns-of-interest MEAN_LOG2_COPY_RATIO_2 --columns-of-interest LOG2_COPY_RATIO \
            -O reproducibility_targets_tmp.tsv.seg -R ${ref_fasta}

        egrep -v "^\@" reproducibility_targets_tmp.tsv.seg > ${group_id}_reproducibility_targets.seg

        echo "Plotting...."
        python /root/run_plot_reproducibility.py \
            ${group_id}_reproducibility_targets.seg \
            ${sample1_name} \
            ${sample2_name} \
            ${group_id}/reproducibility/ \
            ${ploidy} \
            ${true='--cr' false='' is_cr}

        tar zcvf ${group_id}_reproducibility.tar.gz ${group_id}/reproducibility/
    >>>

    runtime {
        docker: "${eval_docker}"
        memory: "4 GB"
        disks: "local-disk 100 HDD"
        preemptible: 2
    }

    output {
        File final_reproducibility_validation_tar_gz = "${group_id}_reproducibility.tar.gz"
        File plot_reproducibility_input_file = "${group_id}_reproducibility_targets.seg"
    }
}
