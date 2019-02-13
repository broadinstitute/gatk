# Workflow to cluster a set of samples
#
#############

import "../cnv_common_tasks.wdl" as CNVTasks

workflow ClusterSamples {

    ##################################
    #### required basic arguments ####
    ##################################
    File intervals_for_clustering
    Array[String]+ normal_bams
    Array[String]+ normal_bais
    File ref_fasta_dict
    File ref_fasta_fai
    File ref_fasta
    String gatk_docker

    ##################################
    #### optional basic arguments ####
    ##################################
    File? gatk4_jar_override
    Int? preemptible_attempts

    ##############################################
    #### optional arguments for CollectCounts ####
    ##############################################
    String? collect_counts_format
    Int? mem_gb_for_collect_counts

    ###############################################
    #### optional arguments for ClusterSamples ####
    ###############################################
    File? clustering_prior_table
    Int? mem_gb_for_cluster_samples

    Array[Pair[String, String]] normal_bams_and_bais = zip(normal_bams, normal_bais)

    scatter (normal_bam_and_bai in normal_bams_and_bais) {
        call CNVTasks.CollectCounts {
            input:
                intervals = PreprocessIntervals.preprocessed_intervals,
                bam = normal_bam_and_bai.left,
                bam_idx = normal_bam_and_bai.right,
                ref_fasta = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,
                ref_fasta_dict = ref_fasta_dict,
                format = collect_counts_format,
                gatk4_jar_override = gatk4_jar_override,
                gatk_docker = gatk_docker,
                mem_gb = mem_gb_for_collect_counts,
                preemptible_attempts = preemptible_attempts
        }
    }

    call ClusterSamples {
        input:
            read_count_files = CollectCounts.counts,
            clustering_prior_table = clustering_prior_table,
            gatk_docker = gatk_docker,
            mem_gb = mem_gb_for_cluster_samples,
            preemptible_attempts = preemptible_attempts
    }

    output {
        clustering_table = ClusterSamples.clustering_table
    }
}

task ClusterSamples {
    Array[File] read_count_files
    File? clustering_prior_table

    # Runtime parameters
    String gatk_docker
    Int? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? cpu
    Int? preemptible_attempts

    Int machine_mem_mb = select_first([mem_gb, 3]) * 1000
    Int command_mem_mb = machine_mem_mb - 500

    command <<<
        #Jack's tool CLI invocation goes here
    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 150]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 8])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File clustering_table = "clustering_table.tsv"
    }
}
