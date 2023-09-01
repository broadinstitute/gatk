version 1.0

import "../cnv_common_tasks.wdl" as CNVTasks
import "https://raw.githubusercontent.com/broadinstitute/gatk-sv/v0.10-beta/wdl/AnnotateChromosome.wdl" as AnnotateVcf

workflow JointCallExomeCNVs {

    # NOTE: We don’t recommend using this for more than ~300 samples.

    ##################################
    #### required basic arguments ####
    ##################################
    input {
      Int num_samples_per_scatter_block
      File intervals
      File? blacklist_intervals

      File contig_ploidy_calls_tar_path_list
      File gcnv_calls_tars_path_list
      File genotyped_intervals_vcf_indexes_path_list
      File genotyped_intervals_vcfs_path_list
      File genotyped_segments_vcf_indexes_path_list
      File genotyped_segments_vcfs_path_list

      #qc arguments
      Int maximum_number_events
      Int maximum_number_pass_events

      Array[File] gcnv_model_tars
      Array[File] calling_configs
      Array[File] denoising_configs
      Array[File] gcnvkernel_version
      Array[File] sharded_interval_lists
      Array[String]? allosomal_contigs
      Int ref_copy_number_autosomal_contigs
      File ref_fasta_dict
      File ref_fasta_fai
      File ref_fasta
      String x_contig_name
      File protein_coding_gtf
      File linc_rna_gtf
      File promoter_bed
      File noncoding_bed
      String gatk_docker
      String gatk_docker_clustering
      String gatk_docker_qual_calc
      String sv_pipeline_docker
    }

    #we do these as FoFNs for Terra compatibility
    Array[File] contig_ploidy_calls_tars = read_lines(contig_ploidy_calls_tar_path_list)
    Array[File] segments_vcfs = read_lines(genotyped_segments_vcfs_path_list)
    Array[File] segments_vcf_indexes = read_lines(genotyped_segments_vcf_indexes_path_list)
    Array[File] intervals_vcfs = read_lines(genotyped_intervals_vcfs_path_list)
    Array[File] intervals_vcf_indexes = read_lines(genotyped_intervals_vcf_indexes_path_list)
    Array[Array[File]] gcnv_calls_tars = read_tsv(gcnv_calls_tars_path_list)
    Array[Array[File]] call_tars_sample_by_shard = transpose(gcnv_calls_tars)

    #create a ped file to use for allosome copy number (e.g. XX, XY)
    call MakePedFile {
      input:
        contig_ploidy_calls_tar = read_lines(contig_ploidy_calls_tar_path_list),
        x_contig_name = x_contig_name
    }

    call CNVTasks.SplitInputArray as SplitSegmentsVcfsList {
        input:
            input_array = segments_vcfs,
            num_inputs_in_scatter_block = num_samples_per_scatter_block,
            gatk_docker = gatk_docker
    }

    call CNVTasks.SplitInputArray as SplitSegmentsIndexesList {
        input:
            input_array = segments_vcf_indexes,
            num_inputs_in_scatter_block = num_samples_per_scatter_block,
            gatk_docker = gatk_docker
    }

    Array[Array[String]] split_segments = SplitSegmentsVcfsList.split_array
    Array[Array[String]] split_segments_indexes = SplitSegmentsIndexesList.split_array

    #for more than num_samples_per_scatter_block, do an intermediate combine first
    if (length(split_segments) > 1) {
      scatter (subarray_index in range(length(split_segments))) {
        call JointSegmentation as ScatterJointSegmentation {
          input:
            segments_vcfs = split_segments[subarray_index],
            segments_vcf_indexes = split_segments_indexes[subarray_index],
            ped_file = MakePedFile.ped_file,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            ref_fasta_dict = ref_fasta_dict,
            gatk_docker = gatk_docker_clustering,
            model_intervals = intervals
        }
      }
    }

    #refine breakpoints over all samples
    call JointSegmentation as GatherJointSegmentation {
      input:
        segments_vcfs = select_first([ScatterJointSegmentation.clustered_vcf, segments_vcfs]),
        segments_vcf_indexes = select_first([ScatterJointSegmentation.clustered_vcf_index, segments_vcfs]),
        ped_file = MakePedFile.ped_file,
        ref_fasta = ref_fasta,
        ref_fasta_fai = ref_fasta_fai,
        ref_fasta_dict = ref_fasta_dict,
        gatk_docker = gatk_docker_clustering,
        model_intervals = intervals
    }

    #recalculate each sample's quality scores based on new breakpoints and filter low QS or high AF events;
    #exclude samples with too many events
    scatter (scatter_index in range(length(segments_vcfs))) {
      call CNVTasks.PostprocessGermlineCNVCalls as RecalcQual {
        input:
              entity_id = sub(sub(basename(intervals_vcfs[scatter_index]), ".vcf.gz", ""), "intervals_output_", ""),
              gcnv_calls_tars = call_tars_sample_by_shard[scatter_index],
              gcnv_model_tars = gcnv_model_tars,
              calling_configs = calling_configs,
              denoising_configs = denoising_configs,
              gcnvkernel_version = gcnvkernel_version,
              sharded_interval_lists = sharded_interval_lists,
              contig_ploidy_calls_tar = read_lines(contig_ploidy_calls_tar_path_list)[0],  #this is always a list of one tar
              allosomal_contigs = allosomal_contigs,
              ref_copy_number_autosomal_contigs = ref_copy_number_autosomal_contigs,
              sample_index = scatter_index,
              maximum_number_events = maximum_number_events,
              maximum_number_pass_events = maximum_number_pass_events,
              intervals_vcf = intervals_vcfs[scatter_index],
              intervals_vcf_index = intervals_vcf_indexes[scatter_index],
              clustered_vcf = GatherJointSegmentation.clustered_vcf,
              clustered_vcf_index = GatherJointSegmentation.clustered_vcf_index,
              gatk_docker = gatk_docker_qual_calc
      }
    }

    #only put samples that passed QC into the combined VCF
    scatter(idx in range(length(RecalcQual.genotyped_segments_vcf))) {
      if (RecalcQual.qc_status_string[idx] == "PASS") {
        String subset = RecalcQual.genotyped_segments_vcf[idx]
        String subset_indexes = RecalcQual.genotyped_segments_vcf_index[idx]
      }
      if (RecalcQual.qc_status_string[idx] != "PASS") {
        String failed = sub(sub(basename(RecalcQual.genotyped_segments_vcf[idx]), ".vcf.gz", ""), "segments_output_", "")
      }
    }
    Array[String] subset_arr = select_all(subset)
    Array[String] subset_index_arr = select_all(subset_indexes)
    Array[String] failed_qc_samples = select_all(failed)

    call FastCombine {
      input:
        input_vcfs = subset_arr,
        input_vcf_indexes = subset_index_arr,
        sv_pipeline_docker = sv_pipeline_docker
    }

    # this annotates any vcf -- for exomes we can do all chromosomes at once
    call AnnotateVcf.AnnotateChromosome as Annotate {
        input:
          prefix = "combined.annotated",
          vcf = FastCombine.combined_vcf,
          protein_coding_gtf = protein_coding_gtf,
          linc_rna_gtf = linc_rna_gtf,
          promoter_bed = promoter_bed,
          noncoding_bed = noncoding_bed,
          sv_pipeline_docker = sv_pipeline_docker
    }

    output {
      File combined_calls = FastCombine.combined_vcf
      File combined_calls_index = FastCombine.combined_vcf_index
      File annotated_vcf = Annotate.annotated_vcf
      File annotated_vcf_index = Annotate.annotated_vcf_idx
      Array[String] sample_qc_status_strings = RecalcQual.qc_status_string
      Array[String] failed_samples = failed_qc_samples
    }
}

task MakePedFile {
  input {
    Array[File] contig_ploidy_calls_tar
    String x_contig_name
  }

  command <<<
    set -e

    while read tar
    do
      mkdir callsDir
      tar -xf $tar -C callsDir

      for sample in $(ls -d -1 callsDir/SAMPLE*)
      do
        sample_name=$(cat $sample/sample_name.txt)
        x_ploidy=$(grep ^~{x_contig_name} $sample/contig_ploidy.tsv | cut -f 2)
        [[ -z "$x_ploidy" ]] && { echo "Chromosome ~{x_contig_name} ploidy call not found for sample " $sample_name; exit 1; }
        printf "%s\t%s\t0\t0\t%s\t0\n" $sample_name $sample_name $x_ploidy >> cohort.ped
      done
      rm -rf callsDir
    done < ~{write_lines(contig_ploidy_calls_tar)}
    >>>

    output {
      File ped_file = "cohort.ped"
    }

    runtime {
      docker: "gatksv/sv-base-mini:b3af2e3"
      memory: "3000 MB"
      disks: "local-disk 100 SSD"
      cpu: 1
      preemptible: 2
    }
}

task JointSegmentation {
  input {
    Array[File] segments_vcfs
    Array[File] segments_vcf_indexes
    File ped_file
    File ref_fasta_dict
    File ref_fasta_fai
    File ref_fasta
    File model_intervals

     # Runtime parameters
    String gatk_docker
    Int? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? cpu
    Int? preemptible_attempts
    }

    parameter_meta {
      segments_vcfs: {localization_optional: true}
      segments_vcf_indexes: {localization_optional: true}
    }

    Int machine_mem_mb = select_first([mem_gb, 2]) * 1000
    Int command_mem_mb = machine_mem_mb - 500

  #NOTE: output has to be gzipped to be read in by pyvcf in the next step
  command <<<
    set -e
    gatk --java-options "-Xmx~{command_mem_mb}m" JointGermlineCNVSegmentation \
    -R ~{ref_fasta} -O clustered.vcf.gz -V ~{sep=' -V ' segments_vcfs} --model-call-intervals ~{model_intervals} -ped ~{ped_file}
    >>>

    output {
      File clustered_vcf = "clustered.vcf.gz"
      File clustered_vcf_index = "clustered.vcf.gz.tbi"
    }

    runtime {
      docker: gatk_docker
      memory: machine_mem_mb + " MB"
      disks: "local-disk " + select_first([disk_space_gb, 40]) + if use_ssd then " SSD" else " HDD"
      cpu: select_first([cpu, 1])
      preemptible: select_first([preemptible_attempts, 2])
    }
}

task FastCombine {
  input {
    Array[File] input_vcfs
    Array[File] input_vcf_indexes
    String sv_pipeline_docker
    Int? preemptible_tries
    Int? disk_size
    Float? mem_gb
  }

  command <<<
  #bcftools gets pissy if the indexes look older than their VCFs
  index_fofn=~{write_lines(input_vcf_indexes)}
  while read index; do touch $index; done < $index_fofn

  bcftools merge -l ~{write_lines(input_vcfs)} -o combined.vcf.gz -O z --threads 4 -m all -0

  tabix combined.vcf.gz
  >>>

  output {
    File combined_vcf = "combined.vcf.gz"
    File combined_vcf_index = "combined.vcf.gz.tbi"
  }

  runtime {
    docker: sv_pipeline_docker
    preemptible: select_first([preemptible_tries, 2])
    memory: select_first([mem_gb, 3.5]) + " GiB"
    cpu: "1"
    disks: "local-disk " + select_first([disk_size, 50]) + " HDD"
  }
}