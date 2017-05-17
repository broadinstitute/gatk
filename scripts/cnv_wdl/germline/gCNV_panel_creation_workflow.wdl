# Workflow for creating a panel of normals for germline CNV pipeline
# Notes:
#
# -Basic sex genotype tab-separated table for homo sapiens must be formatted as follows (Refer to the Javadoc of SexGenotypeTableReader for full description):
#    SAMPLE_NAME  SEX_GENOTYPE
#    sample_name_1 SEX_XX
#    sample_name_2 SEX_XY
#    sample_name_3  SEX_XY
#    sample_name_4  SEX_XX
# Sex genotype identifiers (SEX_XX and SEX_XY in the above example) must match those in the tab-separated germline contig ploidy annotation table. 
# The latter is formatted as follows:
#    CONTIG  CLASS  SEX_XX  SEX_XY
#    1         AUTOSOMAL       2          2
#    2         AUTOSOMAL       2          2
#    ...       ...             ...        ...
#    X         ALLOSOMAL       2          0
#    Y         ALLOSOMAL       1          1
#
# - Input file (normal_bams_list) must contain file paths to bam and bam index files separated by tabs in the following format:
#    normal_bam_1 bam_idx_1
#    normal_bam_2 bam_idx_2
#
# - The target file (targets) is required for the WES workflow and should be a tab-separated file with the column headers:
#    contig    start    stop    name
#   These targets will be padded on both sides by the amount specified by PadTargets.padding (default 250).
#
# - If a target file is not provided, then the WGS workflow will be run instead and the specified value of
#   wgs_bin_size (default 10000) will be used.
#
# - Example invocation:
#    java -jar cromwell.jar run gCNV_panel_creation_workflow.wdl myParameters.json
#   We recommend taking gCNV_cohort_calling_workflow.json as a template json file and modifying it accordingly (please save
#   your modified version with a different filename and do not commit to the gatk-protected repository).
##################

import "cnv_common_tasks.wdl" as CNVTasks

workflow CNVGermlinePanelWorkflow {
  # Workflow input files
  File? targets
  File normal_bams_list
  Array[Array[String]]+ normal_bams = read_tsv(normal_bams_list)
  File sex_genotypes 
  File contig_ploidy_annotations
  File transition_prior_table
  Array[File] copy_number_transition_prior_files
  File ref_fasta
  File ref_fasta_dict
  File ref_fasta_fai
  File gatk_jar

  # Model parameters
  Int num_latents
  # CombineReadCounts name
  String combined_entity_id = "combined_coverage"
  # Sex genotypes file name
  String sex_genotypes_entity_id = "sex_genotypes"
  # PoN output path
  String pon_output_path
  # If no target file is input, then do WGS workflow
  Boolean is_wgs = !defined(targets)

  if (!is_wgs) {
    call CNVTasks.PadTargets {
      input:
        targets = targets,
        gatk_jar = gatk_jar
    }
  }

  scatter (normal_bam in normal_bams) {
    call CNVTasks.CollectCoverage {
      input:
        padded_targets = PadTargets.padded_targets,
        keep_non_autosomes = true,
        bam = normal_bam[0],
        bam_idx = normal_bam[1],
        ref_fasta = ref_fasta,
        ref_fasta_fai = ref_fasta_fai,
        ref_fasta_dict = ref_fasta_dict,
        gatk_jar = gatk_jar,
        transform = "RAW"
    }
  }

  call CombineReadCounts {
    input:
      combined_entity_id = combined_entity_id,
      coverage_file_list = CollectCoverage.coverage,
      gatk_jar = gatk_jar
  }

  call CNVTasks.AnnotateTargets {
    input:
      entity_id = combined_entity_id,
      targets = CollectCoverage.coverage[0],
      ref_fasta = ref_fasta,
      ref_fasta_fai = ref_fasta_fai,
      ref_fasta_dict = ref_fasta_dict,
      gatk_jar = gatk_jar
  }

  call CNVTasks.CorrectGCBias {
    input:
      entity_id = combined_entity_id,
      coverage = CombineReadCounts.combined_coverage,
      annotated_targets = AnnotateTargets.annotated_targets,
      gatk_jar = gatk_jar
  } 

  call GermlineCNVCaller {
    input:
      coverage = CorrectGCBias.corrected_coverage,
      contig_ploidy_annotations = contig_ploidy_annotations,
      sex_genotypes = sex_genotypes,
      transition_prior_table = transition_prior_table,
      copy_number_transition_prior_files = copy_number_transition_prior_files,
      pon_output_path = pon_output_path,
      num_latents = num_latents, 
      gatk_jar = gatk_jar
  }

  output {
    Array[File] posteriors = GermlineCNVCaller.posteriors
    Array[File] model = GermlineCNVCaller.model
    Array[File] segments = GermlineCNVCaller.segments
  }
}

# Combine sample-level coverage files into a single file
task CombineReadCounts {
    String combined_entity_id
    Array[File]+ coverage_file_list
    Int? max_open_files
    File gatk_jar
    Int? mem

    command {
        java -Xmx${default=4 mem}g -jar ${gatk_jar} CombineReadCounts \
            --input ${sep=" --input " coverage_file_list} \
            --maxOpenFiles ${default=100 max_open_files} \
            --output ${combined_entity_id}.tsv
    }

    output {
        File combined_coverage = "${combined_entity_id}.tsv"
    }
}
  
# Learn the coverage model
task GermlineCNVCaller {
    File coverage
    File contig_ploidy_annotations
    File sex_genotypes
    File transition_prior_table
    Array[File] copy_number_transition_prior_files
    String pon_output_path
    Int num_latents
    File gatk_jar
    Int? mem

    command {
        java -Xmx${default=4 mem}g -Ddtype=double -jar ${gatk_jar} GermlineCNVCaller \
            --input ${coverage} \
            --contigAnnotationsTable ${contig_ploidy_annotations} \
            --sexGenotypeTable ${sex_genotypes} \
            --copyNumberTransitionPriorTable ${transition_prior_table} \
            --outputPath ${pon_output_path} \
            --jobType LEARN_AND_CALL \
            --numLatents ${default=5 num_latents} \
            --rddCheckpointing false \
            --disableSpark true
    }

    output {
        Array[File] posteriors = glob("./${pon_output_path}/posteriors_final/*")
        Array[File] model = glob("./${pon_output_path}/model_final/*") 
        Array[File] segments = glob("./${pon_output_path}/posteriors_final/segments/*")  
    }       
}
