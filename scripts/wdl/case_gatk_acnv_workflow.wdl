#
# Case sample workflow for a list of pairs of case-control samples. Includes GATK CNV and ACNV. Supports both WGS and WES samples.
#
# Notes:
#
# - the input file(input_bam_list) must contain a list of tab separated values in the following format(one or more lines must be supplied):
# tumor_entity  tumor_bam  tumor_bam_index  normal_entity  normal_bam  normal_bam_index  <--first input
# tumor_entity  tumor_bam  tumor_bam_index  normal_entity  normal_bam  normal_bam_index  <--second input
# etc...
#
# - set isWGS variable to true or false to specify whether to run a WGS or WES workflow respectively
#
# - file names will use the entity ID specified, but inside the file, the bam SM tag will typically be used.
#
# - target file (which must be in tsv format) is only used with WES workflow, WGS workflow generates its own targets (so user can pass any string as an argument)
#
# - the WGS PoN must be generated with WGS samples
# 
# - THIS SCRIPT SHOULD BE CONSIDERED OF "BETA" QUALITY
#
# - Example invocation:
#    java -jar cromwell.jar run case_gatk_acnv_workflow.wdl myParameters.json
# - See case_gatk_acnv_workflow_template.json for a template json file to modify with your own parameters (please save
#    your modified version with a different filename and do not commit to gatk-protected repo).
#
# - Some call inputs might seem out of place - consult with the comments in task definitions for details
#
#############

workflow case_gatk_acnv_workflow {
    # Workflow input files
    File target_file
    File ref_fasta
    File ref_fasta_dict
    File ref_fasta_fai
    File common_snp_list
    File input_bam_list
    Array[Array[String]] bam_list_array = read_tsv(input_bam_list)
    File PoN
    String gatk_jar

    # Input parameters of the PerformSegmentation tool
    Float seg_param_alpha
    Int seg_param_nperm
    String seg_param_pmethod
    Int seg_param_minWidth
    Int seg_param_kmax
    Int seg_param_nmin
    Float seg_param_eta
    Float seg_param_trim
    String seg_param_undoSplits
    Float seg_param_undoPrune
    Int seg_param_undoSD

    # CalculateTargetCoverage options
    Boolean disable_all_read_filters
    Boolean keep_duplicate_reads
    Boolean disable_sequence_dictionary_validation
    String transform
    String grouping

    # Input parameters of GetBayesianHetCoverage tool
    Float stringency
    Int read_depth_threshold

    # Workflow output directories and other options
    String plots_dir
    String conversion_dir
    Boolean enable_gc_correction
    Boolean isWGS
    Int wgsBinSize

    # Java maximum memory options
    Int calculate_target_coverage_memory
    Int normalize_somatic_read_count_memory
    Int whole_genome_coverage_memory

  call PadTargets {
    input:
        target_file=target_file,
        gatk_jar=gatk_jar,
        isWGS=isWGS,
        mem=1
  }

  scatter (row in bam_list_array) {

    call CalculateTargetCoverage as TumorCalculateTargetCoverage {
      input:
          entity_id=row[0],
          padded_target_file=PadTargets.padded_target_file,
          input_bam=row[1],
          bam_idx=row[2],
          ref_fasta=ref_fasta,
          ref_fasta_fai=ref_fasta_fai,
          ref_fasta_dict=ref_fasta_dict,
          gatk_jar=gatk_jar,
          disable_sequence_dictionary_validation=disable_sequence_dictionary_validation,
          disable_all_read_filters=disable_all_read_filters,
          keep_duplicate_reads=keep_duplicate_reads,
          transform=transform,
          grouping=grouping,
          isWGS=isWGS,
          mem=calculate_target_coverage_memory
    }  

    call WholeGenomeCoverage as TumorWholeGenomeCoverage {
      input:
          entity_id=row[0],
          target_file=PadTargets.padded_target_file,
          input_bam=row[1],
          bam_idx=row[2],
          coverage_file=TumorCalculateTargetCoverage.gatk_coverage_file,
          ref_fasta=ref_fasta,
          ref_fasta_fai=ref_fasta_fai,
          ref_fasta_dict=ref_fasta_dict,
          gatk_jar=gatk_jar,
          isWGS=isWGS,
          wgsBinSize=wgsBinSize,
          mem=whole_genome_coverage_memory
    }

    call AnnotateTargets as TumorAnnotateTargets {
      input:
          entity_id=row[0],
          gatk_jar=gatk_jar,
          target_file=TumorWholeGenomeCoverage.gatk_target_file,
          ref_fasta=ref_fasta,
          ref_fasta_fai=ref_fasta_fai,
          ref_fasta_dict=ref_fasta_dict,
          enable_gc_correction=enable_gc_correction,
          mem=4
    }

    call CorrectGCBias as TumorCorrectGCBias {
      input:
          entity_id=row[0],
          gatk_jar=gatk_jar,
          coverage_file=TumorWholeGenomeCoverage.gatk_coverage_file,
          annotated_targets=TumorAnnotateTargets.annotated_targets,
          enable_gc_correction=enable_gc_correction,
          mem=4
    }

    call NormalizeSomaticReadCounts as TumorNormalizeSomaticReadCounts {
      input:
          entity_id=row[0],
          coverage_file=TumorCorrectGCBias.gatk_cnv_coverage_file_gcbias,
          padded_target_file=TumorWholeGenomeCoverage.gatk_target_file,
          pon=PoN,
          gatk_jar=gatk_jar,
          mem=normalize_somatic_read_count_memory
    }

    call PerformSegmentation as TumorPerformSeg {
      input:
          entity_id=row[0],
          gatk_jar=gatk_jar,
          tn_file=TumorNormalizeSomaticReadCounts.tn_file,
          seg_param_alpha=seg_param_alpha,
          seg_param_nperm=seg_param_nperm,
          seg_param_pmethod=seg_param_pmethod,
          seg_param_minWidth=seg_param_minWidth,
          seg_param_kmax=seg_param_kmax,
          seg_param_nmin=seg_param_nmin,
          seg_param_eta=seg_param_eta,
          seg_param_trim=seg_param_trim,
          seg_param_undoSplits=seg_param_undoSplits,
          seg_param_undoPrune=seg_param_undoPrune,
          seg_param_undoSD=seg_param_undoSD,
          mem=2
    }

    call Caller as TumorCaller {
      input:
          entity_id=row[0],
          gatk_jar=gatk_jar,
          tn_file=TumorNormalizeSomaticReadCounts.tn_file,
          seg_file=TumorPerformSeg.seg_file,
          mem=2
    }

    call BayesianHetPulldownPaired {
      input:
          entity_id_tumor=row[0],
          entity_id_normal=row[3],
          gatk_jar=gatk_jar,
          ref_fasta=ref_fasta,
          ref_fasta_fai=ref_fasta_fai,
          ref_fasta_dict=ref_fasta_dict,
          tumor_bam=row[1],
          tumor_bam_idx=row[2],
          normal_bam=row[4],
          normal_bam_idx=row[5],
          common_snp_list=common_snp_list,
          stringency=stringency,
          read_depth_threshold=read_depth_threshold,
          mem=4
    }

    call AllelicCNV {
      input:
          entity_id=row[0],
          gatk_jar=gatk_jar,
          tumor_hets=BayesianHetPulldownPaired.tumor_hets,
          called_file=TumorCaller.called_file,
          tn_file=TumorNormalizeSomaticReadCounts.tn_file,
          mem=4
    }

    call PlotSegmentedCopyRatio {
      input:
          entity_id=row[0],
          gatk_jar=gatk_jar,
          tn_file=TumorNormalizeSomaticReadCounts.tn_file,
          pre_tn_file=TumorNormalizeSomaticReadCounts.pre_tn_file,
          called_file=TumorCaller.called_file,
          ref_fasta_dict=ref_fasta_dict,
          output_dir="${plots_dir}CopyRatio_Plots/${row[0]}/",
          mem=4
    }

    call PlotACNVResults {
      input:
          entity_id=row[0],
          gatk_jar=gatk_jar,
          tumor_hets=BayesianHetPulldownPaired.tumor_hets,
          acnv_segments=AllelicCNV.acnv_final_segs,
          tn_file=TumorNormalizeSomaticReadCounts.tn_file,
          ref_fasta_dict=ref_fasta_dict,
          output_dir="${plots_dir}ACNV_plots/${row[0]}",
          mem=4
    }

    call ConvertACNVResults {
      input:
          entity_id=row[0],
          gatk_jar=gatk_jar,
          tumor_hets=BayesianHetPulldownPaired.tumor_hets,
          acnv_segments=AllelicCNV.acnv_final_segs,
          tn_file=TumorNormalizeSomaticReadCounts.tn_file,
          output_dir="${conversion_dir}conversion/${row[0]}",
          mem=4
    }

    call CalculateTargetCoverage as NormalCalculateTargetCoverage {
      input:
          entity_id=row[3],
          padded_target_file=PadTargets.padded_target_file,
          input_bam=row[4],
          bam_idx=row[5],
          ref_fasta=ref_fasta,
          ref_fasta_fai=ref_fasta_fai,
          ref_fasta_dict=ref_fasta_dict,
          gatk_jar=gatk_jar,
          disable_sequence_dictionary_validation=disable_sequence_dictionary_validation,
          disable_all_read_filters=disable_all_read_filters,
          keep_duplicate_reads=keep_duplicate_reads,
          transform=transform,
          grouping=grouping,
          isWGS=isWGS,
          mem=calculate_target_coverage_memory
    }

    call WholeGenomeCoverage as NormalWholeGenomeCoverage {
      input:
          entity_id=row[3],
          target_file=PadTargets.padded_target_file,
          input_bam=row[4],
          bam_idx=row[5],
          coverage_file=NormalCalculateTargetCoverage.gatk_coverage_file,
          ref_fasta=ref_fasta,
          ref_fasta_fai=ref_fasta_fai,
          ref_fasta_dict=ref_fasta_dict,
          gatk_jar=gatk_jar,
          isWGS=isWGS,
          wgsBinSize=wgsBinSize,
          mem=whole_genome_coverage_memory
    }

    call AnnotateTargets as NormalAnnotateTargets{
      input:
          entity_id=row[3],
          gatk_jar=gatk_jar,
          target_file=NormalWholeGenomeCoverage.gatk_target_file,
          ref_fasta=ref_fasta,
          ref_fasta_fai=ref_fasta_fai,
          ref_fasta_dict=ref_fasta_dict,
          enable_gc_correction=enable_gc_correction,
          mem=4
    }

    call CorrectGCBias as NormalCorrectGCBias {
      input:
          entity_id=row[3],
          gatk_jar=gatk_jar,
          coverage_file=NormalWholeGenomeCoverage.gatk_coverage_file,
          annotated_targets=NormalAnnotateTargets.annotated_targets,
          enable_gc_correction=enable_gc_correction,
          mem=4
    }

    call NormalizeSomaticReadCounts as NormalNormalizeSomaticReadCounts {
      input:
          entity_id=row[3],
          coverage_file=NormalCorrectGCBias.gatk_cnv_coverage_file_gcbias,
          padded_target_file=NormalWholeGenomeCoverage.gatk_target_file,
          pon=PoN,
          gatk_jar=gatk_jar,
          mem=normalize_somatic_read_count_memory
    }

    call PerformSegmentation as NormalPerformSeg {
      input:
          entity_id=row[3],
          gatk_jar=gatk_jar,
          tn_file=NormalNormalizeSomaticReadCounts.tn_file,
          seg_param_alpha=seg_param_alpha,
          seg_param_nperm=seg_param_nperm,
          seg_param_pmethod=seg_param_pmethod,
          seg_param_minWidth=seg_param_minWidth,
          seg_param_kmax=seg_param_kmax,
          seg_param_nmin=seg_param_nmin,
          seg_param_eta=seg_param_eta,
          seg_param_trim=seg_param_trim,
          seg_param_undoSplits=seg_param_undoSplits,
          seg_param_undoPrune=seg_param_undoPrune,
          seg_param_undoSD=seg_param_undoSD,
          mem=2
    }

    call Caller as NormalCaller {
      input:
          entity_id=row[3],
          gatk_jar=gatk_jar,
          tn_file=NormalNormalizeSomaticReadCounts.tn_file,
          seg_file=NormalPerformSeg.seg_file,
          mem=2
    }
  }
}

# Pad the target file. This was found to help sensitivity and specificity. This step should only be altered
# by advanced users. Note that by changing this, you need to have a PoN that also reflects the change.
task PadTargets {
    File target_file
    Int padding
    String gatk_jar
    Boolean isWGS
    Int mem

    # Note that when isWGS is true, this task is still called by the workflow.
    # In that case, an empty target file is created and passed to the CalculateTargetCoverage 
    # task to satisfy input and output requirements.
    # Regardless of the value of isWGS, the output of PadTargets is passed to WholeGenomeCoverage, 
    # which then outputs the target file accordingly (see below)
    command {
        if [ ${isWGS} = false ]; \
          then java -Xmx${mem}g -jar ${gatk_jar} PadTargets --targets ${target_file} --output targets.padded.tsv \
                --padding ${padding} --help false --version false --verbosity INFO --QUIET false; \
          else touch targets.padded.tsv; \
        fi
    }

    output {
        File padded_target_file = "targets.padded.tsv"
    }

    #runtime {
    #    docker: "gatk-protected/a1"
    #}
}

# Calculate the target proportional coverage
task CalculateTargetCoverage {
    String entity_id
    File padded_target_file
    String transform
    String grouping
    Boolean keep_duplicate_reads
    Boolean disable_all_read_filters
    File input_bam
    File bam_idx
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    String gatk_jar
    Boolean disable_sequence_dictionary_validation
    Boolean isWGS
    Int mem

    # Note that when isWGS is true, this task is still called by the workflow.
    # In that case, an empty coverage file is created and passed to the WholeGenomeCoverage
    # task to satisfy input and output requirements.
    command <<<
        if [ ${isWGS} = false ]
          then
              java -Xmx${mem}g -jar ${gatk_jar} CalculateTargetCoverage --output ${entity_id}.coverage.tsv \
                --groupBy ${grouping} --transform ${transform} --targets ${padded_target_file} --targetInformationColumns FULL \
                --input ${input_bam} --reference ${ref_fasta} --disableToolDefaultReadFilters ${disable_all_read_filters} \
                $(if [ ${keep_duplicate_reads} = true ]; then echo " --disableReadFilter NotDuplicateReadFilter "; else echo ""; fi) \
                --interval_set_rule UNION --interval_padding 0 \
                --secondsBetweenProgressUpdates 10.0 --disableSequenceDictionaryValidation ${disable_sequence_dictionary_validation} \
                --createOutputBamIndex true --help false --version false --verbosity INFO --QUIET false
          else
              touch ${entity_id}.coverage.tsv
        fi
    >>>

    output {
        File gatk_coverage_file = "${entity_id}.coverage.tsv"
    }

    #runtime {
    #    docker: "gatk-protected/a1"
    #}
}

# Calculate coverage on Whole Genome Sequence using Spark.
# This task automatically creates a target output file.
task WholeGenomeCoverage {
    String entity_id
    File coverage_file
    File target_file
    File input_bam
    File bam_idx
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    String gatk_jar
    Boolean isWGS
    Int wgsBinSize
    Int mem

    # If isWGS is set to true, the task produces WGS coverage and targets that are passed to downstream tasks
    # If not, coverage and target files (received from upstream) for WES are passed downstream
    command {
        if [ ${isWGS} = true ]; \
          then java -Xmx${mem}g -jar ${gatk_jar} SparkGenomeReadCounts --outputFile ${entity_id}.coverage.tsv \
                --reference ${ref_fasta} --input ${input_bam} --sparkMaster local[1] --binsize ${wgsBinSize}; \
          else ln -s ${coverage_file} ${entity_id}.coverage.tsv; ln -s ${target_file} ${entity_id}.coverage.tsv.targets.tsv; \
        fi
    }

    output {
        File gatk_coverage_file = "${entity_id}.coverage.tsv"
        File gatk_target_file = "${entity_id}.coverage.tsv.targets.tsv"
    }
}

# Add new columns to an existing target table with various targets
# Note that this task is optional 
task AnnotateTargets {
    String entity_id
    File target_file
    String gatk_jar
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    Boolean enable_gc_correction
    Int mem

    # If GC correction is disabled, then an empty file gets passed downstream
    command {
        if [ ${enable_gc_correction} = true ]; \
          then java -Xmx${mem}g -jar ${gatk_jar} AnnotateTargets --targets ${target_file} --reference ${ref_fasta} --output ${entity_id}.annotated.tsv; \
          else touch ${entity_id}.annotated.tsv; \
        fi
    }

    output {
        File annotated_targets = "${entity_id}.annotated.tsv"
    }
}

# Correct coverage for sample-specific GC bias effects
# Note that this task is optional 
task CorrectGCBias {
    String entity_id
    File coverage_file
    File annotated_targets
    String gatk_jar
    Boolean enable_gc_correction
    Int mem

    # If GC correction is disabled, then the coverage file gets passed downstream unchanged
    command {
        if [ ${enable_gc_correction} = true ]; \
          then java -Xmx${mem}g -jar ${gatk_jar} CorrectGCBias --input ${coverage_file} \
           --output ${entity_id}.gc_corrected_coverage.tsv --targets ${annotated_targets}; \
          else ln -s ${coverage_file} ${entity_id}.gc_corrected_coverage.tsv; \
        fi
    }

    output {
        File gatk_cnv_coverage_file_gcbias = "${entity_id}.gc_corrected_coverage.tsv"
    }
}

# Perform tangent normalization (noise reduction) on the proportional coverage file.
task NormalizeSomaticReadCounts {
    String entity_id
    File coverage_file
    File padded_target_file
    File pon
    String gatk_jar
    Int mem

    command {
        java -Xmx${mem}g -jar ${gatk_jar} NormalizeSomaticReadCounts --input ${coverage_file} \
         --targets ${padded_target_file} --panelOfNormals ${pon} --factorNormalizedOutput ${entity_id}.fnt.tsv --tangentNormalized ${entity_id}.tn.tsv \
         --betaHatsOutput ${entity_id}.betaHats.tsv --preTangentNormalized ${entity_id}.preTN.tsv  --help false --version false --verbosity INFO --QUIET false
    }

    output {
        File tn_file = "${entity_id}.tn.tsv"
        File pre_tn_file = "${entity_id}.preTN.tsv"
        File betahats_file = "${entity_id}.betaHats.tsv"
    }
    #runtime {
    #    docker: "gatk-protected/a1"
    #}
}

# Segment the tangent normalized coverage profile.
task PerformSegmentation {
    String entity_id
    Float seg_param_alpha
    Int seg_param_nperm
    String seg_param_pmethod
    Int seg_param_minWidth
    Int seg_param_kmax
    Int seg_param_nmin
    Float seg_param_eta
    Float seg_param_trim
    String seg_param_undoSplits
    Float seg_param_undoPrune
    Int seg_param_undoSD
    String gatk_jar
    File tn_file
    Int mem

    command {
        java -Xmx${mem}g -jar ${gatk_jar} PerformSegmentation --tangentNormalized ${tn_file} \
         --output ${entity_id}.seg --log2Input true  --alpha ${seg_param_alpha} --nperm ${seg_param_nperm} \
         --pmethod ${seg_param_pmethod} --minWidth ${seg_param_minWidth} --kmax ${seg_param_kmax} \
         --nmin ${seg_param_nmin} --eta ${seg_param_eta} --trim ${seg_param_trim} --undoSplits ${seg_param_undoSplits} \
         --undoPrune ${seg_param_undoPrune} --undoSD ${seg_param_undoSD} --help false --version false \
         --verbosity INFO --QUIET false
    }

    output {
        File seg_file = "${entity_id}.seg"
    }
}

# Make calls (amp, neutral, or deleted) on each segment.
task Caller {
    String entity_id
    String gatk_jar
    File tn_file
    File seg_file
    Int mem

    command {
        java -Xmx${mem}g -jar ${gatk_jar} CallSegments --tangentNormalized ${tn_file} \
         --segments ${seg_file} --output ${entity_id}.called  --legacy false \
         --help false --version false --verbosity INFO --QUIET false
    }

    output {
        File called_file="${entity_id}.called"
    }
}

# Call heterozygous SNPs in the normal and then count the reads in the tumor for each called position.
# Entity IDs can be the same value
task HetPulldown {
    String entity_id_tumor
    String entity_id_normal
    String gatk_jar
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    File tumor_bam
    File tumor_bam_idx
    File normal_bam
    File normal_bam_idx
    File common_snp_list
    Int mem

    command {
        java -Xmx${mem}g -jar ${gatk_jar} GetHetCoverage --reference ${ref_fasta} \
         --normal ${normal_bam} --tumor ${tumor_bam} --snpIntervals ${common_snp_list} \
         --normalHets ${entity_id_normal}.normal.hets.tsv --tumorHets ${entity_id_tumor}.tumor.hets.tsv --pvalueThreshold 0.05 \
         --help false --version false --verbosity INFO --QUIET false --VALIDATION_STRINGENCY LENIENT
    }

    output {
        File normal_hets="${entity_id_normal}.normal.hets.tsv"
        File tumor_hets="${entity_id_tumor}.tumor.hets.tsv"
    }
}

# Call heterozygous SNPs in the normal and then count the reads in the tumor for each called position.
# Entity IDs can be the same value
task BayesianHetPulldownPaired {
    String entity_id_tumor
    String entity_id_normal
    String gatk_jar
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    File tumor_bam
    File tumor_bam_idx
    File normal_bam
    File normal_bam_idx
    File common_snp_list
    Float stringency
    Int read_depth_threshold
    Int mem

    command {
        java -Xmx${mem}g -jar ${gatk_jar} GetBayesianHetCoverage --reference ${ref_fasta} \
         --normal ${normal_bam} --tumor ${tumor_bam} --snpIntervals ${common_snp_list} \
         --normalHets ${entity_id_normal}.normal.hets.tsv --tumorHets ${entity_id_tumor}.tumor.hets.tsv \
         --hetCallingStringency ${stringency} --readDepthThreshold ${read_depth_threshold} \
         --help false --version false --verbosity INFO --QUIET false --VALIDATION_STRINGENCY LENIENT
    }

    output {
        File normal_hets="${entity_id_normal}.normal.hets.tsv"
        File tumor_hets="${entity_id_tumor}.tumor.hets.tsv"
    }
}

# Estimate minor allelic fraction and revise segments
task AllelicCNV {
    String entity_id
    String gatk_jar
    File tumor_hets
    File called_file
    File tn_file
    Int mem

    command {
        java -Xmx${mem}g -jar ${gatk_jar} AllelicCNV --tumorHets ${tumor_hets} \
         --tangentNormalized ${tn_file} --segments ${called_file} --outputPrefix ${entity_id} \
         --smallSegmentThreshold 3 --numSamplesCopyRatio 100 --numBurnInCopyRatio 50 --numSamplesAlleleFraction 100 \
         --numBurnInAlleleFraction 50 --intervalThresholdCopyRatio 5.0 --intervalThresholdAlleleFraction 2.0 \
         --help false --version false --verbosity INFO --QUIET false
    }

    output {
        File acnv_final_segs="${entity_id}-sim-final.seg"
    }
}

# Create plots of coverage data and copy-ratio estimates
task PlotSegmentedCopyRatio {
    String entity_id
    String gatk_jar
    File tn_file
    File pre_tn_file
    File called_file
    File ref_fasta_dict
    String output_dir
    Int mem

    command {
        mkdir -p ${output_dir} && \
        java -Xmx${mem}g -jar ${gatk_jar} PlotSegmentedCopyRatio --tangentNormalized ${tn_file} \
         --preTangentNormalized ${pre_tn_file} --segments ${called_file} \
         -SD ${ref_fasta_dict} \
         --output ${output_dir} --outputPrefix ${entity_id}
    }

    output {
        File segments_plot="${output_dir}/${entity_id}_FullGenome.png"
        File before_after_normalization_plot="${output_dir}/${entity_id}_Before_After.png"
        File before_after_cr_lim_4="${output_dir}/${entity_id}_Before_After_CR_Lim_4.png"
    }
}

# Create plots of allelic-count data and minor-allele-fraction estimates
task PlotACNVResults {
    String entity_id
    String gatk_jar
    File tumor_hets
    File tn_file
    File acnv_segments
    File ref_fasta_dict
    String output_dir
    Int mem

    command {
        mkdir -p ${output_dir} && \
        java -Xmx${mem}g -jar ${gatk_jar} PlotACNVResults --hets ${tumor_hets} \
         --tangentNormalized ${tn_file} --segments ${acnv_segments} \
         -SD ${ref_fasta_dict} \
         --output ${output_dir} --outputPrefix ${entity_id}
    }

    output { 
        File acnv_plot="${output_dir}/${entity_id}_ACNV.png"
    }
}

task ConvertACNVResults {
    String entity_id
    String gatk_jar
    File tumor_hets
    File acnv_segments
    File tn_file
    String output_dir
    Int mem

    command {
        mkdir -p ${output_dir} && \
        java -Xmx${mem}g -jar ${gatk_jar} ConvertACNVResults --tumorHets ${tumor_hets} \
         --segments ${acnv_segments} --tangentNormalized ${tn_file} \
         --outputDir ${output_dir}
    }

    output {
        File allelic_calls_final_acs_segs="${output_dir}/${entity_id}-sim-final.acs.seg"
        File allelic_calls_final_cnb_called_segs="${output_dir}/${entity_id}-sim-final.cnb_called.seg"
        File allelic_calls_final_cnv_segs="${output_dir}/${entity_id}-sim-final.cnv.seg"
        File allelic_calls_final_titan_hets="${output_dir}/${entity_id}-sim-final.titan.het.tsv"
        File allelic_calls_final_titan_tn="${output_dir}/${entity_id}-sim-final.titan.tn.tsv"
    }
}
