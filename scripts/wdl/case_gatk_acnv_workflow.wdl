#
# Case sample workflow for a single pair of case-control samples.  Includes GATK CNV and ACNV.
#
# Notes:
#
# - the input file(input_bam_list) must contain a list of tab separated values in the following format(one or more lines must be supplied):
# tumor_entity  tumor_bam  tumor_bam_index  normal_entity  normal_bam  normal_bam_index  <--first input
# tumor_entity  tumor_bam  tumor_bam_index  normal_entity  normal_bam  normal_bam_index  <--second input
# etc...
#
# - file names will use the entity ID specified, but inside the file, the bam SM tag will typically be used.
#
# - THIS SCRIPT SHOULD BE CONSIDERED OF "BETA" QUALITY
#
# - For jni_lib values, see README.  For example, Ubuntu: /usr/lib/jni
# - Example invocation:
#    java -jar cromwell.jar run case_gatk_acnv_workflow.wdl myParameters.json
# - See case_gatk_acnv_workflow_template.json for a template json file to modify with your own parameters (please save
#    your modified version with a different filename and do not commit to gatk-protected repo).
#
###########

workflow case_gatk_acnv_workflow { 
    File target_bed
    File ref_fasta
    File ref_fasta_dict
    File ref_fasta_fai
    File common_snp_list
    File input_bam_list
    Array[Array[String]] bam_list_array = read_tsv(input_bam_list)
    File gatk_jar
    File PoN
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
    String is_disable_reference_validation
    String jni_lib    
    String plots_dir
    String call_cnloh_dir
    Boolean enable_gc_correction


  call PadTargets {
    input:
        target_bed=target_bed,
        gatk_jar=gatk_jar,
        mem=1
  }


  scatter (row in bam_list_array) {

    call AnnotateTargets {
      input:
          entity_id=row[0],
          gatk_jar=gatk_jar,
          padded_target_bed=PadTargets.padded_target_bed,
          ref_fasta=ref_fasta,
          ref_fasta_fai=ref_fasta_fai,
          ref_fasta_dict=ref_fasta_dict,
          enable_gc_correction=enable_gc_correction,
          mem=4
    }

    call CalculateTargetCoverage as TumorCalculateTargetCoverage {
      input:
          entity_id=row[0],
          padded_target_bed=PadTargets.padded_target_bed,
          input_bam=row[1],
          ref_fasta=ref_fasta,
          ref_fasta_fai=ref_fasta_fai,
          ref_fasta_dict=ref_fasta_dict,
          gatk_jar=gatk_jar,
          is_disable_reference_validation=is_disable_reference_validation,
          mem=2
    }  

    call CorrectGCBias as TumorCorrectGCBias {
      input:
          entity_id=row[0], 
          gatk_jar=gatk_jar,
          coverage_file=TumorCalculateTargetCoverage.gatk_cnv_coverage_file,
          annotated_targets=AnnotateTargets.annotated_targets,
          enable_gc_correction=enable_gc_correction,
          mem=4
    }

    call NormalizeSomaticReadCounts as TumorNormalizeSomaticReadCounts {
      input:
          entity_id=row[0], 
          coverage_file=TumorCorrectGCBias.gatk_cnv_coverage_file_gcbias,
          padded_target_bed=PadTargets.padded_target_bed,
          pon=PoN,
          gatk_jar=gatk_jar,
          mem=2,
          jni_lib=jni_lib
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

    call HetPulldown {
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
          mem = 4
    }

    call AllelicCNV {
      input:
          entity_id=row[0],
          gatk_jar=gatk_jar,
          tumor_hets=HetPulldown.tumor_hets,
          called_file=TumorCaller.called_file,
          tn_file=TumorNormalizeSomaticReadCounts.tn_file,
          mem=4
    }

    call PlotSegmentedCopyRatio {
      input:
          gatk_jar=gatk_jar,
          tn_file=TumorNormalizeSomaticReadCounts.tn_file,
          pre_tn_file=TumorNormalizeSomaticReadCounts.pre_tn_file,
          called_file=TumorCaller.called_file,
          output_dir="${plots_dir}CopyRatio_Plots/${row[0]}",
          mem=4
    }

    call PlotACNVResults {
      input:
          gatk_jar=gatk_jar,
          tumor_hets=HetPulldown.tumor_hets,
          acnv_segments=AllelicCNV.acnv_final_segs,
          tn_file=TumorNormalizeSomaticReadCounts.tn_file,
          output_dir="${plots_dir}ACNV_plots/${row[0]}",
          mem=4
    }

    call CNLoHAndSplitsCaller {
      input:
          entity_id=row[0],
          gatk_jar=gatk_jar,
          tumor_hets=HetPulldown.tumor_hets,
          acnv_segments=AllelicCNV.acnv_final_segs,
          tn_file=TumorNormalizeSomaticReadCounts.tn_file,
          output_dir="${call_cnloh_dir}CallCNLoH/${row[0]}",
          mem=4
    }

    call CalculateTargetCoverage as NormalCalculateTargetCoverage {
      input:
          entity_id=row[3],
          padded_target_bed=PadTargets.padded_target_bed,
          input_bam=row[4],
          ref_fasta=ref_fasta,
          ref_fasta_fai=ref_fasta_fai,
          ref_fasta_dict=ref_fasta_dict,
          gatk_jar=gatk_jar,
          is_disable_reference_validation=is_disable_reference_validation,
          mem=2
    }

    call CorrectGCBias as NormalCorrectGCBias {
      input:
          entity_id=row[3],
          gatk_jar=gatk_jar,
          coverage_file=NormalCalculateTargetCoverage.gatk_cnv_coverage_file,
          annotated_targets=AnnotateTargets.annotated_targets,
          enable_gc_correction=enable_gc_correction,
          mem=4
    }

    call NormalizeSomaticReadCounts as NormalNormalizeSomaticReadCounts {
      input:
          entity_id=row[3],
          coverage_file=NormalCorrectGCBias.gatk_cnv_coverage_file_gcbias,
          padded_target_bed=PadTargets.padded_target_bed,
          pon=PoN,
          gatk_jar=gatk_jar,
          mem=2,
          jni_lib=jni_lib
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

# Pad the target file.  This was found to help sensitivity and specificity.  This step should only be altered
#  by advanced users.  Note that by changing this, you need to have a PoN that also reflects the change.
task PadTargets {
    File target_bed
    Int padding = 250
    File gatk_jar
    Int mem

    command {
        java -Xmx${mem}g -jar ${gatk_jar} PadTargets  --targets ${target_bed} --output targets.padded.tsv \
         --padding ${padding}  --help false --version false --verbosity INFO --QUIET false
    }

    output {
        File padded_target_bed = "targets.padded.tsv"
    }
    #runtime {
    #    docker: "gatk-protected/a1"
    #}
}

# Calculate the target proportional coverage
task CalculateTargetCoverage {
    String entity_id
    File padded_target_bed
    String transform = "PCOV"
    String grouping = "SAMPLE"
    Boolean keep_duplicate_reads = true
    Boolean disable_all_read_filters = false
    File input_bam
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    File gatk_jar
    String is_disable_reference_validation
    Int mem

    command {
        java -Xmx${mem}g -jar ${gatk_jar} CalculateTargetCoverage --output ${entity_id}.coverage.tsv \
         --groupBy ${grouping} --transform ${transform} --targets ${padded_target_bed} --targetInformationColumns FULL \
         --keepduplicatereads ${keep_duplicate_reads} --input ${input_bam} --reference ${ref_fasta} \
         --disable_all_read_filters ${disable_all_read_filters} --interval_set_rule UNION --interval_padding 0 \
         --secondsBetweenProgressUpdates 10.0 --disableSequenceDictionaryValidation ${is_disable_reference_validation} \
         --createOutputBamIndex true --help false --version false --verbosity INFO --QUIET false
    }

    output {
        File gatk_cnv_coverage_file = "${entity_id}.coverage.tsv"
    }

    #runtime {
    #    docker: "gatk-protected/a1"
    #}
}

# Add new columns to an existing target table with various targets
task AnnotateTargets {
    String entity_id
    File padded_target_bed
    File gatk_jar
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    Boolean enable_gc_correction
    Int mem

    command {
        if [ ${enable_gc_correction} = true ]; \
          then java -Xmx${mem}g -jar ${gatk_jar} AnnotateTargets --targets ${padded_target_bed} --reference ${ref_fasta} --output ${entity_id}.annotated.tsv; \
          else touch ${entity_id}.annotated.tsv; \
        fi
    }

    output {
        File annotated_targets = "${entity_id}.annotated.tsv"
    }
}

# Correct coverage for sample-specific GC bias effects
task CorrectGCBias {
    String entity_id
    File coverage_file
    File annotated_targets
    File gatk_jar
    Boolean enable_gc_correction
    Int mem

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
    File padded_target_bed
    File pon
    File gatk_jar
    String jni_lib
    Int mem

    command {
        java -Xmx${mem}g -Djava.library.path=${jni_lib} -jar ${gatk_jar} NormalizeSomaticReadCounts  --input ${coverage_file} \
         --targets ${padded_target_bed} --panelOfNormals ${pon} --factorNormalizedOutput ${entity_id}.fnt.tsv --tangentNormalized ${entity_id}.tn.tsv \
         --betaHatsOutput ${entity_id}.betaHats.tsv --preTangentNormalized  ${entity_id}.preTN.tsv  --help false --version false --verbosity INFO --QUIET false
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
#  Many parameters are unexposed.
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
    File gatk_jar
    File tn_file
    Int mem

    command {
        java -Xmx${mem}g -jar ${gatk_jar} PerformSegmentation  --targets ${tn_file} \
         --output ${entity_id}.seg --log2Input true  --alpha ${seg_param_alpha} --nperm ${seg_param_nperm} \
         --pmethod ${seg_param_pmethod} --minWidth ${seg_param_minWidth} --kmax ${seg_param_kmax} \
         --nmin ${seg_param_nmin} --eta ${seg_param_eta} --trim ${seg_param_trim} --undoSplits ${seg_param_undoSplits} --undoPrune ${seg_param_undoPrune} --undoSD ${seg_param_undoSD} --help false --version false \
         --verbosity INFO --QUIET false
    }

    output {
        File seg_file = "${entity_id}.seg"
    }
}

# Make calls (amp, neutral, or deleted) on each segment.
task Caller {
    String entity_id
    File gatk_jar
    File tn_file
    File seg_file
    Int mem

    command {
        java -Xmx${mem}g -jar ${gatk_jar} CallSegments  --targets ${tn_file} \
         --segments ${seg_file} --output ${entity_id}.called  --legacy false \
         --help false --version false --verbosity INFO --QUIET false
    }

    output {
        File called_file="${entity_id}.called"
    }
}

# Call heterozygous SNPs in the normal and then count the reads in the tumor for each position called position.
# entity IDs can be the same value
task HetPulldown {
    String entity_id_tumor
    String entity_id_normal
    File gatk_jar
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
        java -Xmx${mem}g -jar ${gatk_jar} GetHetCoverage  --reference ${ref_fasta} \
         --normal ${normal_bam} --tumor ${tumor_bam} --snpIntervals ${common_snp_list}  \
         --normalHets ${entity_id_normal}.normal.hets.tsv --tumorHets ${entity_id_tumor}.tumor.hets.tsv --pvalueThreshold 0.05 \
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
    File gatk_jar
    File tumor_hets
    File called_file
    File tn_file
    Int mem

    command {
        java -Xmx${mem}g -jar ${gatk_jar} AllelicCNV  --tumorHets ${tumor_hets} \
         --tangentNormalized ${tn_file} --segments ${called_file} --outputPrefix ${entity_id} \
         --smallSegmentThreshold 3 --numSamplesCopyRatio 100 --numBurnInCopyRatio 50 --numSamplesAlleleFraction 100 \
         --numBurnInAlleleFraction 50 --intervalThresholdCopyRatio 5.0 --intervalThresholdAlleleFraction 2.0  \
         --help false --version false --verbosity INFO --QUIET false
    }

    output {
        File acnv_final_segs="${entity_id}-sim-final.seg"
    }
}


#Create plots of copy number variant data
task PlotSegmentedCopyRatio { 
    File gatk_jar
    File tn_file
    File pre_tn_file
    File called_file
    String output_dir
    Int mem

    command {
        mkdir -p ${output_dir} && \
        java -Xmx${mem}g -jar ${gatk_jar} PlotSegmentedCopyRatio --tangentNormalized ${tn_file} \
         --preTangentNormalized ${pre_tn_file} --segments ${called_file} \
         --output ${output_dir}
    }

    output {
        #corresponding output images
    }
}

#Create plots of allele fraction data
task PlotACNVResults {
    File gatk_jar
    File tumor_hets
    File tn_file
    File acnv_segments
    String output_dir
    Int mem

    command {
        mkdir -p ${output_dir} && \
        java -Xmx${mem}g -jar ${gatk_jar} PlotACNVResults --hets ${tumor_hets} \
         --tangentNormalized ${tn_file} --segments ${acnv_segments} \
         --output ${output_dir}
    }

    output {
        #corresponding output images
    }
}

#Calling ACNV segments as CNLoH and balanced
task CNLoHAndSplitsCaller {
    String entity_id
    File gatk_jar
    File tumor_hets
    File acnv_segments
    File tn_file
    String output_dir
    Int mem

    command {
        mkdir -p ${output_dir} && \
        java -Xmx${mem}g -jar ${gatk_jar} CallCNLoHAndSplits --tumorHets ${tumor_hets} \
         --segments ${acnv_segments} --tangentNormalized ${tn_file} \
         --outputDir ${output_dir}
    }

    output {
        File cnloh_final_acs_segs="${output_dir}/${entity_id}-sim-final.acs.seg" 
        File cnloh_final_cnb_called_segs="${output_dir}/${entity_id}-sim-final.cnb_called.seg"
        File cnloh_final_cnv_segs="${output_dir}/${entity_id}-sim-final.cnv.seg"
        File cnloh_final_titan_hets="${output_dir}/${entity_id}-sim-final.titan.het.tsv"
        File cnloh_final_titan_tn="${output_dir}/${entity_id}-sim-final.titan.tn.tsv"
    }
}
