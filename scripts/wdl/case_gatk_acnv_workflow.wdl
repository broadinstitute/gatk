#
# Case sample workflow for a single pair of case-control samples.  Includes GATK CNV and ACNV.
#
# Notes:
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
# - Plotting tasks are not included yet. (issue #519)
# - Many segmenter parameters are still unexposed. (issue #520)
# - CNLoH/Balanced Segment calling not integrated (issue #528)
#
###########

workflow case_gatk_acnv_workflow {
    String wf_entity_id_tumor
    String wf_entity_id_normal
    File target_bed
    File ref_fasta
    File ref_fasta_dict
    File ref_fasta_fai
    File common_snp_list
    File tumor_bam
    File tumor_bam_idx
    File normal_bam
    File normal_bam_idx
    File jar_file
    File PoN
    String is_disable_reference_validation
    String jni_lib    
    String plots_dir
    Boolean enable_gc_correction

  call PadTargets {
    input:
        target_bed=target_bed,
        jar_file=jar_file,
        mem=1
  }

  call CalculateTargetCoverage as TumorCalculateTargetCoverage {
    input:
        entity_id=wf_entity_id_tumor,
        padded_target_bed=PadTargets.padded_target_bed,
        input_bam=tumor_bam,
        ref_fasta=ref_fasta,
        ref_fasta_fai=ref_fasta_fai,
        ref_fasta_dict=ref_fasta_dict,
        jar_file=jar_file,
        is_disable_reference_validation=is_disable_reference_validation,
        mem=2
  }

  call AnnotateTargets {
    input:
        entity_id=wf_entity_id_tumor,
        jar_file=jar_file,
        padded_target_bed=PadTargets.padded_target_bed,
        ref_fasta=ref_fasta,
        ref_fasta_fai=ref_fasta_fai,
        ref_fasta_dict=ref_fasta_dict,
        enable_gc_correction=enable_gc_correction,
        mem=4
  }

  call CorrectGCBias as TumorCorrectGCBias {
    input:
        entity_id=wf_entity_id_tumor,
        jar_file=jar_file,
        coverage_file=TumorCalculateTargetCoverage.gatk_cnv_coverage_file,
        annotated_targets=AnnotateTargets.annotated_targets,
        enable_gc_correction=enable_gc_correction,
        mem=4
  }

  call NormalizeSomaticReadCounts as TumorNormalizeSomaticReadCounts {
    input:
        entity_id=wf_entity_id_tumor,
        coverage_file=TumorCorrectGCBias.gatk_cnv_coverage_file_gcbias,
        padded_target_bed=PadTargets.padded_target_bed,
        pon=PoN,
        jar_file=jar_file,
        mem=2,
	    jni_lib=jni_lib
  }

  call PerformSegmentation as TumorPerformSeg {
    input:
        entity_id=wf_entity_id_tumor,
        jar_file=jar_file,
        tn_file=TumorNormalizeSomaticReadCounts.tn_file,
        mem=2
  }

  call Caller as TumorCaller {
    input:
        entity_id=wf_entity_id_tumor,
        jar_file=jar_file,
        tn_file=TumorNormalizeSomaticReadCounts.tn_file,
        seg_file=TumorPerformSeg.seg_file,
        mem=2
  }

  call CalculateTargetCoverage as NormalCalculateTargetCoverage {
    input:
        entity_id=wf_entity_id_normal,
        padded_target_bed=PadTargets.padded_target_bed,
        input_bam=tumor_bam,
        ref_fasta=ref_fasta,
        ref_fasta_fai=ref_fasta_fai,
        ref_fasta_dict=ref_fasta_dict,
        jar_file=jar_file,
        is_disable_reference_validation=is_disable_reference_validation,
        mem=2
  }

  call CorrectGCBias as NormalCorrectGCBias {
    input:
        entity_id=wf_entity_id_normal,
        jar_file=jar_file,
        coverage_file=NormalCalculateTargetCoverage.gatk_cnv_coverage_file,
        annotated_targets=AnnotateTargets.annotated_targets,
        enable_gc_correction=enable_gc_correction,
        mem=4
  }

  call NormalizeSomaticReadCounts as NormalNormalizeSomaticReadCounts {
    input:
        entity_id=wf_entity_id_normal,
        coverage_file=NormalCorrectGCBias.gatk_cnv_coverage_file_gcbias,
        padded_target_bed=PadTargets.padded_target_bed,
        pon=PoN,
        jar_file=jar_file,
        mem=2,
	    jni_lib=jni_lib
  }

  call PerformSegmentation as NormalPerformSeg {
    input:
        entity_id=wf_entity_id_normal,
        jar_file=jar_file,
        tn_file=NormalNormalizeSomaticReadCounts.tn_file,
        mem=2
  }

  call Caller as NormalCaller {
    input:
        entity_id=wf_entity_id_normal,
        jar_file=jar_file,
        tn_file=NormalNormalizeSomaticReadCounts.tn_file,
        seg_file=NormalPerformSeg.seg_file,
        mem=2
  }

  call HetPulldown {
    input:
        entity_id_tumor=wf_entity_id_tumor,
        entity_id_normal=wf_entity_id_normal,
        jar_file=jar_file,
        ref_fasta=ref_fasta,
        ref_fasta_fai=ref_fasta_fai,
        ref_fasta_dict=ref_fasta_dict,
        tumor_bam=tumor_bam,
        tumor_bam_idx=tumor_bam_idx,
        normal_bam=normal_bam,
        normal_bam_idx=normal_bam_idx,
        common_snp_list=common_snp_list,
        mem = 4
  }

  call AllelicCNV {
    input:
        entity_id=wf_entity_id_tumor,
        jar_file=jar_file,
        tumor_hets=HetPulldown.tumor_hets,
        called_file=TumorCaller.called_file,
        tn_file=TumorNormalizeSomaticReadCounts.tn_file,
        mem=4
  }

  call PlotSegmentedCopyRatio {
    input:
        jar_file=jar_file,
        tn_file=TumorNormalizeSomaticReadCounts.tn_file,
        pre_tn_file=TumorNormalizeSomaticReadCounts.pre_tn_file,
        called_file=TumorCaller.called_file,
        output_dir="${plots_dir}CopyRatio_Plots",
        mem=4
  }

  call PlotACNVResults {
    input:
        jar_file=jar_file,
        tumor_hets=HetPulldown.tumor_hets,
        acnv_seg=AllelicCNV.acnv_final_segs,
        tn_file=TumorNormalizeSomaticReadCounts.tn_file,
        output_dir="${plots_dir}ACNV_plots",
        mem=4
  }
}

# Pad the target file.  This was found to help sensitivity and specificity.  This step should only be altered
#  by advanced users.  Note that by changing this, you need to have a PoN that also reflects the change.
task PadTargets {
    File target_bed
    Int padding = 250
    File jar_file
    Int mem

    command {
        java -Xmx${mem}g -jar ${jar_file} PadTargets  --targets ${target_bed} --output targets.padded.tsv \
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
    File jar_file
    String is_disable_reference_validation
    Int mem

    command {
        java -Xmx${mem}g -jar ${jar_file} CalculateTargetCoverage --output ${entity_id}.coverage.tsv \
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
    File jar_file
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    Boolean enable_gc_correction
    Int mem

    command {
        if [ ${enable_gc_correction} = true ]; \
          then java -Xmx${mem}g -jar ${jar_file} AnnotateTargets --targets ${padded_target_bed} --reference ${ref_fasta} --output ${entity_id}.annotated.tsv; \
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
    File jar_file
    Boolean enable_gc_correction
    Int mem

    command {
        if [ ${enable_gc_correction} = true ]; \
          then java -Xmx${mem}g -jar ${jar_file} CorrectGCBias --input ${coverage_file} \
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
    File jar_file
    String jni_lib
    Int mem

    command {
        java -Xmx${mem}g -Djava.library.path=${jni_lib} -jar ${jar_file} NormalizeSomaticReadCounts  --input ${coverage_file} \
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
    File jar_file
    File tn_file
    Int mem

    command {
        java -Xmx${mem}g -jar ${jar_file} PerformSegmentation  --targets ${tn_file} \
         --output ${entity_id}.seg --log2Input true  --alpha 0.01 --nperm 10000 --pmethod HYBRID --minWidth 2 --kmax 25 \
         --nmin 200 --eta 0.05 --trim 0.025 --undoSplits NONE --undoPrune 0.05 --undoSD 3 --help false --version false \
         --verbosity INFO --QUIET false
    }

    output {
        File seg_file = "${entity_id}.seg"
    }
}

# Make calls (amp, neutral, or deleted) on each segment.
task Caller {
    String entity_id
    File jar_file
    File tn_file
    File seg_file
    Int mem

    command {
        java -Xmx${mem}g -jar ${jar_file} CallSegments  --targets ${tn_file} \
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
    File jar_file
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
        java -Xmx${mem}g -jar ${jar_file} GetHetCoverage  --reference ${ref_fasta} \
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
    File jar_file
    File tumor_hets
    File called_file
    File tn_file
    Int mem

    command {
        java -Xmx${mem}g -jar ${jar_file} AllelicCNV  --tumorHets ${tumor_hets} \
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
    File jar_file
    File tn_file
    File pre_tn_file
    File called_file
    String output_dir
    Int mem

    command {
        mkdir -p ${output_dir} && \
	java -Xmx${mem}g -jar ${jar_file} PlotSegmentedCopyRatio --tangentNormalized ${tn_file} \
	 --preTangentNormalized ${pre_tn_file} --segments ${called_file} \
	 --output ${output_dir}
    }

    output {
        #corresponding output images
    }
}

#Create plots of allele fraction data
task PlotACNVResults {
    File jar_file
    File tumor_hets
    File tn_file
    File acnv_seg
    String output_dir
    Int mem

    command {
        mkdir -p ${output_dir} && \
	java -Xmx${mem}g -jar ${jar_file} PlotACNVResults --hets ${tumor_hets} \
         --tangentNormalized ${tn_file} --segments ${acnv_seg} \
         --output ${output_dir}
    }

    output {
        #corresponding output images
    }
}
