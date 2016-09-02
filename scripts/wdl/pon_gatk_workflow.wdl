#
# A workflow for creating a Panel of Normals given a list of normal samples. Supports both WGS and WES samples. This was tested on a3c7368 commit of gatk-protected.
#
# Notes:
#
# - Input file (normal_bam_list) must contain file paths to bam and bam index files separated by tabs in the following format:
#   normal_bam_1	bam_idx_1
#   normal_bam_2	bam_idx_2
#   ...
# 
# - set isWGS variable to true or false to specify whether to run a WGS or WES workflow respectively
#
# - target file (which must be in tsv format) is only used with WES workflow - WGS workflow generates its own targets (so user can pass any string as an argument)
#
# - THIS SCRIPT SHOULD BE CONSIDERED OF "BETA" QUALITY
#
# - Example invocation:
#    java -jar cromwell.jar run pon_gatk_workflow.wdl myParameters.json
# - See pon_gatk_workflow_template.json for a template json file to modify with your own parameters (please save
#    your modified version with a different filename and do not commit to gatk-protected repo).
#
# - Some call inputs might seem out of place - consult with the comments in task definitions for details
#
#############

workflow pon_gatk_workflow {
    # Workflow input files
    File normal_bam_list
    Array[Array[String]] bam_file_names = read_tsv(normal_bam_list)
    String target_file 
    File ref_fasta
    File ref_fasta_dict
    File ref_fasta_fai
    File gatk_jar

    # Workflow options 
    Boolean enable_gc_correction
    Boolean isWGS
    Boolean disable_sequence_dictionary_validation
    String combined_entity_id = "combined_read_counts"
    Int max_open_files
    Int wgsBinSize
    Boolean noQC

    # PoN name
    String pon_entity_id
    
  call PadTargets {
    input:
        target_file=target_file,
        gatk_jar=gatk_jar,
        isWGS=isWGS,
        mem=1
  }

  scatter (row in bam_file_names) {

    call GetBamFileName {
      input:
          input_bam=row[0]
    }

    call CalculateTargetCoverage {
      input:
          entity_id=GetBamFileName.name,
          padded_target_file=PadTargets.padded_target_file,
          input_bam=row[0],
          bam_idx=row[1],
          ref_fasta=ref_fasta,
          ref_fasta_fai=ref_fasta_fai,
          ref_fasta_dict=ref_fasta_dict,
          gatk_jar=gatk_jar,
          disable_sequence_dictionary_validation=disable_sequence_dictionary_validation,
          isWGS=isWGS,
          mem=2
    }

    call WholeGenomeCoverage {
      input:
          entity_id=GetBamFileName.name,
          target_file=PadTargets.padded_target_file,
          input_bam=row[0],
          bam_idx=row[1],
          coverage_file=CalculateTargetCoverage.gatk_coverage_file,
          ref_fasta=ref_fasta,
          ref_fasta_fai=ref_fasta_fai,
          ref_fasta_dict=ref_fasta_dict,
          gatk_jar=gatk_jar,
          isWGS=isWGS,
          wgsBinSize=wgsBinSize,
          mem=4
    }
  }

  call AggregateTargetFiles {
    input:
        target_file_list=WholeGenomeCoverage.gatk_target_file
  }

  call CombineReadCounts {
    input:
        combined_entity_id=combined_entity_id,
        coverage_file_list=WholeGenomeCoverage.gatk_coverage_file,
        gatk_jar=gatk_jar,
        max_open_files=max_open_files,
        mem=4
  }

  call AnnotateTargets {
    input:
        entity_id=combined_entity_id,
        gatk_jar=gatk_jar,
        target_file=AggregateTargetFiles.target_file,
        ref_fasta=ref_fasta,
        ref_fasta_fai=ref_fasta_fai,
        ref_fasta_dict=ref_fasta_dict,
        enable_gc_correction=enable_gc_correction,
        mem=4
  }

  call CorrectGCBias {
    input:
        entity_id=combined_entity_id,
        gatk_jar=gatk_jar,
        coverage_file=CombineReadCounts.combined_read_counts,
        annotated_targets=AnnotateTargets.annotated_targets,
        enable_gc_correction=enable_gc_correction,
        mem=4
  }

  call CreatePanelOfNormals {
    input:
        pon_entity_id=pon_entity_id,
        read_counts_file=CorrectGCBias.coverage_file_gcbias_corrected,
        gatk_jar=gatk_jar,
        noQC=noQC,
        mem=4
  }
}

# Pad the target file. This was found to help sensitivity and specificity. This step should only be altered
# by advanced users. Note that by changing this, you need to have a PoN that also reflects the change.
task PadTargets {
    String target_file
    Int padding = 250
    File gatk_jar
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

# Helper task to get the name of the given bam file
task GetBamFileName {
    File input_bam

    command <<<
        echo $(basename "${input_bam}" .bam)
     >>>

    output {
        String name=read_string(stdout())
    }
}

# Calculate the target proportional coverage
task CalculateTargetCoverage {
    String entity_id
    File padded_target_file
    String transform = "PCOV"
    String grouping = "SAMPLE"
    Boolean keep_duplicate_reads = true
    Boolean disable_all_read_filters = false
    File input_bam
    File bam_idx
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    File gatk_jar
    Boolean disable_sequence_dictionary_validation
    Boolean isWGS
    Int mem

    # Note that when isWGS is true, this task is still called by the workflow.
    # In that case, an empty coverage file is created and passed to the WholeGenomeCoverage 
    # task to satisfy input and output requirements.
    command {
        if [ ${isWGS} = false ]; \
          then java -Xmx${mem}g -jar ${gatk_jar} CalculateTargetCoverage --output ${entity_id}.coverage.tsv \
                --groupBy ${grouping} --transform ${transform} --targets ${padded_target_file} --targetInformationColumns FULL \
                --keepduplicatereads ${keep_duplicate_reads} --input ${input_bam} --reference ${ref_fasta} \
                --disable_all_read_filters ${disable_all_read_filters} --interval_set_rule UNION --interval_padding 0 \
                --secondsBetweenProgressUpdates 10.0 --disableSequenceDictionaryValidation ${disable_sequence_dictionary_validation} \
                --createOutputBamIndex true --help false --version false --verbosity INFO --QUIET false; \
          else touch ${entity_id}.coverage.tsv; \
        fi
    }

    output {
        File gatk_coverage_file = "${entity_id}.coverage.tsv"
    }

    #runtime {
    #    docker: "gatk-protected/a1"
    #}
}

# Calculate coverage on Whole Genome Sequence using Spark. This task automatically creates an
# output file, containing targets
task WholeGenomeCoverage {
    String entity_id
    File coverage_file 
    File target_file
    File input_bam
    File bam_idx
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    File gatk_jar
    Boolean isWGS
    Int wgsBinSize
    Int mem

    # If isWGS is set to true, the task produces WGS coverage and targets that are passed to downstream tasks
    # If not, coverage and target files (received from upstream) for WES are passed downstream
    command {
        if [ ${isWGS} = true ]; \
          then java -Xmx${mem}g -jar ${gatk_jar} SparkGenomeReadCounts --outputFile ${entity_id}.coverage.tsv \
                --reference ${ref_fasta} --input ${input_bam} --binsize ${wgsBinSize}; \
          else ln -s ${coverage_file} ${entity_id}.coverage.tsv; ln -s ${target_file} ${entity_id}.coverage.tsv.targets.tsv; \
        fi
    }

    output {
        File gatk_coverage_file = "${entity_id}.coverage.tsv"
        File gatk_target_file = "${entity_id}.coverage.tsv.targets.tsv"
    }
}

# Helper task to aggregate the array of targets that is output by the scatter block into a single file
task AggregateTargetFiles {
    Array[File]+ target_file_list

    command <<<
        ln -s ${target_file_list[0]} targets.tsv
    >>>

    output {
        File target_file = "targets.tsv"
    }
}

# Combine a set of read-count input files into a single multicolumn output file
task CombineReadCounts {
    String combined_entity_id
    Array[File] coverage_file_list
    Int max_open_files
    File gatk_jar
    Int mem

    command {
        java -Xmx${mem}g -jar ${gatk_jar} CombineReadCounts --input ${sep=" --input " coverage_file_list} \
         --maxOpenFiles ${max_open_files} --output ${combined_entity_id}.tsv
    }

    output {
        File combined_read_counts = "${combined_entity_id}.tsv"
    }
}

# Add new columns to an existing target table with various targets
# Note that this task is optional
task AnnotateTargets {
    String entity_id
    File target_file
    File gatk_jar
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
    File gatk_jar
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
        File coverage_file_gcbias_corrected = "${entity_id}.gc_corrected_coverage.tsv"
    }
}

# Create panel of normals given a collection of read counts for control samples
task CreatePanelOfNormals {
    String pon_entity_id
    File read_counts_file
    File gatk_jar
    Boolean noQC
    Int mem

    command {
        # If there are no removed samples the output file still needs to be created
        touch "${pon_entity_id}.pon.removed_samples.txt"
        java -Xmx${mem}g -jar ${gatk_jar} CreatePanelOfNormals --extremeColumnMedianCountPercentileThreshold 2.5 \
         --truncatePercentileThreshold 0.1 --input ${read_counts_file} --output ${pon_entity_id}.pon \
         --noQC ${noQC}
    }

    output {
        File output_pon = "${pon_entity_id}.pon"
        File removed_samples = "${pon_entity_id}.pon.removed_samples.txt"
        File target_weights = "${pon_entity_id}.pon.target_weights.txt"
    }
}
