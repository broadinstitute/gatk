# Subworkflow for running GATK CNV on a single BAM. Supports both WGS and WES.
#
# Notes:
#
# - The padded target file (padded_targets) is required for the WES workflow and should be a TSV file with the column headers:
#    contig    start    stop    name
#
# - If a target file is not provided, then the WGS workflow will be run instead and the specified value of
#   wgs_bin_size (default 10000) will be used.
#
#############

import "cnv_common_tasks.wdl" as CNVTasks

workflow CNVSomaticCopyRatioBAMWorkflow {
    # Workflow input files
    File? padded_targets
    File bam
    File bam_idx
    File ref_fasta
    File ref_fasta_dict
    File ref_fasta_fai
    File cnv_panel_of_normals
    String gatk_jar

    # If no padded target file is input, then do WGS workflow
    Boolean is_wgs = select_first([padded_targets, ""]) == ""

    call CNVTasks.CollectCoverage {
        input:
            padded_targets = padded_targets,
            bam = bam,
            bam_idx = bam_idx,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            ref_fasta_dict = ref_fasta_dict,
            gatk_jar = gatk_jar
    }

    call CNVTasks.AnnotateTargets {
        input:
            entity_id = CollectCoverage.entity_id,
            targets = CollectCoverage.coverage,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            ref_fasta_dict = ref_fasta_dict,
            gatk_jar = gatk_jar
    }

    call CNVTasks.CorrectGCBias {
        input:
            entity_id = CollectCoverage.entity_id,
            coverage = CollectCoverage.coverage,
            annotated_targets = AnnotateTargets.annotated_targets,
            gatk_jar = gatk_jar
    }

    call NormalizeSomaticReadCounts {
        input:
            entity_id = CollectCoverage.entity_id,
            coverage = CorrectGCBias.corrected_coverage,
            padded_targets = AnnotateTargets.annotated_targets,
            cnv_panel_of_normals = cnv_panel_of_normals,
            gatk_jar = gatk_jar
    }

    call PerformSegmentation {
        input:
            entity_id = CollectCoverage.entity_id,
            tn_coverage = NormalizeSomaticReadCounts.tn_coverage,
            is_wgs = is_wgs,
            gatk_jar = gatk_jar
    }

    call CallSegments {
        input:
            entity_id = CollectCoverage.entity_id,
            tn_coverage = NormalizeSomaticReadCounts.tn_coverage,
            segments = PerformSegmentation.segments,
            gatk_jar = gatk_jar
    }

    call PlotSegmentedCopyRatio  {
        input:
            entity_id = CollectCoverage.entity_id,
            tn_coverage = NormalizeSomaticReadCounts.tn_coverage,
            pre_tn_coverage = NormalizeSomaticReadCounts.pre_tn_coverage,
            called_segments = CallSegments.called_segments,
            ref_fasta_dict = ref_fasta_dict,
            gatk_jar = gatk_jar
    }

    output {
        String entity_id = CollectCoverage.entity_id
        File tn_coverage = NormalizeSomaticReadCounts.tn_coverage
        File called_segments = CallSegments.called_segments
    }
}

# Perform tangent normalization (noise reduction) on the proportional coverage file
task NormalizeSomaticReadCounts {
    String entity_id
    File coverage
    File padded_targets
    File cnv_panel_of_normals
    String gatk_jar
    Int? mem

    command {
        java -Xmx${default=4 mem}g -jar ${gatk_jar} NormalizeSomaticReadCounts \
            --input ${coverage} \
            --targets ${padded_targets} \
            --panelOfNormals ${cnv_panel_of_normals} \
            --tangentNormalized ${entity_id}.tn.tsv \
            --factorNormalizedOutput ${entity_id}.fnt.tsv \
            --preTangentNormalized ${entity_id}.preTN.tsv \
            --betaHatsOutput ${entity_id}.betaHats.tsv
    }

    output {
        File tn_coverage = "${entity_id}.tn.tsv"
        File fnt_coverage = "${entity_id}.fnt.tsv"
        File pre_tn_coverage = "${entity_id}.preTN.tsv"
        File betahats_coverage = "${entity_id}.betaHats.tsv"
    }
}

# Segment the tangent-normalized coverage profile
task PerformSegmentation {
    String entity_id
    File tn_coverage
    Boolean is_wgs
    Float? seg_param_alpha
    Float? seg_param_eta
    Int? seg_param_kmax
    Int? seg_param_minWidth
    Int? seg_param_nmin
    Int? seg_param_nperm
    String? seg_param_pmethod
    Float? seg_param_trim
    Float? seg_param_undoPrune
    Int? seg_param_undoSD
    String? seg_param_undoSplits
    String gatk_jar
    Int? mem

    command <<<
        if [ ${is_wgs} = true ]
            then
                java -Xmx${default=4 mem}g -jar ${gatk_jar} PerformSegmentation \
                    --tangentNormalized ${tn_coverage} \
                    --log2Input true \
                    --alpha ${default="0.01" seg_param_alpha} \
                    --eta ${default="0.05" seg_param_eta} \
                    --kmax ${default=25 seg_param_kmax} \
                    --minWidth ${default=2 seg_param_minWidth} \
                    --nmin ${default=200 seg_param_nmin} \
                    --nperm ${default=10000 seg_param_nperm} \
                    --pmethod ${default="HYBRID" seg_param_pmethod} \
                    --trim ${default="0.025" seg_param_trim} \
                    --undoPrune ${default="0.05" seg_param_undoPrune} \
                    --undoSD ${default=3 seg_param_undoSD} \
                    --undoSplits ${default="SDUNDO" seg_param_undoSplits} \
                    --output ${entity_id}.seg
            else
                java -Xmx${default=4 mem}g -jar ${gatk_jar} PerformSegmentation \
                    --tangentNormalized ${tn_coverage} \
                    --log2Input true \
                    --alpha ${default="0.01" seg_param_alpha} \
                    --eta ${default="0.05" seg_param_eta} \
                    --kmax ${default=25 seg_param_kmax} \
                    --minWidth ${default=2 seg_param_minWidth} \
                    --nmin ${default=200 seg_param_nmin} \
                    --nperm ${default=10000 seg_param_nperm} \
                    --pmethod ${default="HYBRID" seg_param_pmethod} \
                    --trim ${default="0.025" seg_param_trim} \
                    --undoPrune ${default="0.05" seg_param_undoPrune} \
                    --undoSD ${default=2 seg_param_undoSD} \
                    --undoSplits ${default="NONE" seg_param_undoSplits} \
                    --output ${entity_id}.seg
        fi
    >>>

    output {
        File segments = "${entity_id}.seg"
    }
}

# Make calls (amplified, neutral, or deleted) on each segment
task CallSegments {
    String entity_id
    File tn_coverage
    File segments
    String gatk_jar
    Int? mem

    command {
        java -Xmx${default=4 mem}g -jar ${gatk_jar} CallSegments \
            --tangentNormalized ${tn_coverage} \
            --segments ${segments} \
            --legacy false \
            --output ${entity_id}.called
    }

    output {
        File called_segments = "${entity_id}.called"
    }
}

# Create plots of coverage data and copy-ratio estimates
task PlotSegmentedCopyRatio {
    String entity_id
    File tn_coverage
    File pre_tn_coverage
    File called_segments
    File ref_fasta_dict
    String? output_dir
    String gatk_jar
    Int? mem

    # If optional output_dir not specified, use "."
    String output_dir_ = select_first([output_dir, "."])

    command {
        mkdir -p ${output_dir_}; \
        java -Xmx${default=4 mem}g -jar ${gatk_jar} PlotSegmentedCopyRatio \
            --tangentNormalized ${tn_coverage} \
            --preTangentNormalized ${pre_tn_coverage} \
            --segments ${called_segments} \
            -SD ${ref_fasta_dict} \
            --output ${output_dir_} \
            --outputPrefix ${entity_id}
    }

    output {
        File segments_plot = "${output_dir_}/${entity_id}_FullGenome.png"
        File before_after_normalization_plot = "${output_dir_}/${entity_id}_Before_After.png"
        File before_after_cr_lim_4 = "${output_dir_}/${entity_id}_Before_After_CR_Lim_4.png"
    }
}