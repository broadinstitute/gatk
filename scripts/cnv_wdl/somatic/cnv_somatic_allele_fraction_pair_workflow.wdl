# Subworkflow for running GATK ACNV on tumor/normal or tumor-only cases. Supports both WGS and WES.
#
# Notes:
#
# - A normal BAM (normal_bam) is requrired for the tumor/normal workflow.  If not provided, the tumor-only workflow
#   will be run.
#
# - ACNV is run only on the tumor.
#
# - If is_wgs is true, then more aggressive segment-merging parameters will be used for ACNV by default.
#   These parameter choices can be overrided manually in the json parameter file for either WES or WGS, if desired.  
#
#############

import "cnv_somatic_tasks.wdl" as CNVSomatic

workflow CNVSomaticAlleleFractionPairWorkflow {
    # Workflow input files
    File common_sites
    File tumor_bam
    File tumor_bam_idx
    File? normal_bam
    File? normal_bam_idx
    File tumor_tn_coverage
    File tumor_called_segments
    File ref_fasta
    File ref_fasta_dict
    File ref_fasta_fai
    String gatk_jar

    # If WGS, different segment-merging parameters will be used for ACNV
    Boolean is_wgs

    call GetBayesianHetCoverage {
        input:
            common_sites = common_sites,
            tumor_bam = tumor_bam,
            tumor_bam_idx = tumor_bam_idx,
            normal_bam = normal_bam,
            normal_bam_idx = normal_bam_idx,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            ref_fasta_dict = ref_fasta_dict,
            gatk_jar = gatk_jar
    }

    call AllelicCNV {
        input:
            entity_id = GetBayesianHetCoverage.tumor_entity_id,
            hets = GetBayesianHetCoverage.tumor_hets,
            tn_coverage = tumor_tn_coverage,
            called_segments = tumor_called_segments,
            is_wgs = is_wgs,
            gatk_jar = gatk_jar
    }

    call PlotACNVResults {
        input:
            entity_id = GetBayesianHetCoverage.tumor_entity_id,
            hets = GetBayesianHetCoverage.tumor_hets,
            tn_coverage = tumor_tn_coverage,
            acnv_segments = AllelicCNV.acnv_segments,
            ref_fasta_dict = ref_fasta_dict,
            gatk_jar = gatk_jar
    }

    call ConvertACNVResults {
        input:
            entity_id = GetBayesianHetCoverage.tumor_entity_id,
            hets = GetBayesianHetCoverage.tumor_hets,
            tn_coverage = tumor_tn_coverage,
            acnv_segments = AllelicCNV.acnv_segments,
            gatk_jar = gatk_jar
    }

    output {
        String tumor_entity_id = GetBayesianHetCoverage.tumor_entity_id
        File tumor_hets = GetBayesianHetCoverage.tumor_hets
        File acnv_segments = AllelicCNV.acnv_segments
    }
}

# Call heterozygous SNPs in the normal and then collect allele counts in the tumor for each called position
# Entity IDs can be the same value and are derived from the BAMs if not provided
# If the normal BAM is not provided, then tumor-only mode is run
task GetBayesianHetCoverage {
    File common_sites
    File tumor_bam
    File tumor_bam_idx
    File? normal_bam
    File? normal_bam_idx
    Float? stringency
    Int? read_depth_threshold
    File ref_fasta
    File ref_fasta_dict
    File ref_fasta_fai
    String gatk_jar
    Int? mem

    # Sample names are derived from the bam filenames
    String tumor_base_filename = sub(sub(sub(tumor_bam, "gs://", ""), "[/]*.*/", ""), "\\.bam$", "")
    String normal_base_filename = sub(sub(sub(select_first([normal_bam, ""]), "gs://", ""), "[/]*.*/", ""), "\\.bam$", "") # Evaluates to "" if tumor-only and is not used

    # TODO the below if-else block can be cleaned up once Cromwell support for optionals/conditionals improves
    Boolean is_tumor_only = normal_base_filename == ""

    command <<<
        if [ ${is_tumor_only} = true ]
            then
                java -Xmx${default=4 mem}g -jar ${gatk_jar} GetBayesianHetCoverage \
                    --reference ${ref_fasta} \
                    --tumor ${tumor_bam} \
                    --tumorHets ${tumor_base_filename + ".tumor.hets.tsv"} \
                    --snpIntervals ${common_sites} \
                    --hetCallingStringency ${default=30 stringency} \
                    --readDepthThreshold ${default=15 read_depth_threshold} \
                    --VALIDATION_STRINGENCY LENIENT
            else
                java -Xmx${default=4 mem}g -jar ${gatk_jar} GetBayesianHetCoverage \
                    --reference ${ref_fasta} \
                    --tumor ${tumor_bam} \
                    --tumorHets ${tumor_base_filename + ".tumor.hets.tsv"} \
                    ${"--normal " + normal_bam} \
                    --normalHets ${normal_base_filename + ".normal.hets.tsv"} \
                    --snpIntervals ${common_sites} \
                    --hetCallingStringency ${default=30 stringency} \
                    --readDepthThreshold ${default=15 read_depth_threshold} \
                    --VALIDATION_STRINGENCY LENIENT
        fi
    >>>

    output {
        String tumor_entity_id = tumor_base_filename
        File tumor_hets = "${tumor_base_filename}.tumor.hets.tsv"
    }
}

# Estimate minor-allele fraction and revise copy-ratio segments
task AllelicCNV {
    String entity_id
    File hets
    File tn_coverage
    File called_segments
    Boolean is_wgs
    String gatk_jar
    Boolean? use_all_copy_ratio_segments
    Int? max_num_iterations_snp_seg
    Int? small_segment_threshold
    Int? num_samples_copy_ratio
    Int? num_burn_in_copy_ratio
    Int? num_samples_allele_fraction
    Int? num_burn_in_allele_fraction
    Float? interval_threshold_copy_ratio
    Float? interval_threshold_allele_fraction
    Int? num_iterations_sim_seg_per_fit
    Int? max_num_iterations_sim_seg
    Int? mem

    command <<<
        if [ ${is_wgs} = true ]
            then
                java -Xmx${default=4 mem}g -jar ${gatk_jar} AllelicCNV \
                    --tumorHets ${hets} \
                    --tangentNormalized ${tn_coverage} \
                    --segments ${called_segments} \
                    --outputPrefix ${entity_id} \
                    --useAllCopyRatioSegments ${default=false use_all_copy_ratio_segments} \
                    --maxNumIterationsSNPSeg ${default=5 max_num_iterations_snp_seg} \
                    --smallSegmentThreshold ${default=3 small_segment_threshold} \
                    --numSamplesCopyRatio ${default=50 num_samples_copy_ratio} \
                    --numBurnInCopyRatio ${default=25 num_burn_in_copy_ratio} \
                    --numSamplesAlleleFraction ${default=50 num_samples_allele_fraction} \
                    --numBurnInAlleleFraction ${default=25 num_burn_in_allele_fraction} \
                    --intervalThresholdCopyRatio ${default="4.0" interval_threshold_copy_ratio} \
                    --intervalThresholdAlleleFraction ${default="2.0" interval_threshold_allele_fraction} \
                    --numIterationsSimSegPerFit ${default=5 num_iterations_sim_seg_per_fit} \
                    --maxNumIterationsSimSeg ${default=25 max_num_iterations_sim_seg}
            else
                java -Xmx${default=4 mem}g -jar ${gatk_jar} AllelicCNV \
                    --tumorHets ${hets} \
                    --tangentNormalized ${tn_coverage} \
                    --segments ${called_segments} \
                    --outputPrefix ${entity_id} \
                    --useAllCopyRatioSegments ${default=false use_all_copy_ratio_segments} \
                    --maxNumIterationsSNPSeg ${default=25 max_num_iterations_snp_seg} \
                    --smallSegmentThreshold ${default=3 small_segment_threshold} \
                    --numSamplesCopyRatio ${default=100 num_samples_copy_ratio} \
                    --numBurnInCopyRatio ${default=50 num_burn_in_copy_ratio} \
                    --numSamplesAlleleFraction ${default=100 num_samples_allele_fraction} \
                    --numBurnInAlleleFraction ${default=50 num_burn_in_allele_fraction} \
                    --intervalThresholdCopyRatio ${default="4.0" interval_threshold_copy_ratio} \
                    --intervalThresholdAlleleFraction ${default="2.0" interval_threshold_allele_fraction} \
                    --numIterationsSimSegPerFit ${default=1 num_iterations_sim_seg_per_fit} \
                    --maxNumIterationsSimSeg ${default=10 max_num_iterations_sim_seg}
        fi
    >>>

    output {
        File acnv_segments = "${entity_id}-sim-final.seg"
    }
}

# Create plots of allele-count data and minor-allele-fraction estimates
task PlotACNVResults {
    String entity_id
    File hets
    File tn_coverage
    File acnv_segments
    File ref_fasta_dict
    String? output_dir
    String gatk_jar
    Int? mem

    # If optional output_dir not specified, use "."
    String output_dir_ = select_first([output_dir, "."])

    command {
        mkdir -p ${output_dir_}; \
        java -Xmx${default=4 mem}g -jar ${gatk_jar} PlotACNVResults \
            --hets ${hets} \
            --tangentNormalized ${tn_coverage} \
            --segments ${acnv_segments} \
            -SD ${ref_fasta_dict} \
            --output ${output_dir_} \
            --outputPrefix ${entity_id}
    }

    output {
        File acnv_plot = "${output_dir_}/${entity_id}_ACNV.png"
    }
}

# Convert ACNV results to GATK CNV, TITAN, and ACS formats
task ConvertACNVResults {
    String entity_id
    File hets
    File tn_coverage
    File acnv_segments
    String? output_dir
    String gatk_jar
    Int? mem

    # If optional output_dir not specified, use "."
    String output_dir_ = select_first([output_dir, "."])

    command {
        mkdir -p ${output_dir_}; \
        java -Xmx${default=4 mem}g -jar ${gatk_jar} ConvertACNVResults \
            --tumorHets ${hets} \
            --tangentNormalized ${tn_coverage} \
            --segments ${acnv_segments} \
            --outputDir ${output_dir_}
    }

    output {
        File acs_segments = "${output_dir_}/${entity_id}-sim-final.acs.seg"
        File cnb_called_segments = "${output_dir_}/${entity_id}-sim-final.cnb_called.seg"
        File cnv_segments = "${output_dir_}/${entity_id}-sim-final.cnv.seg"
        File titan_hets = "${output_dir_}/${entity_id}-sim-final.titan.het.tsv"
        File titan_tn_coverage = "${output_dir_}/${entity_id}-sim-final.titan.tn.tsv"
    }
}