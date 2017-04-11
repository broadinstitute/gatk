# Workflow for creating a GATK CNV Panel of Normals given a list of normal samples. Supports both WGS and WES.
#
# Notes:
#
# - Input file (normal_bams_list) must contain file paths to bam and bam index files separated by tabs in the following format:
#    normal_bam_1    bam_idx_1
#    normal_bam_2    bam_idx_2
#    ...
#
# - The target file (targets) is required for the WES workflow and should be a TSV file with the column headers:
#    contig    start    stop    name
#   These targets will be padded on both sides by the amount specified by PadTargets.padding (default 250).
#
# - If a target file is not provided, then the WGS workflow will be run instead and the specified value of
#   wgs_bin_size (default 10000) will be used.
#
# - Example invocation:
#    java -jar cromwell.jar run cnv_somatic_panel_workflow.wdl myParameters.json
#   See cnv_somatic_panel_workflow_template.json for a template json file to modify with your own parameters (please save
#   your modified version with a different filename and do not commit to the gatk-protected repository).
#
#############

import "cnv_somatic_tasks.wdl" as CNVSomatic

workflow CNVSomaticPanelWorkflow {
    # Workflow input files
    File? targets
    File normal_bams_list
    Array[Array[String]]+ normal_bams = read_tsv(normal_bams_list)
    File ref_fasta
    File ref_fasta_dict
    File ref_fasta_fai
    File gatk_jar

    # CombineReadCounts name
    String combined_entity_id = "combined_coverage"
    # PoN name
    String pon_entity_id

    # If no target file is input, then do WGS workflow
    Boolean is_wgs = select_first([targets, ""]) == ""

    if (!is_wgs) {
        call CNVSomatic.PadTargets {
            input:
                targets = targets,
                gatk_jar = gatk_jar
        }
    }

    scatter (normal_bam in normal_bams) {
        call CNVSomatic.CollectCoverage {
            input:
                padded_targets = PadTargets.padded_targets,
                bam = normal_bam[0],
                bam_idx = normal_bam[1],
                ref_fasta = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,
                ref_fasta_dict = ref_fasta_dict,
                gatk_jar = gatk_jar
        }
    }

    call CombineReadCounts {
        input:
            combined_entity_id = combined_entity_id,
            coverage_file_list = CollectCoverage.coverage,
            gatk_jar = gatk_jar
    }

    call CNVSomatic.AnnotateTargets {
        input:
            entity_id = combined_entity_id,
            targets = CollectCoverage.coverage[0],
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            ref_fasta_dict = ref_fasta_dict,
            gatk_jar = gatk_jar
    }

    call CNVSomatic.CorrectGCBias {
        input:
            entity_id = combined_entity_id,
            coverage = CombineReadCounts.combined_coverage,
            annotated_targets = AnnotateTargets.annotated_targets,
            gatk_jar = gatk_jar
    }

    call CreatePanelOfNormals {
        input:
            pon_entity_id = pon_entity_id,
            corrected_coverage = CorrectGCBias.corrected_coverage,
            gatk_jar = gatk_jar
    }

    output {
        File output_pon = CreatePanelOfNormals.output_pon
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

# Create panel of normals from combined, GC-corrected coverage profiles
task CreatePanelOfNormals {
    String pon_entity_id
    File corrected_coverage
    Boolean? no_qc
    File gatk_jar
    Int? mem

    command {
        # If there are no removed samples the output file still needs to be created
        touch "${pon_entity_id}.pon.removed_samples.txt" ; \
        java -Xmx${default=4 mem}g -jar ${gatk_jar} CreatePanelOfNormals \
            --input ${corrected_coverage} \
            --extremeColumnMedianCountPercentileThreshold 2.5 \
            --truncatePercentileThreshold 0.1 \
            --noQC ${default="false" no_qc} \
            --output ${pon_entity_id}.pon
    }

    output {
        File output_pon = "${pon_entity_id}.pon"
        File removed_samples = "${pon_entity_id}.pon.removed_samples.txt"
        File target_weights = "${pon_entity_id}.pon.target_weights.txt"
    }
}