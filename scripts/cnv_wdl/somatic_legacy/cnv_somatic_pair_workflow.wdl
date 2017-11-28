# Workflow for running GATK CNV (and optionally, ACNV) on tumor/normal or tumor-only cases. Supports both WGS and WES.
#
# Notes:
#
# - The target file (targets) is required for the WES workflow and should be a TSV file with the column headers:
#    contig    start    stop    name
#   These targets will be padded on both sides by the amount specified by PadTargets.padding (default 250).
#
# - If a target file is not provided, then the WGS workflow will be run instead and the specified value of
#   wgs_bin_length (default 1000) will be used.
#
# - A normal BAM (normal_bam) is requrired for the tumor/normal workflow.  If not provided, the tumor-only workflow
#   will be run.
#
# - The sites file (common_sites) is required for the ACNV workflow and should be a Picard interval list.
#   If not provided, the ACNV workflow will not be run.
#
# - Example invocation:
#    java -jar cromwell.jar run cnv_somatic_pair_workflow.wdl myParameters.json
#   See cnv_somatic_pair_workflow_template.json for a template json file to modify with your own parameters (please save
#   your modified version with a different filename and do not commit to the gatk repository).
#
#############

import "cnv_common_tasks.wdl" as CNVTasks
import "cnv_somatic_copy_ratio_bam_workflow.wdl" as CopyRatio
import "cnv_somatic_allele_fraction_pair_workflow.wdl" as AlleleFraction
import "cnv_somatic_oncotate.wdl" as Oncotate

workflow CNVSomaticPairWorkflow {
    # Workflow input files
    File? targets
    File? common_sites
    File tumor_bam
    File tumor_bam_idx
    File? normal_bam
    File? normal_bam_idx
    File ref_fasta
    File ref_fasta_dict
    File ref_fasta_fai
    File cnv_panel_of_normals
    String gatk_jar

    # If no target file is input, then do WGS workflow
    Boolean is_wgs = select_first([targets, ""]) == ""
    # If no sites file is input, then do not do ACNV workflow
    Boolean is_cnv_only = select_first([common_sites, ""]) == ""
    # If no normal BAM is input, then do tumor-only workflow
    Boolean is_tumor_only = select_first([normal_bam, ""]) == ""

    Boolean is_run_oncotator = false

    # docker images
    String gatk_docker
    String oncotator_docker="broadinstitute/oncotator:1.9.3.0-eval-gatk-protected"

    if (!is_wgs) {
        call CNVTasks.PadTargets {
            input:
                # The task will fail if targets is not defined when it gets here, but that should not be allowed to happen.
                targets = select_first([targets, ""]),
                gatk_jar = gatk_jar,
                gatk_docker = gatk_docker
        }
    }

    call CopyRatio.CNVSomaticCopyRatioBAMWorkflow as TumorCopyRatioWorkflow {
        input:
            padded_targets = PadTargets.padded_targets,
            bam = tumor_bam,
            bam_idx = tumor_bam_idx,
            ref_fasta = ref_fasta,
            ref_fasta_dict = ref_fasta_dict,
            ref_fasta_fai = ref_fasta_fai,
            cnv_panel_of_normals = cnv_panel_of_normals,
            gatk_jar = gatk_jar,
            gatk_docker = gatk_docker
    }

    if (!is_tumor_only) {
        call CopyRatio.CNVSomaticCopyRatioBAMWorkflow as NormalCopyRatioWorkflow {
            input:
                padded_targets = PadTargets.padded_targets,
                bam = select_first([normal_bam, ""]),
                bam_idx = select_first([normal_bam_idx, ""]),
                ref_fasta = ref_fasta,
                ref_fasta_dict = ref_fasta_dict,
                ref_fasta_fai = ref_fasta_fai,
                cnv_panel_of_normals = cnv_panel_of_normals,
                gatk_jar = gatk_jar,
                gatk_docker = gatk_docker
        }
    }

    if (!is_cnv_only) {
        call AlleleFraction.CNVSomaticAlleleFractionPairWorkflow as TumorAlleleFractionWorkflow {
            input:
                common_sites = select_first([common_sites, ""]),
                tumor_bam = tumor_bam,
                tumor_bam_idx = tumor_bam_idx,
                normal_bam = normal_bam,    # If no normal BAM is input, tumor-only GetBayesianHetCoverage will be run
                normal_bam_idx = normal_bam_idx,
                tumor_tn_coverage = TumorCopyRatioWorkflow.tn_coverage,
                tumor_called_segments = TumorCopyRatioWorkflow.called_segments,
                ref_fasta = ref_fasta,
                ref_fasta_dict = ref_fasta_dict,
                ref_fasta_fai = ref_fasta_fai,
                gatk_jar = gatk_jar,
                gatk_docker = gatk_docker,
                is_wgs = is_wgs
        }
    }

    if (is_run_oncotator) {
        call Oncotate.CNVOncotateCalledSegments as OncotateCalledCNVWorkflow {
            input:
                 called_file=TumorCopyRatioWorkflow.called_segments,
                 oncotator_docker=oncotator_docker
        }
    }

    output {
        String tumor_entity_id = TumorCopyRatioWorkflow.entity_id
        File tumor_tn_coverage = TumorCopyRatioWorkflow.tn_coverage
        File tumor_called_segments = TumorCopyRatioWorkflow.called_segments
        String? normal_entity_id = NormalCopyRatioWorkflow.entity_id
        File? normal_tn_coverage = NormalCopyRatioWorkflow.tn_coverage
        File? normal_called_segments = NormalCopyRatioWorkflow.called_segments
        File? tumor_hets = TumorAlleleFractionWorkflow.tumor_hets
        File? tumor_acnv_segments = TumorAlleleFractionWorkflow.acnv_segments
        File? oncotated_called_file = OncotateCalledCNVWorkflow.oncotated_called_file
    }
}