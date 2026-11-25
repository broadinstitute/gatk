version 1.0

## Copyright Broad Institute, 2017
##
## This WDL workflow runs GATK4 Mutect 2 on a single tumor-normal pair or on a single tumor sample,
## and performs additional filtering.
##
## Main requirements/expectations :
## - One analysis-ready BAM or CRAM file for each sample, with its index beside it (see tumor_reads below)
##
## Description of inputs:
##
## ** Runtime **
## gatk_docker: docker image to use for GATK 4 Mutect2
## preemptible: how many preemptions to tolerate before switching to a non-preemptible machine (on Google)
## max_retries: how many times to retry failed tasks -- very important on the cloud when there are transient errors
## gatk_override: (optional) local file or Google bucket path to a GATK 4 java jar file to be used instead of the GATK 4 jar
##                in the docker image.  This must be supplied when running in an environment that does not support docker
##                (e.g. SGE cluster on a Broad on-prem VM)
##
## ** Workflow options **
## intervals: genomic intervals (will be used for scatter)
## scatter_count: number of parallel jobs to generate when scattering over intervals
## m2_extra_args, m2_extra_filtering_args: additional arguments for Mutect2 calling and filtering (optional)
## split_intervals_extra_args: additional arguments for splitting intervals before scattering (optional)
## run_orientation_bias_mixture_model_filter: (optional) if true, filter orientation bias sites with the read orientation artifact mixture model.
##
## ** Primary inputs **
## ref_fasta: reference genome FASTA.  Its companion files are not inputs: GATK finds them by name in the FASTA's directory
##     and accepts no other names, so they must exist there under these names:
##     - its index, the FASTA's path plus .fai
##     - its sequence dictionary, the FASTA's path with its extension (.fasta, .fa, .fna, .fas or .txt, optionally
##       followed by .gz) replaced by .dict
##     - for a bgzipped FASTA (.gz), its block index, the FASTA's path plus .gzi
##     For example Homo_sapiens_assembly38.fasta, Homo_sapiens_assembly38.fasta.fai and Homo_sapiens_assembly38.dict.
##     The workflow passes these files to tasks along with the FASTA, so that backends that localize files put them
##     beside it.  (Earlier versions took the index and dictionary as the inputs ref_fai and ref_dict.)
## tumor_reads: the tumor sample's reads, as an array whose first element is a BAM or CRAM file and whose optional second
##     element is its index, eg ["gs://bucket/tumor.bam"] or ["/data/tumor.bam", "/data/tumor.bai"].  Tasks use
##     tumor_reads[0] as the reads, and take the whole array so that the index travels with them.  GATK finds the index
##     only by name in the reads' directory, so it must be there under one of these names, whether or not the array
##     includes it:
##     - for x.bam: x.bai, x.bam.bai, x.csi or x.bam.csi
##     - for x.cram: x.crai or x.cram.crai
##     Whether to include the index depends on the backend:
##     - On Terra and other Google Cloud backends, the tasks stream the reads from the bucket rather than localizing them,
##       and GATK finds the index there, so the BAM or CRAM alone suffices.  In a Terra configuration, an attribute
##       expression such as this.tumor_bam gives a one-element array, and [this.tumor_bam, this.tumor_bai] a two-element
##       one.
##     - On backends that localize files, such as a local Cromwell or a cluster with a shared filesystem, include the
##       index, so that it is localized beside the reads.  Otherwise Mutect2 may fail for lack of an index.
##     Including the index always works.  (Earlier versions took it as the separate input tumor_reads_index.)
## normal_reads: the normal sample's reads, likewise, or [] for tumor-only calling.  In a Terra configuration,
##     this.normal_bam or [this.normal_bam, this.normal_bai] gives [] if the attributes are unset.  (Earlier versions took
##     the index as the separate input normal_reads_index, and the reads as an optional input that was left unset for
##     tumor-only calling.)
##
## ** Primary resources ** (optional but strongly recommended)
## pon: optional panel of normals in VCF format containing probable technical artifacts (false positves)
## germline_resource: optional database of known germline variants (see http://gnomad.broadinstitute.org/downloads)
## common_variants_vcf: VCF of common variants with allele frequencies, at whose sites Mutect2 tabulates pileup summaries,
##                      for calculating contamination
##
## ** VCF indices ** (usually unnecessary)
## pon_idx, germline_resource_idx, common_variants_vcf_idx, etc: each VCF's index is assumed to be the VCF's path plus .tbi
##     for a bgzipped VCF (.vcf.gz or .vcf.bgz) or .idx for a plain .vcf.  Specify an index only to override that, for
##     example for an index with a non-standard extension.
##
## ** Secondary resources ** (for optional tasks)
## realignment_index_bundle: resource for FilterAlignmentArtifacts, which runs if and only if it is specified.  Generated by BwaMemIndexImageCreator.
##
## Outputs :
## - One VCF file and its index with primary filtering applied; secondary filtering if requested; a bamout.bam
##   file of reassembled reads if requested
## - The F1R2 counts of the scatters (f1r2_counts), a tar.gz for each, if run_orientation_bias_mixture_model_filter
##   is true, for learning a read orientation model from them, even when skip_filtering is true
##
## Cromwell version support
## - Successfully tested on v34
##
## LICENSING :
## This script is released under the GATK source code license (Apache 2.0) (see LICENSE in
## https://github.com/broadinstitute/gatk). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script.  Please see the docker
## pages at https://hub.docker.com/r/broadinstitute/* for detailed licensing information
## pertaining to the included programs.


struct Runtime {
    String gatk_docker
    File? gatk_override
    Int max_retries
    Int preemptible
    Int cpu
    Int machine_mem
    Int command_mem
    Int disk
    Int boot_disk_size
}

workflow Mutect2 {
    input {
        # basic inputs
        File? intervals
        File? masked_intervals
        # the reference FASTA, whose index and dictionary must be beside it under the names GATK expects (see the
        # derivation of reference below)
        File ref_fasta
        # The tumor's reads: a BAM or CRAM, optionally followed by its index, eg ["gs://bucket/tumor.bam"] on Terra or
        # ["/data/tumor.bam", "/data/tumor.bai"] on a backend that localizes files, which needs the index.  The index must
        # be beside the reads under a name GATK looks for, such as tumor.bai or tumor.bam.bai (see tumor_reads in the
        # header).  Tasks use tumor_reads[0] as the reads.
        Array[File]+ tumor_reads
        # the normal's reads, likewise, or [] for tumor-only calling
        Array[File] normal_reads

        # VCFs: optional but usually recommended resources (germline_resource, pon, common_variants_vcf), alleles to
        # force-call (gga_vcf), and truth VCFs for Permutect datasets
        File? germline_resource
        File? pon
        File? gga_vcf
        File? common_variants_vcf
        File? permutect_training_dataset_truth_vcf
        File? permutect_test_dataset_truth_vcf

        # Indices of the VCFs above, usually unnecessary on both local and cloud backends.  Each VCF's index is assumed to
        # be the VCF's path plus .tbi for a bgzipped VCF (.vcf.gz or .vcf.bgz) or .idx for a plain .vcf, and is passed to
        # tasks along with the VCF.  Specify an index only to override that, for example for a non-standard extension.
        File? germline_resource_idx
        File? pon_idx
        File? gga_vcf_idx
        File? common_variants_vcf_idx
        File? permutect_training_dataset_truth_vcf_idx
        File? permutect_test_dataset_truth_vcf_idx

        # DEPRECATED: the former names of germline_resource, germline_resource_idx, common_variants_vcf and
        # common_variants_vcf_idx, kept so that existing configurations still work, and to be removed in a future release.
        # Each is used only if the input that replaced it is not given.  WDL 1.0 has no formal deprecation annotation.
        File? gnomad
        File? gnomad_idx
        File? variants_for_contamination
        File? variants_for_contamination_idx

        # extra arguments
        String? m2_extra_args
        String? m2_extra_filtering_args
        String? split_intervals_extra_args

        # Ranges of population allele frequency of the common variants: the sites where Mutect2 tabulates pileup
        # summaries, and the subset that CalculateContamination uses.  The maximum for the pileups defaults to the
        # maximum for contamination.
        Float pileups_min_population_af = 0.01
        Float? pileups_max_population_af
        Float contamination_min_population_af = 0.01
        Float contamination_max_population_af = 0.2

        # additional modes and outputs
        File? realignment_index_bundle
        String? realignment_extra_args
        File? dragstr_model
        Boolean run_orientation_bias_mixture_model_filter = false
        Boolean make_bamout = false
        Boolean compress_vcfs = false
        Boolean make_permutect_training_dataset = false
        Boolean make_permutect_test_dataset = false
        Boolean skip_filtering = false
        # skip CalculateContamination even when filtering, eg when another tool estimates contamination from the pileups
        Boolean skip_calculate_contamination = false


        # runtime
        String gatk_docker
        File? gatk_override
        String basic_bash_docker = "ubuntu:16.04"
        Int scatter_count
        Int preemptible = 2
        Int max_retries = 1
        Int small_task_cpu = 2
        Int small_task_mem = 4
        Int small_task_disk = 100
        Int boot_disk_size = 12
        Int learn_read_orientation_mem = 8000
        Int filter_alignment_artifacts_mem = 9000
        String? gcs_project_for_requester_pays

        # Use as a last resort to increase the disk given to every task in case of ill behaving data
        Int emergency_extra_disk = 0
    }

    # the inputs, or their deprecated former names if they are not given
    File? germline_resource_or_alias = if defined(germline_resource) then germline_resource else gnomad
    File? germline_resource_idx_or_alias = if defined(germline_resource_idx) then germline_resource_idx else gnomad_idx
    File? common_variants_vcf_or_alias = if defined(common_variants_vcf) then common_variants_vcf else variants_for_contamination
    File? common_variants_vcf_idx_or_alias = if defined(common_variants_vcf_idx) then common_variants_vcf_idx else variants_for_contamination_idx

    # Each VCF's index, unless given, is the VCF's path plus .tbi for a bgzipped VCF (.vcf.gz or .vcf.bgz) or .idx for a
    # plain .vcf.  Tasks take both, so that a backend that localizes files puts each index beside its VCF, while on Terra
    # both are streamed.  WDL has no user-defined functions, so we scatter this rule over the VCFs; a scatter without calls
    # only evaluates declarations.  The arrays hold paths, with "" for a file that isn't given, because Cromwell can't
    # evaluate an Array[File?] with undefined elements.  Keep the two arrays and the indices below in the same order.
    Array[String] vcf_paths = ["~{germline_resource_or_alias}", "~{pon}", "~{gga_vcf}", "~{common_variants_vcf_or_alias}",
        "~{permutect_training_dataset_truth_vcf}", "~{permutect_test_dataset_truth_vcf}"]
    Array[String] given_vcf_idx_paths = ["~{germline_resource_idx_or_alias}", "~{pon_idx}", "~{gga_vcf_idx}", "~{common_variants_vcf_idx_or_alias}",
        "~{permutect_training_dataset_truth_vcf_idx}", "~{permutect_test_dataset_truth_vcf_idx}"]
    scatter (n in range(length(vcf_paths))) {
        String vcf_idx_path = if given_vcf_idx_paths[n] != "" then given_vcf_idx_paths[n]
            else if vcf_paths[n] == "" then ""
            else vcf_paths[n] + (if sub(vcf_paths[n], "\\.b?gz$", "") != vcf_paths[n] then ".tbi" else ".idx")
    }

    # WDL 1.0 has no literal for an undefined optional, but a declaration inside a conditional is undefined outside it
    # when the condition is false
    if (false) {
        File NONE = ""
    }
    File? germline_resource_idx_or_conventional = if vcf_idx_path[0] == "" then NONE else vcf_idx_path[0]
    File? pon_idx_or_conventional = if vcf_idx_path[1] == "" then NONE else vcf_idx_path[1]
    File? gga_vcf_idx_or_conventional = if vcf_idx_path[2] == "" then NONE else vcf_idx_path[2]
    File? common_variants_vcf_idx_or_conventional = if vcf_idx_path[3] == "" then NONE else vcf_idx_path[3]
    File? permutect_training_dataset_truth_vcf_idx_or_conventional = if vcf_idx_path[4] == "" then NONE else vcf_idx_path[4]
    File? permutect_test_dataset_truth_vcf_idx_or_conventional = if vcf_idx_path[5] == "" then NONE else vcf_idx_path[5]

    # The reference FASTA and its companion files, which GATK finds by name in the FASTA's directory, accepting no other
    # names:
    # - its index, the FASTA's path plus .fai
    # - its sequence dictionary, the FASTA's path with its extension (.fasta, .fa, .fna, .fas or .txt, optionally followed
    #   by .gz) replaced by .dict
    # - for a bgzipped FASTA (.gz), its block index, the FASTA's path plus .gzi
    # Since only these names work, the workflow derives them rather than taking them as inputs.  Tasks take the files as
    # one array, so that a backend that localizes files puts them beside the FASTA, while on Terra the tasks that mark the
    # array localization_optional stream all of them from the bucket.  Keep the order: tasks use reference[0] as the
    # FASTA, and MergeMutect2Outputs takes reference[2], the dictionary.
    String ref_fasta_path = "~{ref_fasta}"
    Array[String] ref_gzi = if sub(ref_fasta_path, "\\.gz$", "") != ref_fasta_path then [ref_fasta_path + ".gzi"] else []
    Array[File] reference = flatten([[ref_fasta_path, ref_fasta_path + ".fai",
        sub(ref_fasta_path, "\\.(fasta|fas|fa|fna|txt)(\\.gz)?$", "") + ".dict"], ref_gzi])

    # Disk sizes used for dynamic sizing
    Int ref_size = ceil(size(reference, "GB"))
    # each sample's reads, including the index if the array has it
    Int tumor_reads_size = ceil(size(tumor_reads, "GB"))
    Int germline_resource_vcf_size = if defined(germline_resource_or_alias) then ceil(size(germline_resource_or_alias, "GB")) else 0
    Int normal_reads_size = ceil(size(normal_reads, "GB"))

    # This is added to every task as padding, should increase if systematically you need more disk for every call
    Int disk_pad = 10 + emergency_extra_disk

    Runtime standard_runtime = {"gatk_docker": gatk_docker, "gatk_override": gatk_override,
            "max_retries": max_retries, "preemptible": preemptible, "cpu": small_task_cpu,
            "machine_mem": small_task_mem * 1000, "command_mem": small_task_mem * 1000 - 500,
            "disk": small_task_disk + disk_pad, "boot_disk_size": boot_disk_size}

    Int m2_output_size = tumor_reads_size / scatter_count
    #TODO: do we need to change this disk size now that NIO is always going to happen (for the google backend only)
    Int m2_per_scatter_size = (tumor_reads_size + normal_reads_size) + ref_size + germline_resource_vcf_size + m2_output_size + disk_pad

    Boolean calculate_contamination = defined(common_variants_vcf_or_alias) && !skip_filtering && !skip_calculate_contamination
    Float pileups_max_af = select_first([pileups_max_population_af, contamination_max_population_af])

    call SplitIntervals {
        input:
            intervals = intervals,
            masked_intervals = masked_intervals,
            reference = reference,
            scatter_count = scatter_count,
            split_intervals_extra_args = split_intervals_extra_args,
            runtime_params = standard_runtime
    }

    scatter (subintervals in SplitIntervals.interval_files ) {
        call M2 {
            input:
                intervals = subintervals,
                reference = reference,
                tumor_reads = tumor_reads,
                normal_reads = normal_reads,
                pon = pon,
                pon_idx = pon_idx_or_conventional,
                germline_resource = germline_resource_or_alias,
                germline_resource_idx = germline_resource_idx_or_conventional,
                preemptible = preemptible,
                max_retries = max_retries,
                m2_extra_args = m2_extra_args,
                common_variants_vcf = common_variants_vcf_or_alias,
                common_variants_vcf_idx = common_variants_vcf_idx_or_conventional,
                pileups_min_population_af = pileups_min_population_af,
                pileups_max_population_af = pileups_max_af,
                dragstr_model = dragstr_model,
                make_bamout = make_bamout,
                run_ob_filter = run_orientation_bias_mixture_model_filter,
                compress_vcfs = compress_vcfs,
                gga_vcf = gga_vcf,
                gga_vcf_idx = gga_vcf_idx_or_conventional,
                make_permutect_training_dataset = make_permutect_training_dataset,
                make_permutect_test_dataset = make_permutect_test_dataset,
                permutect_training_dataset_truth_vcf = permutect_training_dataset_truth_vcf,
                permutect_training_dataset_truth_vcf_idx = permutect_training_dataset_truth_vcf_idx_or_conventional,
                permutect_test_dataset_truth_vcf = permutect_test_dataset_truth_vcf,
                permutect_test_dataset_truth_vcf_idx = permutect_test_dataset_truth_vcf_idx_or_conventional,
                gatk_override = gatk_override,
                gatk_docker = gatk_docker,
                disk_space = m2_per_scatter_size,
                gcs_project_for_requester_pays = gcs_project_for_requester_pays
        }
    }

    Int merged_vcf_size = ceil(size(M2.unfiltered_vcf, "GB"))
    Int merged_bamout_size = ceil(size(M2.output_bamOut, "GB"))

    if (run_orientation_bias_mixture_model_filter && (!skip_filtering)) {
        call LearnReadOrientationModel {
            input:
                f1r2_tar_gz = M2.f1r2_counts,
                runtime_params = standard_runtime,
                mem = learn_read_orientation_mem
        }
    }

    call MergeMutect2Outputs {
        input:
            input_vcfs = M2.unfiltered_vcf,
            input_vcf_indices = M2.unfiltered_vcf_idx,
            compress_vcfs = compress_vcfs,
            stats = M2.stats,
            tumor_pileups = flatten(M2.tumor_pileups),
            normal_pileups = flatten(M2.normal_pileups),
            ref_dict = reference[2],
            runtime_params = standard_runtime
    }

    if (make_bamout) {
        call MergeBamOuts {
            input:
                reference = reference,
                bam_outs = M2.output_bamOut,
                runtime_params = standard_runtime,
                disk_space = ceil(merged_bamout_size * 10) + disk_pad,
        }
    }

    # Mutect2 tabulates pileup summaries at the common variants, which MergeMutect2Outputs merges for the outputs and for
    # CalculateContamination, whose contamination and minor allele fraction segments FilterMutectCalls uses
    if (calculate_contamination) {
        call CalculateContamination {
            input:
                tumor_pileups = select_first([MergeMutect2Outputs.merged_tumor_pileups]),
                normal_pileups = MergeMutect2Outputs.merged_normal_pileups,
                min_population_af = contamination_min_population_af,
                max_population_af = contamination_max_population_af,
                runtime_params = standard_runtime
        }
    }

    if (make_permutect_training_dataset) {
        call Concatenate as ConcatenatePermutectTrainingData {
            input:
                input_files = M2.permutect_training_dataset,
                gatk_docker = gatk_docker
        }
    }

    if (make_permutect_test_dataset) {
        call Concatenate as ConcatenatePermutectTestData {
            input:
                input_files = M2.permutect_test_dataset,
                gatk_docker = gatk_docker
        }
    }

    if (!skip_filtering) {
        call Filter {
            input:
                reference = reference,
                intervals = intervals,
                unfiltered_vcf = MergeMutect2Outputs.merged_vcf,
                unfiltered_vcf_idx = MergeMutect2Outputs.merged_vcf_idx,
                compress_vcfs = compress_vcfs,
                mutect_stats = MergeMutect2Outputs.merged_stats,
                contamination_table = CalculateContamination.contamination_table,
                maf_segments = CalculateContamination.maf_segments,
                artifact_priors_tar_gz = LearnReadOrientationModel.artifact_prior_table,
                m2_extra_filtering_args = m2_extra_filtering_args,
                runtime_params = standard_runtime,
                disk_space = ceil(size(MergeMutect2Outputs.merged_vcf, "GB") * 4) + disk_pad
        }

        if (defined(realignment_index_bundle)) {
            call FilterAlignmentArtifacts {
                input:
                    reference = reference,
                    reads = tumor_reads,
                    realignment_index_bundle = select_first([realignment_index_bundle]),
                    realignment_extra_args = realignment_extra_args,
                    compress_vcfs = compress_vcfs,
                    input_vcf = Filter.filtered_vcf,
                    input_vcf_idx = Filter.filtered_vcf_idx,
                    runtime_params = standard_runtime,
                    mem = filter_alignment_artifacts_mem,
                    gcs_project_for_requester_pays = gcs_project_for_requester_pays
            }
        }
    }

    output {
        File output_vcf = select_first([FilterAlignmentArtifacts.filtered_vcf, Filter.filtered_vcf, MergeMutect2Outputs.merged_vcf])
        File output_vcf_idx = select_first([FilterAlignmentArtifacts.filtered_vcf_idx, Filter.filtered_vcf_idx, MergeMutect2Outputs.merged_vcf_idx])
        File? filtering_stats = Filter.filtering_stats
        File mutect_stats = MergeMutect2Outputs.merged_stats
        File? contamination_table = CalculateContamination.contamination_table

        File? bamout = MergeBamOuts.merged_bam_out
        File? bamout_index = MergeBamOuts.merged_bam_out_index
        File? maf_segments = CalculateContamination.maf_segments
        File? read_orientation_model_params = LearnReadOrientationModel.artifact_prior_table
        File? permutect_training_dataset = ConcatenatePermutectTrainingData.concatenated
        File? permutect_test_dataset = ConcatenatePermutectTestData.concatenated
        File permutect_contigs_table = select_first(M2.permutect_contigs_table)
        File permutect_read_groups_table = select_first(M2.permutect_read_groups_table)
        File? tumor_pileups = MergeMutect2Outputs.merged_tumor_pileups
        File? normal_pileups = MergeMutect2Outputs.merged_normal_pileups
        # without the orientation bias filter, M2's F1R2 counts are empty placeholder files
        Array[File] f1r2_counts = if run_orientation_bias_mixture_model_filter then M2.f1r2_counts else []
    }
}

task SplitIntervals {
    input {
      File? intervals
      File? masked_intervals
      # the reference FASTA, then its .fai index, its .dict dictionary and, if it's bgzipped, its .gzi block index
      Array[File] reference
      Int scatter_count
      String? split_intervals_extra_args

      # runtime
      Runtime runtime_params
    }

    command {
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}

        mkdir interval-files
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" SplitIntervals \
            -R ~{reference[0]} \
            ~{"-L " + intervals} \
            ~{"-XL " + masked_intervals} \
            -scatter ~{scatter_count} \
            -O interval-files \
            ~{split_intervals_extra_args}
        cp interval-files/*.interval_list .
    }

    runtime {
        docker: runtime_params.gatk_docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.machine_mem + " MB"
        disks: "local-disk " + runtime_params.disk + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    output {
        Array[File] interval_files = glob("*.interval_list")
    }
}

task M2 {
    input {
        File? intervals
        # the reference FASTA, then its .fai index, its .dict dictionary and, if it's bgzipped, its .gzi block index
        Array[File] reference
        # Each sample's reads: a BAM or CRAM, optionally followed by its index, which GATK finds by name beside the reads
        # (see tumor_reads in the workflow's header).  The command uses tumor_reads[0] and normal_reads[0] as the reads,
        # and there is no normal if normal_reads is empty.
        Array[File]+ tumor_reads
        Array[File] normal_reads
        File? pon
        File? pon_idx
        File? germline_resource
        File? germline_resource_idx
        String? m2_extra_args
        Boolean? make_bamout
        Boolean? run_ob_filter
        Boolean compress_vcfs
        File? gga_vcf
        File? gga_vcf_idx
        File? common_variants_vcf
        File? common_variants_vcf_idx
        Float pileups_min_population_af = 0.01
        Float pileups_max_population_af = 0.2
        File? dragstr_model

        File? gatk_override

        String? gcs_project_for_requester_pays

        Boolean make_permutect_training_dataset = false
        Boolean make_permutect_test_dataset = false
        File? permutect_training_dataset_truth_vcf
        File? permutect_training_dataset_truth_vcf_idx
        File? permutect_test_dataset_truth_vcf
        File? permutect_test_dataset_truth_vcf_idx

        # runtime
        String gatk_docker
        Int? mem
        Int? preemptible
        Int? max_retries
        Int? disk_space
        Int? cpu
        Boolean use_ssd = false
    }

    String output_vcf = "output" + if compress_vcfs then ".vcf.gz" else ".vcf"
    String output_vcf_idx = output_vcf + if compress_vcfs then ".tbi" else ".idx"

    String output_stats = output_vcf + ".stats"

    Boolean has_normal = length(normal_reads) > 0

    # Mutect2 tabulates the pileup summaries of the tumor and normal at the common variants as it traverses the reads
    String pileup_summaries_args = if defined(common_variants_vcf) then
        "--tumor-pileup-summaries tumor-pileups.table" +
        (if has_normal then " --normal-pileup-summaries normal-pileups.table" else "") +
        " --pileup-summaries-min-population-af ~{pileups_min_population_af} --pileup-summaries-max-population-af ~{pileups_max_population_af}"
        else ""

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1000 else 3500
    Int command_mem = machine_mem - 500

    parameter_meta{
        intervals: {localization_optional: true}
        reference: {localization_optional: true}
        tumor_reads: {localization_optional: true}
        normal_reads: {localization_optional: true}
        pon: {localization_optional: true}
        pon_idx: {localization_optional: true}
        germline_resource: {localization_optional: true}
        germline_resource_idx: {localization_optional: true}
        gga_vcf: {localization_optional: true}
        gga_vcf_idx: {localization_optional: true}
        common_variants_vcf: {localization_optional: true}
        common_variants_vcf_idx: {localization_optional: true}
        permutect_training_dataset_truth_vcf: {localization_optional: true}
        permutect_training_dataset_truth_vcf_idx: {localization_optional: true}
    }

    command <<<
        set -e

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

        # We need to create these files regardless, even if they stay empty
        touch bamout.bam
        touch f1r2.tar.gz
        touch training-dataset.txt
        touch test-dataset.txt
        touch contigs.table
        touch read-groups.table

        # The normal's BAM or CRAM, or "" without a normal.  Cromwell evaluates every placeholder, even in a part of the
        # script that doesn't run, so the condition keeps it from indexing an empty array.  (A declaration such as
        # String normal = normal_reads[0] wouldn't do, because it would hold the path from before localization.)
        normal_bam_or_cram="~{if has_normal then normal_reads[0] else ""}"

        if [[ -n "$normal_bam_or_cram" ]]; then
            gatk --java-options "-Xmx~{command_mem}m" GetSampleName -R ~{reference[0]} -I $normal_bam_or_cram -O normal_names.txt -encode \
            ~{"--gcs-project-for-requester-pays " + gcs_project_for_requester_pays}
            # add "-normal " to the start of each line and " " to the end, then remove newlines
            # to get -normal sample1 -normal sample2 etc
            normal_sample_line=`awk '{ print "-normal", $0 }' normal_names.txt | tr '\n' ' '`
            normal_input="-I $normal_bam_or_cram"
        fi

        gatk --java-options "-Xmx~{command_mem}m" Mutect2 \
            -R ~{reference[0]} \
            -I ~{tumor_reads[0]} \
            $normal_input \
            $normal_sample_line \
            ~{"--germline-resource " + germline_resource} \
            ~{"-pon " + pon} \
            ~{"-L " + intervals} \
            ~{"--alleles " + gga_vcf} \
            ~{"--dragstr-params-path " + dragstr_model} \
            -O "~{output_vcf}" \
            ~{true='--bam-output bamout.bam' false='' make_bamout} \
            ~{true='--f1r2-tar-gz f1r2.tar.gz' false='' run_ob_filter} \
            ~{"--common-variants-for-pileup-summaries " + common_variants_vcf} \
            ~{pileup_summaries_args} \
            ~{true='--permutect-training-dataset training-dataset.txt' false='' make_permutect_training_dataset} \
            ~{true='--permutect-test-dataset test-dataset.txt' false='' make_permutect_test_dataset} \
            ~{"--permutect-training-truth " + permutect_training_dataset_truth_vcf} \
            ~{"--permutect-test-truth " + permutect_test_dataset_truth_vcf} \
            ~{m2_extra_args} \
            ~{"--gcs-project-for-requester-pays " + gcs_project_for_requester_pays}
    >>>

    runtime {
        docker: gatk_docker
        bootDiskSizeGb: 12
        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, 100]) + if use_ssd then " SSD" else " HDD"
        preemptible: select_first([preemptible, 10])
        maxRetries: select_first([max_retries, 0])
        cpu: select_first([cpu, 1])
    }

    output {
        File unfiltered_vcf = "~{output_vcf}"
        File unfiltered_vcf_idx = "~{output_vcf_idx}"
        File output_bamOut = "bamout.bam"
        File stats = "~{output_stats}"
        File f1r2_counts = "f1r2.tar.gz"
        Array[File] tumor_pileups = glob("*tumor-pileups.table")
        Array[File] normal_pileups = glob("*normal-pileups.table")
        File permutect_training_dataset = "training-dataset.txt"
        File permutect_test_dataset = "test-dataset.txt"
        File permutect_contigs_table = "contigs.table"
        File permutect_read_groups_table = "read-groups.table"
    }
}

# Merges the scattered outputs of M2 other than the bamouts, which need a bigger disk: the VCFs, the stats and, when
# there are common variants, the tumor's and normal's pileup summaries.  A merged table of pileups exists only if
# there are tables to merge.
task MergeMutect2Outputs {
    input {
      Array[File] input_vcfs
      Array[File] input_vcf_indices
      Boolean compress_vcfs
      Array[File]+ stats
      Array[File] tumor_pileups = []
      Array[File] normal_pileups = []
      File ref_dict
      Runtime runtime_params
    }

    String output_vcf = if compress_vcfs then "merged.vcf.gz" else "merged.vcf"
    String output_vcf_idx = output_vcf + if compress_vcfs then ".tbi" else ".idx"

    # GatherPileupSummaries streams the pileups from the bucket on Terra
    parameter_meta {
        tumor_pileups: {localization_optional: true}
        normal_pileups: {localization_optional: true}
    }

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}

        # using MergeVcfs instead of GatherVcfs so we can create indices
        # WARNING 2015-10-28 15:01:48 GatherVcfs  Index creation not currently supported when gathering block compressed VCFs.
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" MergeVcfs -I ~{sep=' -I ' input_vcfs} -O ~{output_vcf}

        gatk --java-options "-Xmx~{runtime_params.command_mem}m" MergeMutectStats \
            -stats ~{sep=" -stats " stats} -O merged.stats

        if [ ~{length(tumor_pileups)} -gt 0 ]; then
            gatk --java-options "-Xmx~{runtime_params.command_mem}m" GatherPileupSummaries \
                --sequence-dictionary ~{ref_dict} -I ~{sep=' -I ' tumor_pileups} -O tumor-pileups.tsv
        fi

        if [ ~{length(normal_pileups)} -gt 0 ]; then
            gatk --java-options "-Xmx~{runtime_params.command_mem}m" GatherPileupSummaries \
                --sequence-dictionary ~{ref_dict} -I ~{sep=' -I ' normal_pileups} -O normal-pileups.tsv
        fi
    >>>

    runtime {
        docker: runtime_params.gatk_docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.machine_mem + " MB"
        disks: "local-disk " + runtime_params.disk + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    output {
        File merged_vcf = "~{output_vcf}"
        File merged_vcf_idx = "~{output_vcf_idx}"
        File merged_stats = "merged.stats"
        File? merged_tumor_pileups = "tumor-pileups.tsv"
        File? merged_normal_pileups = "normal-pileups.tsv"
    }
}

task MergeBamOuts {
    input {
      # the reference FASTA, then its .fai index, its .dict dictionary and, if it's bgzipped, its .gzi block index
      Array[File] reference
      Array[File]+ bam_outs
      Runtime runtime_params
      Int? disk_space   #override to request more disk than default small task params
    }

    command <<<
        # This command block assumes that there is at least one file in bam_outs.
        #  Do not call this task if len(bam_outs) == 0
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" GatherBamFiles \
            -I ~{sep=" -I " bam_outs} -O unsorted.out.bam -R ~{reference[0]}

        # We must sort because adjacent scatters may have overlapping (padded) assembly regions, hence
        # overlapping bamouts

        gatk --java-options "-Xmx~{runtime_params.command_mem}m" SortSam -I unsorted.out.bam \
            -O bamout.bam --SORT_ORDER coordinate -VALIDATION_STRINGENCY LENIENT
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" BuildBamIndex -I bamout.bam -VALIDATION_STRINGENCY LENIENT
    >>>

    runtime {
        docker: runtime_params.gatk_docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, runtime_params.disk]) + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    output {
        File merged_bam_out = "bamout.bam"
        File merged_bam_out_index = "bamout.bai"
    }
}

# Learning step of the orientation bias mixture model, which is the recommended orientation bias filter as of September 2018
task LearnReadOrientationModel {
    input {
      Array[File] f1r2_tar_gz
      Runtime runtime_params
      Int? mem  #override memory
    }

    Int machine_mem = select_first([mem, runtime_params.machine_mem])
    Int command_mem = machine_mem - 1000

    command {
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}

        gatk --java-options "-Xmx~{command_mem}m" LearnReadOrientationModel \
            -I ~{sep=" -I " f1r2_tar_gz} \
            -O "artifact-priors.tar.gz"
    }

    runtime {
        docker: runtime_params.gatk_docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: machine_mem + " MB"
        disks: "local-disk " + runtime_params.disk + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    output {
        File artifact_prior_table = "artifact-priors.tar.gz"
    }

}

task CalculateContamination {
    input {
      String? intervals
      File tumor_pileups
      File? normal_pileups
      Float min_population_af = 0.01
      Float max_population_af = 0.2
      Runtime runtime_params
    }

    command {
        set -e

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}

        gatk --java-options "-Xmx~{runtime_params.command_mem}m" CalculateContamination -I ~{tumor_pileups} \
        -O contamination.table --tumor-segmentation segments.table ~{"-matched " + normal_pileups} \
        --minimum-population-allele-frequency ~{min_population_af} --maximum-population-allele-frequency ~{max_population_af}
    }

    runtime {
        docker: runtime_params.gatk_docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.machine_mem + " MB"
        disks: "local-disk " + runtime_params.disk + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    output {
        File contamination_table = "contamination.table"
        File maf_segments = "segments.table"
    }
}

task Filter {
    input {
        File? intervals
        # the reference FASTA, then its .fai index, its .dict dictionary and, if it's bgzipped, its .gzi block index
        Array[File] reference
        File unfiltered_vcf
        File unfiltered_vcf_idx
        Boolean compress_vcfs
        File? mutect_stats
        File? artifact_priors_tar_gz
        File? contamination_table
        File? maf_segments
        String? m2_extra_filtering_args

        Runtime runtime_params
        Int? disk_space
    }

    String output_vcf = if compress_vcfs then "filtered.vcf.gz" else "filtered.vcf"
    String output_vcf_idx = output_vcf + if compress_vcfs then ".tbi" else ".idx"

    parameter_meta{
      reference: {localization_optional: true}
    }

    command {
        set -e

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}

        gatk --java-options "-Xmx~{runtime_params.command_mem}m" FilterMutectCalls -V ~{unfiltered_vcf} \
            -R ~{reference[0]} \
            -O ~{output_vcf} \
            ~{"--contamination-table " + contamination_table} \
            ~{"--tumor-segmentation " + maf_segments} \
            ~{"--ob-priors " + artifact_priors_tar_gz} \
            ~{"-stats " + mutect_stats} \
            --filtering-stats filtering.stats \
            ~{m2_extra_filtering_args}
    }

    runtime {
        docker: runtime_params.gatk_docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, runtime_params.disk]) + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    output {
        File filtered_vcf = "~{output_vcf}"
        File filtered_vcf_idx = "~{output_vcf_idx}"
        File filtering_stats = "filtering.stats"
    }
}

task FilterAlignmentArtifacts {
    input {
      # the reference FASTA, then its .fai index, its .dict dictionary and, if it's bgzipped, its .gzi block index
      Array[File] reference
      File input_vcf
      File input_vcf_idx
      # the reads: a BAM or CRAM, which the command uses as reads[0], optionally followed by its index, which GATK finds
      # by name beside the reads (see tumor_reads in the workflow's header)
      Array[File]+ reads
      Boolean compress_vcfs
      File realignment_index_bundle
      String? realignment_extra_args
      String? gcs_project_for_requester_pays
      Runtime runtime_params
      Int mem
    }

    String output_vcf = if compress_vcfs then "filtered.vcf.gz" else "filtered.vcf"
    String output_vcf_idx = output_vcf +  if compress_vcfs then ".tbi" else ".idx"

    Int machine_mem = mem
    Int command_mem = machine_mem - 500

    parameter_meta{
      reference: {localization_optional: true}
      input_vcf: {localization_optional: true}
      input_vcf_idx: {localization_optional: true}
      reads: {localization_optional: true}
    }

    command {
        set -e

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}

        gatk --java-options "-Xmx~{command_mem}m" FilterAlignmentArtifacts \
            -R ~{reference[0]} \
            -V ~{input_vcf} \
            -I ~{reads[0]} \
            --bwa-mem-index-image ~{realignment_index_bundle} \
            ~{realignment_extra_args} \
            -O ~{output_vcf} \
            ~{"--gcs-project-for-requester-pays " + gcs_project_for_requester_pays}
    }

    runtime {
        docker: runtime_params.gatk_docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: machine_mem + " MB"
        disks: "local-disk " + runtime_params.disk + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    output {
        File filtered_vcf = "~{output_vcf}"
        File filtered_vcf_idx = "~{output_vcf_idx}"
    }
}

task Concatenate {
    input {
        Array[File] input_files
        Int? mem
        String gatk_docker
    }

    Int machine_mem = if defined(mem) then mem * 1000 else 7000

    command {
        cat ~{sep=' ' input_files} > output.txt
    }

    runtime {
        docker: gatk_docker
        bootDiskSizeGb: 12
        memory: machine_mem + " MB"
        disks: "local-disk 100 HDD"
        preemptible: 1
        maxRetries: 1
        cpu: 2
    }

    output {
        File concatenated = "output.txt"
    }
}

