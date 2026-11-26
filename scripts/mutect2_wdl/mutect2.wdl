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
## <task>.preemptible, eg M2.preemptible: each task's own number of preemptible VMs to try on Google before a regular
##     one (see above the Runtime struct for the defaults)
## max_retries: how many times to retry failed tasks -- very important on the cloud when there are transient errors
## <task>.mem, <task>.disk_space, eg Filter.mem: each task sizes its memory and disk from its inputs (see above the
##     Runtime struct), and these override that: memory in MB (in GB for M2.mem), disk in GB
## emergency_extra_disk: GB added to every task's disk, as a last resort for ill-behaved data
## gatk_override: (optional) local file or Google bucket path to a GATK 4 java jar file to be used instead of the GATK 4 jar
##                in the docker image.  This must be supplied when running in an environment that does not support docker
##                (e.g. SGE cluster on a Broad on-prem VM)
##
## ** Workflow options **
## intervals: files of the genomic intervals to call, such as .interval_list or .bed files, or [] to call the whole
##     genome.  Tasks give GATK each file with its own -L, and GATK takes the union of their intervals, which the
##     workflow splits into scatter_count parts.  In a Terra configuration, an attribute expression such as
##     workspace.intervals gives a one-element array.  (Earlier versions took a single, optional file.)
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
## ** Primary resources ** (optional but strongly recommended, so each must be given deliberately, if only as [] to go
## without it)
## pon: panel of normals in VCF format containing probable technical artifacts (false positives)
## germline_resource: database of known germline variants with population allele frequencies (see
##     http://gnomad.broadinstitute.org/downloads)
## common_variants_vcf: VCF of common variants with allele frequencies, at whose sites Mutect2 tabulates pileup summaries,
##                      for calculating contamination
##
## ** Other VCFs ** (genuinely optional, so they may be left unset)
## gga_vcf: alleles to force-call
## permutect_training_dataset_truth_vcf, permutect_test_dataset_truth_vcf: truth VCFs that label the Permutect datasets
##
## ** VCFs and their indices **
## Each VCF input above is an array of the VCF and, optionally, its index: [vcf] or [vcf, idx], or [] for none, eg
##     ["gs://bucket/af-only-gnomad.hg38.vcf.gz"].  GATK finds a VCF's index only by name beside the VCF: the VCF's path
##     plus .tbi for a bgzipped VCF (.vcf.gz or .vcf.bgz), or plus .idx for a plain .vcf, such as
##     af-only-gnomad.hg38.vcf.gz.tbi.  Given [vcf], the workflow derives that path and passes the index to tasks along
##     with the VCF, so that backends that localize files put it beside the VCF, while on Terra the tasks stream both.
##     So [vcf] suffices on every backend, and an index given as the second element must be that same file.  In a Terra
##     configuration, an attribute expression such as workspace.gnomad gives a one-element array, and [workspace.gnomad,
##     workspace.gnomad_idx] a two-element one.  (Earlier versions took each index as a separate input, such as
##     germline_resource_idx, and also accepted the former names gnomad and variants_for_contamination of
##     germline_resource and common_variants_vcf.)
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


# Every task but M2 runs on 1 CPU, since its tools are single-threaded.  Mutect2's native PairHMM can use more threads,
# but a larger scatter_count is a cheaper way to speed up M2 than more CPUs per scatter.
# Every task requests 4 GB of boot disk, beyond the 30 GB that Cromwell on GCP Batch adds for the OS and the docker image
# (Terra records 42 GB for a request of 12).  That leaves room for the GATK and Permutect images, 2.4 and 3.9 GB
# compressed, and tasks write their files to the working disk, not the boot disk.
# Each task's preemptible input is how many preemptible (spot) VMs Google tries before a regular one: 3 for tasks that
# take minutes, 2 for tasks that can take an hour or two, and 1 for M2, whose scatters can take hours on genomes and
# hold up the workflow when they're preempted.
# Each task sizes its memory and its disk from the sizes of its inputs, by rules fit to a benchmark on Terra in November
# 2026 of exomes, genomes and FFPE samples, in which each JVM's heap after garbage collection and each disk's peak use
# were logged.  The JVMs get the machine's memory less 500 MB, and a disk is at least 10 GB.  Every disk is an SSD: the
# disks are small, an SSD localizes files several times faster than a standard disk of any size (half a minute against
# 1.6 minutes for the 3.2 GB hg38 reference), and the price of a small SSD is a few percent of the VM's.  Backends other
# than Google's ignore disk sizes and types.
struct Runtime {
    String gatk_docker
    File? gatk_override
    Int max_retries
    # GB added to every disk
    Int emergency_extra_disk
}

workflow Mutect2 {
    input {
        # basic inputs
        # files of the genomic intervals to call, whose union the workflow calls, or [] for the whole genome
        Array[File] intervals
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

        # VCFs, each an array of the VCF and, optionally, its index: [vcf] or [vcf, idx], or [] for none (see "VCFs and
        # their indices" in the header).  These resources are optional but strongly recommended, so each must be given,
        # if only as [] to go without it.
        Array[File] germline_resource
        Array[File] pon
        Array[File] common_variants_vcf
        # These are genuinely optional and may be left unset: alleles to force-call, and truth VCFs for Permutect
        # datasets
        Array[File]? gga_vcf
        Array[File]? permutect_training_dataset_truth_vcf
        Array[File]? permutect_test_dataset_truth_vcf

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
        # each task has its own preemptible input, and sizes its memory and disk (see above the Runtime struct)
        Int max_retries = 1
        String? gcs_project_for_requester_pays

        # Use as a last resort to increase the disk given to every task in case of ill behaving data
        Int emergency_extra_disk = 0
    }

    # Each VCF with its index: [] without the VCF, the array as given if it has the index, and otherwise the VCF and the
    # index where GATK looks for it, the VCF's path plus .tbi for a bgzipped VCF (.vcf.gz or .vcf.bgz) or .idx for a
    # plain .vcf.  Tasks take both, so that a backend that localizes files puts the index beside the VCF, while on Terra
    # both are streamed.  WDL has no user-defined functions, so we scatter this rule over the VCFs; a scatter without
    # calls only evaluates declarations.  Keep given_vcfs and the indices into vcf_with_index below in the same order.
    Array[Array[File]] given_vcfs = [germline_resource, pon, common_variants_vcf, select_first([gga_vcf, []]),
        select_first([permutect_training_dataset_truth_vcf, []]), select_first([permutect_test_dataset_truth_vcf, []])]
    scatter (given_vcf in given_vcfs) {
        String vcf_path = if length(given_vcf) > 0 then "~{given_vcf[0]}" else ""
        Array[File] vcf_with_index = if length(given_vcf) == 1
            then [vcf_path, vcf_path + (if sub(vcf_path, "\\.b?gz$", "") != vcf_path then ".tbi" else ".idx")]
            else given_vcf
    }
    Array[File] germline_resource_with_index = vcf_with_index[0]
    Array[File] pon_with_index = vcf_with_index[1]
    Array[File] common_variants_vcf_with_index = vcf_with_index[2]
    Array[File] gga_vcf_with_index = vcf_with_index[3]
    Array[File] permutect_training_dataset_truth_vcf_with_index = vcf_with_index[4]
    Array[File] permutect_test_dataset_truth_vcf_with_index = vcf_with_index[5]

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

    Runtime standard_runtime = {"gatk_docker": gatk_docker, "gatk_override": gatk_override,
            "max_retries": max_retries, "emergency_extra_disk": emergency_extra_disk}

    # M2 localizes the reference and streams its other inputs, so its disk holds the reference and its outputs: at most
    # 0.25 GB on Terra in November 2026, but for a bamout, which this bounds by the tumor reads in a scatter's share of
    # the genome.  The workflow sizes M2's disk because the task doesn't know the scatter count.
    Int m2_output_size = if make_bamout then ceil(size(tumor_reads, "GB")) / scatter_count + 1 else 1
    Int m2_disk_size = ceil(size(reference, "GB")) + m2_output_size + 5 + emergency_extra_disk

    Boolean calculate_contamination = length(common_variants_vcf) > 0 && !skip_filtering
        && !skip_calculate_contamination
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
                pon = pon_with_index,
                germline_resource = germline_resource_with_index,
                max_retries = max_retries,
                m2_extra_args = m2_extra_args,
                common_variants_vcf = common_variants_vcf_with_index,
                pileups_min_population_af = pileups_min_population_af,
                pileups_max_population_af = pileups_max_af,
                dragstr_model = dragstr_model,
                make_bamout = make_bamout,
                run_ob_filter = run_orientation_bias_mixture_model_filter,
                compress_vcfs = compress_vcfs,
                gga_vcf = gga_vcf_with_index,
                make_permutect_training_dataset = make_permutect_training_dataset,
                make_permutect_test_dataset = make_permutect_test_dataset,
                permutect_training_dataset_truth_vcf = permutect_training_dataset_truth_vcf_with_index,
                permutect_test_dataset_truth_vcf = permutect_test_dataset_truth_vcf_with_index,
                gatk_override = gatk_override,
                gatk_docker = gatk_docker,
                disk_space = m2_disk_size,
                gcs_project_for_requester_pays = gcs_project_for_requester_pays
        }
    }

    if (run_orientation_bias_mixture_model_filter && (!skip_filtering)) {
        call LearnReadOrientationModel {
            input:
                f1r2_tar_gz = M2.f1r2_counts,
                runtime_params = standard_runtime
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
                bam_outs = M2.output_bamOut,
                runtime_params = standard_runtime
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
                runtime_params = standard_runtime
        }
    }

    if (make_permutect_test_dataset) {
        call Concatenate as ConcatenatePermutectTestData {
            input:
                input_files = M2.permutect_test_dataset,
                runtime_params = standard_runtime
        }
    }

    if (!skip_filtering) {
        call Filter {
            input:
                reference = reference,
                unfiltered_vcf = MergeMutect2Outputs.merged_vcf,
                unfiltered_vcf_idx = MergeMutect2Outputs.merged_vcf_idx,
                compress_vcfs = compress_vcfs,
                mutect_stats = MergeMutect2Outputs.merged_stats,
                contamination_table = CalculateContamination.contamination_table,
                maf_segments = CalculateContamination.maf_segments,
                artifact_priors_tar_gz = LearnReadOrientationModel.artifact_prior_table,
                m2_extra_filtering_args = m2_extra_filtering_args,
                runtime_params = standard_runtime
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
      # files of genomic intervals, whose union the task splits, or none for the whole genome
      Array[File] intervals
      File? masked_intervals
      # the reference FASTA, then its .fai index, its .dict dictionary and, if it's bgzipped, its .gzi block index
      Array[File] reference
      Int scatter_count
      String? split_intervals_extra_args

      # runtime
      Runtime runtime_params
      Int preemptible = 3
      Int? mem
      Int? disk_space
    }

    # memory in MB and disk in GB, from the interval files (a heap of 320 MB for 17 MB of hg38 intervals)
    Float intervals_mb = size(intervals, "MB") + size(masked_intervals, "MB")
    Int machine_mem = select_first([mem, 1500 + ceil(30 * intervals_mb)])
    Int disk_needed = ceil(3 * intervals_mb / 1000) + 5 + runtime_params.emergency_extra_disk
    Int disk = select_first([disk_space, if disk_needed > 10 then disk_needed else 10])

    # SplitIntervals needs only the reference's dictionary and index, so on Terra GATK streams them rather than the
    # whole FASTA being copied to the disk
    parameter_meta {
      reference: {localization_optional: true}
    }

    command {
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}

        mkdir interval-files
        gatk --java-options "-Xmx~{machine_mem - 500}m" SplitIntervals \
            -R ~{reference[0]} \
            ~{if length(intervals) > 0 then "-L " else ""}~{sep=" -L " intervals} \
            ~{"-XL " + masked_intervals} \
            -scatter ~{scatter_count} \
            -O interval-files \
            ~{split_intervals_extra_args}
        cp interval-files/*.interval_list .
    }

    runtime {
        docker: runtime_params.gatk_docker
        bootDiskSizeGb: 4
        memory: machine_mem + " MB"
        disks: "local-disk " + disk + " SSD"
        preemptible: preemptible
        maxRetries: runtime_params.max_retries
        cpu: 1
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
        # Each VCF with its index, or [] for none (see "VCFs and their indices" in the workflow's header).  The command
        # uses element 0 as the VCF, and GATK finds the index by name beside it.
        Array[File] pon
        Array[File] germline_resource
        String? m2_extra_args
        Boolean? make_bamout
        Boolean? run_ob_filter
        Boolean compress_vcfs
        Array[File] gga_vcf
        Array[File] common_variants_vcf
        Float pileups_min_population_af = 0.01
        Float pileups_max_population_af = 0.2
        File? dragstr_model

        File? gatk_override

        String? gcs_project_for_requester_pays

        Boolean make_permutect_training_dataset = false
        Boolean make_permutect_test_dataset = false
        Array[File] permutect_training_dataset_truth_vcf
        Array[File] permutect_test_dataset_truth_vcf

        # runtime
        String gatk_docker
        Int? mem
        Int preemptible = 1
        Int? max_retries
        Int? disk_space
        Int? cpu
        # The CPU platform Google gives the scatters.  On Terra in November 2026, Intel Ice Lake ran scatters in 0.51
        # to 0.66 times the time on the default CPUs (2.3 GHz Haswell) for 0.97 to 1.03 times the cost: Cromwell gives
        # it an N2 machine, with at least 2 CPUs, both of which the native PairHMM uses, and its AVX-512 instructions.
        # "Intel Skylake" gives a 1-CPU N1 machine with AVX-512 at the default's price, for 0.8 to 0.9 times the time.
        String cpu_platform = "Intel Ice Lake"
    }

    String output_vcf = "output" + if compress_vcfs then ".vcf.gz" else ".vcf"
    String output_vcf_idx = output_vcf + if compress_vcfs then ".tbi" else ".idx"

    String output_stats = output_vcf + ".stats"

    Boolean has_normal = length(normal_reads) > 0

    # Mutect2 tabulates the pileup summaries of the tumor and normal at the common variants as it traverses the reads
    String pileup_summaries_args = if length(common_variants_vcf) > 0 then
        "--tumor-pileup-summaries tumor-pileups.table" +
        (if has_normal then " --normal-pileup-summaries normal-pileups.table" else "") +
        " --pileup-summaries-min-population-af ~{pileups_min_population_af} --pileup-summaries-max-population-af ~{pileups_max_population_af}"
        else ""

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1000 else 3500
    Int command_mem = machine_mem - 500

    # The reference is localized rather than streamed.  GATK reads a streamed FASTA one 128 kB request at a time, and
    # on Terra in November 2026 that cost each scatter about 4 to 7 minutes of waiting, while localizing the reference
    # to an SSD took half a minute.  Once GATK streams FASTAs efficiently, this may be reversed: mark the reference
    # localization_optional again, leave it out of the disk size, and drop the SSD.
    parameter_meta{
        intervals: {localization_optional: true}
        tumor_reads: {localization_optional: true}
        normal_reads: {localization_optional: true}
        pon: {localization_optional: true}
        germline_resource: {localization_optional: true}
        gga_vcf: {localization_optional: true}
        common_variants_vcf: {localization_optional: true}
        permutect_training_dataset_truth_vcf: {localization_optional: true}
        permutect_test_dataset_truth_vcf: {localization_optional: true}
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

        # Each VCF array is empty or holds the VCF and its index, so an argument indexes it only if it isn't empty
        gatk --java-options "-Xmx~{command_mem}m" Mutect2 \
            -R ~{reference[0]} \
            -I ~{tumor_reads[0]} \
            $normal_input \
            $normal_sample_line \
            ~{if length(germline_resource) > 0 then "--germline-resource " + germline_resource[0] else ""} \
            ~{if length(pon) > 0 then "-pon " + pon[0] else ""} \
            ~{"-L " + intervals} \
            ~{if length(gga_vcf) > 0 then "--alleles " + gga_vcf[0] else ""} \
            ~{"--dragstr-params-path " + dragstr_model} \
            -O "~{output_vcf}" \
            ~{true='--bam-output bamout.bam' false='' make_bamout} \
            ~{true='--f1r2-tar-gz f1r2.tar.gz' false='' run_ob_filter} \
            ~{if length(common_variants_vcf) > 0 then "--common-variants-for-pileup-summaries " + common_variants_vcf[0] else ""} \
            ~{pileup_summaries_args} \
            ~{true='--permutect-training-dataset training-dataset.txt' false='' make_permutect_training_dataset} \
            ~{true='--permutect-test-dataset test-dataset.txt' false='' make_permutect_test_dataset} \
            ~{if length(permutect_training_dataset_truth_vcf) > 0 then "--permutect-training-truth " + permutect_training_dataset_truth_vcf[0] else ""} \
            ~{if length(permutect_test_dataset_truth_vcf) > 0 then "--permutect-test-truth " + permutect_test_dataset_truth_vcf[0] else ""} \
            ~{m2_extra_args} \
            ~{"--gcs-project-for-requester-pays " + gcs_project_for_requester_pays}
    >>>

    runtime {
        docker: gatk_docker
        bootDiskSizeGb: 4
        memory: machine_mem + " MB"
        # an SSD, which localizes the reference (3.2 GB for hg38) in half a minute, against more than a minute and a
        # half for a standard disk of any size
        disks: "local-disk " + select_first([disk_space, 100]) + " SSD"
        preemptible: preemptible
        maxRetries: select_first([max_retries, 0])
        cpu: select_first([cpu, 1])
        cpuPlatform: cpu_platform
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
      Int preemptible = 3
      Int? mem
      Int? disk_space
    }

    String output_vcf = if compress_vcfs then "merged.vcf.gz" else "merged.vcf"
    String output_vcf_idx = output_vcf + if compress_vcfs then ".tbi" else ".idx"

    # memory in MB and disk in GB, from the VCFs and pileups (for a genome's 0.3 to 0.4 GB of VCFs and 0.54 GB of
    # pileups, a heap of 1 GB and 1.3 GB of disk)
    Float vcfs_gb = size(input_vcfs, "GB") + size(input_vcf_indices, "GB")
    Float pileups_gb = size(tumor_pileups, "GB") + size(normal_pileups, "GB")
    Int machine_mem = select_first([mem, 2000 + ceil(2000 * (vcfs_gb + pileups_gb))])
    Int command_mem = machine_mem - 500
    Int disk_needed = ceil(2 * vcfs_gb + pileups_gb) + 5 + runtime_params.emergency_extra_disk
    Int disk = select_first([disk_space, if disk_needed > 10 then disk_needed else 10])

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
        gatk --java-options "-Xmx~{command_mem}m" MergeVcfs -I ~{sep=' -I ' input_vcfs} -O ~{output_vcf}

        gatk --java-options "-Xmx~{command_mem}m" MergeMutectStats \
            -stats ~{sep=" -stats " stats} -O merged.stats

        if [ ~{length(tumor_pileups)} -gt 0 ]; then
            gatk --java-options "-Xmx~{command_mem}m" GatherPileupSummaries \
                --sequence-dictionary ~{ref_dict} -I ~{sep=' -I ' tumor_pileups} -O tumor-pileups.tsv
        fi

        if [ ~{length(normal_pileups)} -gt 0 ]; then
            gatk --java-options "-Xmx~{command_mem}m" GatherPileupSummaries \
                --sequence-dictionary ~{ref_dict} -I ~{sep=' -I ' normal_pileups} -O normal-pileups.tsv
        fi
    >>>

    runtime {
        docker: runtime_params.gatk_docker
        bootDiskSizeGb: 4
        memory: machine_mem + " MB"
        disks: "local-disk " + disk + " SSD"
        preemptible: preemptible
        maxRetries: runtime_params.max_retries
        cpu: 1
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
      Array[File]+ bam_outs
      Runtime runtime_params
      Int preemptible = 2
      Int? mem
      Int? disk_space
    }

    # memory in MB, enough for the 500,000 reads SortSam keeps in memory, and disk in GB, from the bamouts (sorting
    # used 4.2 times their size)
    Int machine_mem = select_first([mem, 4000])
    Int command_mem = machine_mem - 500
    Int disk_needed = ceil(5 * size(bam_outs, "GB")) + 5 + runtime_params.emergency_extra_disk
    Int disk = select_first([disk_space, if disk_needed > 10 then disk_needed else 10])

    command <<<
        # This command block assumes that there is at least one file in bam_outs.
        #  Do not call this task if len(bam_outs) == 0
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}
        gatk --java-options "-Xmx~{command_mem}m" GatherBamFiles \
            -I ~{sep=" -I " bam_outs} -O unsorted.out.bam

        # We must sort because adjacent scatters may have overlapping (padded) assembly regions, hence
        # overlapping bamouts.  SortSam spills to the working disk, which is sized for it, rather than to /tmp on the
        # boot disk.

        gatk --java-options "-Xmx~{command_mem}m" SortSam -I unsorted.out.bam \
            -O bamout.bam --SORT_ORDER coordinate -VALIDATION_STRINGENCY LENIENT --TMP_DIR .
        gatk --java-options "-Xmx~{command_mem}m" BuildBamIndex -I bamout.bam -VALIDATION_STRINGENCY LENIENT
    >>>

    runtime {
        docker: runtime_params.gatk_docker
        bootDiskSizeGb: 4
        memory: machine_mem + " MB"
        disks: "local-disk " + disk + " SSD"
        preemptible: preemptible
        maxRetries: runtime_params.max_retries
        cpu: 1
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
      Int preemptible = 3
      Int? mem
      Int? disk_space
    }

    # memory in MB and disk in GB, from the F1R2 counts (for a genome's 80 MB, a heap of 150 MB; they extract to 5 times
    # their size)
    Float f1r2_gb = size(f1r2_tar_gz, "GB")
    Int machine_mem = select_first([mem, 2000 + ceil(4000 * f1r2_gb)])
    Int command_mem = machine_mem - 500
    Int disk_needed = ceil(6 * f1r2_gb) + 5 + runtime_params.emergency_extra_disk
    Int disk = select_first([disk_space, if disk_needed > 10 then disk_needed else 10])

    command {
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}

        gatk --java-options "-Xmx~{command_mem}m" LearnReadOrientationModel \
            -I ~{sep=" -I " f1r2_tar_gz} \
            -O "artifact-priors.tar.gz"
    }

    runtime {
        docker: runtime_params.gatk_docker
        bootDiskSizeGb: 4
        memory: machine_mem + " MB"
        disks: "local-disk " + disk + " SSD"
        preemptible: preemptible
        maxRetries: runtime_params.max_retries
        cpu: 1
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
      Int preemptible = 3
      Int? mem
      Int? disk_space
    }

    # memory in MB and disk in GB, from the pileups (for a genome's 0.43 to 0.54 GB of tumor and normal pileups, a heap
    # of up to 2.2 GB)
    Float pileups_gb = size(tumor_pileups, "GB") + size(normal_pileups, "GB")
    Int machine_mem = select_first([mem, 2000 + ceil(10000 * pileups_gb)])
    Int command_mem = machine_mem - 500
    Int disk_needed = ceil(2 * pileups_gb) + 5 + runtime_params.emergency_extra_disk
    Int disk = select_first([disk_space, if disk_needed > 10 then disk_needed else 10])

    command {
        set -e

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}

        gatk --java-options "-Xmx~{command_mem}m" CalculateContamination -I ~{tumor_pileups} \
        -O contamination.table --tumor-segmentation segments.table ~{"-matched " + normal_pileups} \
        --minimum-population-allele-frequency ~{min_population_af} --maximum-population-allele-frequency ~{max_population_af}
    }

    runtime {
        docker: runtime_params.gatk_docker
        bootDiskSizeGb: 4
        memory: machine_mem + " MB"
        disks: "local-disk " + disk + " SSD"
        preemptible: preemptible
        maxRetries: runtime_params.max_retries
        cpu: 1
    }

    output {
        File contamination_table = "contamination.table"
        File maf_segments = "segments.table"
    }
}

task Filter {
    input {
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
        Int preemptible = 2
        Int? mem
        Int? disk_space
    }

    # memory in MB, from the VCF (for a 0.4 GB VCF, a heap of 0.9 GB), and disk in GB, for the reference, the VCF and
    # the outputs, twice its size
    Float vcf_gb = size(unfiltered_vcf, "GB")
    Int machine_mem = select_first([mem, 2000 + ceil(4000 * vcf_gb)])
    Int command_mem = machine_mem - 500
    Int disk_needed = ceil(size(reference, "GB") + 3 * vcf_gb) + 5 + runtime_params.emergency_extra_disk
    Int disk = select_first([disk_space, if disk_needed > 10 then disk_needed else 10])

    String output_vcf = if compress_vcfs then "filtered.vcf.gz" else "filtered.vcf"
    String output_vcf_idx = output_vcf + if compress_vcfs then ".tbi" else ".idx"

    # The reference is localized rather than streamed.  With read orientation priors FilterMutectCalls reads reference
    # bases around every SNV on each of its passes, and GATK reads a streamed FASTA one 128 kB request at a time: on
    # Terra in November 2026 the task took 52 to 128 minutes per sample with its CPUs 11% to 18% busy.  This may be
    # reversed (reference: {localization_optional: true}, the reference left out of the disk size and a standard disk)
    # once GATK streams FASTAs efficiently.

    command {
        set -e

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}

        gatk --java-options "-Xmx~{command_mem}m" FilterMutectCalls -V ~{unfiltered_vcf} \
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
        bootDiskSizeGb: 4
        memory: machine_mem + " MB"
        disks: "local-disk " + disk + " SSD"
        preemptible: preemptible
        maxRetries: runtime_params.max_retries
        cpu: 1
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
      Int preemptible = 2
      Int? mem
      Int? disk_space
    }

    String output_vcf = if compress_vcfs then "filtered.vcf.gz" else "filtered.vcf"
    String output_vcf_idx = output_vcf +  if compress_vcfs then ".tbi" else ".idx"

    # Memory in MB for the BWA-MEM index image, which BWA maps into memory outside the JVM's heap, and the heap (at most
    # 0.6 GB), and disk in GB for the reference, the image and the output VCF, which is no bigger than the input VCF
    Float image_gb = size(realignment_index_bundle, "GB")
    Int machine_mem = select_first([mem, 3500 + ceil(1000 * image_gb)])
    Int command_mem = machine_mem - ceil(1000 * image_gb) - 1000
    Int disk_needed = ceil(size(reference, "GB") + image_gb + size(input_vcf, "GB")) + 5
        + runtime_params.emergency_extra_disk
    Int disk = select_first([disk_space, if disk_needed > 10 then disk_needed else 10])

    # The reference is localized rather than streamed: FilterAlignmentArtifacts assembles around every variant, and
    # GATK reads a streamed FASTA one 128 kB request at a time, so on Terra in November 2026 it spent about half its
    # time waiting on exomes.  This may be reversed (localization_optional, the reference left out of the disk size
    # and a standard disk) once GATK streams FASTAs efficiently.
    parameter_meta{
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
        bootDiskSizeGb: 4
        memory: machine_mem + " MB"
        disks: "local-disk " + disk + " SSD"
        preemptible: preemptible
        maxRetries: runtime_params.max_retries
        cpu: 1
    }

    output {
        File filtered_vcf = "~{output_vcf}"
        File filtered_vcf_idx = "~{output_vcf_idx}"
    }
}

task Concatenate {
    input {
        Array[File] input_files
        Runtime runtime_params
        Int preemptible = 3
        Int? mem
        Int? disk_space
    }

    # memory in MB, which cat hardly uses, and disk in GB for the files and their concatenation
    Int machine_mem = select_first([mem, 2000])
    Int disk_needed = ceil(2 * size(input_files, "GB")) + 5 + runtime_params.emergency_extra_disk
    Int disk = select_first([disk_space, if disk_needed > 10 then disk_needed else 10])

    command {
        cat ~{sep=' ' input_files} > output.txt
    }

    runtime {
        docker: runtime_params.gatk_docker
        bootDiskSizeGb: 4
        memory: machine_mem + " MB"
        disks: "local-disk " + disk + " SSD"
        preemptible: preemptible
        maxRetries: runtime_params.max_retries
        cpu: 1
    }

    output {
        File concatenated = "output.txt"
    }
}

