version 1.0

#  Create a Mutect2 panel of normals
#
#  Description of inputs
#  intervals: files of the genomic intervals to call, whose union the workflow calls, or [] for the whole genome
#  ref_fasta: reference genome FASTA, with its index (.fai) and sequence dictionary (.dict) beside it under the names
#  GATK expects, which mutect2.wdl describes
#  normal_bams, normal_bais: the normals' BAMs or CRAMs and their indices, in the same order.  mutect2.wdl takes each
#  BAM or CRAM and its index as its tumor_reads.
#  germline_resource: a germline resource with population allele frequencies, such as gnomAD's, with which
#  CreateSomaticPanelOfNormals excludes germline variants from the panel: an array of the VCF and, optionally, its
#  index, [vcf] or [vcf, idx], or [] to go without one.  mutect2.wdl's "VCFs and their indices" describes the index.
#  (Earlier versions took the VCF and its index as the inputs gnomad and gnomad_idx.)
#  scatter_count: number of parallel jobs when scattering over intervals
#  pon_name: the resulting panel of normals is {pon_name}.vcf
#  m2_extra_args: additional command line parameters for Mutect2.  This should not involve --max-mnp-distance,
#  which the wdl hard-codes to 0 because GenomicsDBImport can't handle MNPs

import "mutect2.wdl" as m2
workflow Mutect2_Panel {
  input {
    Array[File] intervals
    File ref_fasta
    Int scatter_count
    Array[String] normal_bams
    Array[String] normal_bais
    Array[File] germline_resource
    String? m2_extra_args
    String? create_pon_extra_args
    Boolean compress = false
    String pon_name

    Int min_contig_size = 1000000
    Int create_panel_scatter_count = 24

    String? gcs_project_for_requester_pays

    # runtime
    String gatk_docker
    File? gatk_override
    String basic_bash_docker = "ubuntu:16.04"

    Int max_retries = 2
    # GB for MergeVCFs' memory and for CreatePanel's and MergeVCFs' disks; the Mutect2 workflow's tasks size their own
    Int small_task_mem = 4
    Int small_task_disk = 100
  }

  # The reference FASTA and its companion files, which GATK finds by name in the FASTA's directory: its index (.fai),
  # its sequence dictionary (.dict) and, if it's bgzipped, its block index (.gzi).  See mutect2.wdl, which derives them
  # the same way.
  String ref_fasta_path = "~{ref_fasta}"
  Array[String] ref_gzi = if sub(ref_fasta_path, "\\.gz$", "") != ref_fasta_path then [ref_fasta_path + ".gzi"] else []
  Array[File] reference = flatten([[ref_fasta_path, ref_fasta_path + ".fai",
      sub(ref_fasta_path, "\\.(fasta|fas|fa|fna|txt)(\\.gz)?$", "") + ".dict"], ref_gzi])

  # the germline resource with its index, derived as mutect2.wdl derives the indices of its VCFs: given only the VCF,
  # its path plus .tbi if it's bgzipped or .idx if not
  String germline_resource_path = if length(germline_resource) > 0 then "~{germline_resource[0]}" else ""
  Array[File] germline_resource_with_index = if length(germline_resource) == 1
      then [germline_resource_path, germline_resource_path +
          (if sub(germline_resource_path, "\\.b?gz$", "") != germline_resource_path then ".tbi" else ".idx")]
      else germline_resource

  Runtime standard_runtime = {"gatk_docker": gatk_docker, "gatk_override": gatk_override,
            "max_retries": max_retries, "emergency_extra_disk": 0}

    scatter (normal_bam in zip(normal_bams, normal_bais)) {
        call m2.Mutect2 {
            input:
                intervals = intervals,
                ref_fasta = ref_fasta,
                tumor_reads = [normal_bam.left, normal_bam.right],
                normal_reads = [],
                # each normal's calls, without the resources that would remove artifacts and germline variants from them
                germline_resource = [],
                pon = [],
                common_variants_vcf = [],
                scatter_count = scatter_count,
                m2_extra_args = select_first([m2_extra_args, ""]) + "--max-mnp-distance 0",
                gatk_override = gatk_override,
                gatk_docker = gatk_docker,
                max_retries = max_retries,
                gcs_project_for_requester_pays = gcs_project_for_requester_pays
        }
    }

    call m2.SplitIntervals {
        input:
            # the whole genome
            intervals = [],
            reference = reference,
            scatter_count = create_panel_scatter_count,
            split_intervals_extra_args = "--dont-mix-contigs --min-contig-size " + min_contig_size,
            runtime_params = standard_runtime
    }

    scatter (subintervals in SplitIntervals.interval_files ) {
            call CreatePanel {
                input:
                    input_vcfs = Mutect2.output_vcf,
                    intervals = subintervals,
                    reference = reference,
                    germline_resource = germline_resource_with_index,
                    output_vcf_name = pon_name,
                    create_pon_extra_args = create_pon_extra_args,
                    runtime_params = standard_runtime,
                    disk = small_task_disk
            }
    }

    call MergeVCFs {
        input:
            input_vcfs = CreatePanel.output_vcf,
            input_vcf_indices = CreatePanel.output_vcf_index,
            compress_vcfs = compress,
            runtime_params = standard_runtime,
            mem = small_task_mem,
            disk = small_task_disk
    }

    output {
        File pon = MergeVCFs.merged_vcf
        File pon_idx = MergeVCFs.merged_vcf_idx
        Array[File] normal_calls = Mutect2.output_vcf
        Array[File] normal_calls_idx = Mutect2.output_vcf_idx
    }
}

task CreatePanel {
    input {
      File intervals
      Array[String] input_vcfs
      # the reference FASTA, then its .fai index, its .dict dictionary and, if it's bgzipped, its .gzi block index
      Array[File] reference
      String output_vcf_name
      # the germline resource with its index, or [] for none
      Array[File] germline_resource
      String? create_pon_extra_args

      # runtime
      Runtime runtime_params
      Int preemptible = 2
      Int disk
    }

    Int machine_mem = 8
    Int command_mem = machine_mem - 1

        parameter_meta{
          germline_resource: {localization_optional: true}
        }

    command {
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}

        gatk GenomicsDBImport --genomicsdb-workspace-path pon_db -R ~{reference[0]} -V ~{sep=' -V ' input_vcfs} -L ~{intervals}

        gatk --java-options "-Xmx~{command_mem}g"  CreateSomaticPanelOfNormals -R ~{reference[0]} \
            ~{if length(germline_resource) > 0 then "--germline-resource " + germline_resource[0] else ""} \
            -V gendb://pon_db -O ~{output_vcf_name}.vcf ~{create_pon_extra_args}
    }

    runtime {
        docker: runtime_params.gatk_docker
        bootDiskSizeGb: 4
        memory: machine_mem + " GB"
        disks: "local-disk " + disk + " HDD"
        preemptible: preemptible
        maxRetries: runtime_params.max_retries
        cpu: 1
    }

    output {
        File output_vcf = "~{output_vcf_name}.vcf"
        File output_vcf_index = "~{output_vcf_name}.vcf.idx"
    }
}

task MergeVCFs {
    input {
      Array[File] input_vcfs
      Array[File] input_vcf_indices
      Boolean compress_vcfs
      Runtime runtime_params
      Int preemptible = 3
      # memory and disk in GB
      Int mem
      Int disk
    }

    String output_vcf = if compress_vcfs then "merged.vcf.gz" else "merged.vcf"
    String output_vcf_idx = output_vcf + if compress_vcfs then ".tbi" else ".idx"

    # using MergeVcfs instead of GatherVcfs so we can create indices
    # WARNING 2015-10-28 15:01:48 GatherVcfs  Index creation not currently supported when gathering block compressed VCFs.
    command {
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}
        gatk --java-options "-Xmx~{mem * 1000 - 500}m" MergeVcfs -I ~{sep=' -I ' input_vcfs} -O ~{output_vcf}
    }

    runtime {
        docker: runtime_params.gatk_docker
        bootDiskSizeGb: 4
        memory: mem + " GB"
        disks: "local-disk " + disk + " HDD"
        preemptible: preemptible
        maxRetries: runtime_params.max_retries
        cpu: 1
    }

    output {
        File merged_vcf = "~{output_vcf}"
        File merged_vcf_idx = "~{output_vcf_idx}"
    }
}
