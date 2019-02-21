import "AlignmentPipeline.wdl" as AlignAndMarkDuplicates
import "MitochondriaCalling.wdl" as MutectAndFilter

workflow MitochondriaPipeline {

  meta {
    description: "Takes in fully aligned hg38 bam and outputs VCF of SNP/Indel calls on the mitochondria."
  }
  parameter_meta {
    wgs_aligned_inpu_bam: "Full WGS hg38 bam or cram"
    autosomal_coverage: "Median coverage of full input bam"
    out_vcf: "Final VCF of mitochondrial SNPs and INDELs"
  }
  File wgs_aligned_input_bam_or_cram
  Int? autosomal_coverage

  File MT_with_numts_intervals

  # Using an older version of the default Mutect LOD cutoff. This value can be changed and is only useful at low depths
  # to catch sites that would not get caught by the LOD divided by depth filter.
  Float lod_cutoff = 6.3
  # Read length used for optimization only. If this is too small CollectWgsMetrics might fail, but the results are not
  # affected by this number. Default is 151.
  Int? max_read_length

  # Full reference is only requred if starting with a CRAM (BAM doesn't need these files)
  File? ref_fasta
  File? ref_fasta_index
  File? ref_dict

  File mt_dict
  File mt_fasta
  File mt_fasta_index
  File mt_amb
  File mt_ann
  File mt_bwt
  File mt_pac
  File mt_sa
  File blacklisted_sites
  File blacklisted_sites_index

  #Shifted reference is used for calling the control region (edge of mitochondria reference).
  #This solves the problem that BWA doesn't support alignment to circular contigs.
  File mt_shifted_dict
  File mt_shifted_fasta
  File mt_shifted_fasta_index
  File mt_shifted_amb
  File mt_shifted_ann
  File mt_shifted_bwt
  File mt_shifted_pac
  File mt_shifted_sa
  File blacklisted_sites_shifted
  File blacklisted_sites_shifted_index

  File shift_back_chain

  File control_region_shifted_reference_interval_list
  File non_control_region_interval_list

  File? gatk_override
  String? m2_extra_args

  #Optional runtime arguments
  Int? preemptible_tries

  call SubsetBam {
    input:
      input_bam = wgs_aligned_input_bam_or_cram,
      interval_list = MT_with_numts_intervals,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      preemptible_tries = preemptible_tries
  }

  call AddOriginalAlignmentTags {
    input:
      input_bam = SubsetBam.output_bam,
      input_bam_index = SubsetBam.output_bai,
      gatk_override = gatk_override,
      preemptible_tries = preemptible_tries
  }

  call RevertSam {
    input:
      input_bam = AddOriginalAlignmentTags.output_bam,
      preemptible_tries = preemptible_tries
  }

  call AlignAndMarkDuplicates.AlignmentPipeline as AlignToMt {
    input:
      input_bam = RevertSam.unmapped_bam,
      mt_dict = mt_dict,
      mt_fasta = mt_fasta,
      mt_fasta_index = mt_fasta_index,
      mt_amb = mt_amb,
      mt_ann = mt_ann,
      mt_bwt = mt_bwt,
      mt_pac = mt_pac,
      mt_sa = mt_sa,
      preemptible_tries = preemptible_tries
  }

  call AlignAndMarkDuplicates.AlignmentPipeline as AlignToShiftedMt {
    input:
      input_bam = RevertSam.unmapped_bam,
      mt_dict = mt_shifted_dict,
      mt_fasta = mt_shifted_fasta,
      mt_fasta_index = mt_shifted_fasta_index,
      mt_amb = mt_shifted_amb,
      mt_ann = mt_shifted_ann,
      mt_bwt = mt_shifted_bwt,
      mt_pac = mt_shifted_pac,
      mt_sa = mt_shifted_sa,
      preemptible_tries = preemptible_tries
  }

  call CollectWgsMetrics {
    input:
      input_bam = AlignToMt.mt_aligned_bam,
      input_bam_index = AlignToMt.mt_aligned_bam,
      ref_fasta = mt_fasta,
      ref_fasta_index = mt_fasta_index,
      read_length = max_read_length,
      coverage_cap = 100000,
      preemptible_tries = preemptible_tries
  }

  call GetContamination {
    input:
      input_bam = AlignToMt.mt_aligned_bam,
      input_bam_index = AlignToMt.mt_aligned_bai,
      ref_fasta = mt_fasta,
      ref_fasta_index = mt_fasta_index,
      preemptible_tries = preemptible_tries
  }

  call MutectAndFilter.MitochondriaCalling as CallAndFilterMt {
    input:
      input_bam = AlignToMt.mt_aligned_bam,
      input_bam_index = AlignToMt.mt_aligned_bai,
      ref_fasta = mt_fasta,
      ref_fasta_index = mt_fasta_index,
      ref_dict = mt_dict,
      lod_cutoff = lod_cutoff,
      gatk_override = gatk_override,
      # Everything is called except the control region.
      m2_extra_args = select_first([m2_extra_args, ""]) + " -L chrM:576-16024 ",
      blacklisted_sites = blacklisted_sites,
      blacklisted_sites_index = blacklisted_sites_index,
      mean_coverage = CollectWgsMetrics.mean_coverage,
      autosomal_coverage = autosomal_coverage,
      contamination = GetContamination.minor_level,
      max_read_length = max_read_length,
      preemptible_tries = preemptible_tries
  }

  call MutectAndFilter.MitochondriaCalling as CallAndFilterShiftedMt {
    input:
      input_bam = AlignToShiftedMt.mt_aligned_bam,
      input_bam_index = AlignToShiftedMt.mt_aligned_bai,
      ref_fasta = mt_shifted_fasta,
      ref_fasta_index = mt_shifted_fasta_index,
      ref_dict = mt_shifted_dict,
      lod_cutoff = lod_cutoff,
      gatk_override = gatk_override,
      # Interval correspondes to control region in the shifted reference
      m2_extra_args = select_first([m2_extra_args, ""]) + " -L chrM:8025-9144 ",
      blacklisted_sites = blacklisted_sites_shifted,
      blacklisted_sites_index = blacklisted_sites_shifted_index,
      mean_coverage = CollectWgsMetrics.mean_coverage,
      autosomal_coverage = autosomal_coverage,
      contamination = GetContamination.minor_level,
      max_read_length = max_read_length,
      preemptible_tries = preemptible_tries
  }

  call LiftoverAndCombineVcfs {
    input:
      shifted_vcf = CallAndFilterShiftedMt.vcf,
      vcf = CallAndFilterMt.vcf,
      ref_fasta = mt_fasta,
      ref_fasta_index = mt_fasta_index,
      ref_dict = mt_dict,
      shift_back_chain = shift_back_chain,
      preemptible_tries = preemptible_tries
  }

  # This is a temporary task to handle "joint calling" until Mutect2 can produce a GVCF.
  # This proivdes coverage at each base so low coverage sites can be considered ./. rather than 0/0.
  call CoverageAtEveryBase {
    input:
      input_bam_regular_ref = AlignToMt.mt_aligned_bam,
      input_bam_regular_ref_index = AlignToMt.mt_aligned_bai,
      input_bam_shifted_ref = AlignToShiftedMt.mt_aligned_bam,
      input_bam_shifted_ref_index = AlignToShiftedMt.mt_aligned_bai,
      shift_back_chain = shift_back_chain,
      control_region_shifted_reference_interval_list = control_region_shifted_reference_interval_list,
      non_control_region_interval_list = non_control_region_interval_list,
      ref_fasta = mt_fasta,
      ref_fasta_index = mt_fasta_index,
      ref_dict = mt_dict,
      shifted_ref_fasta = mt_shifted_fasta,
      shifted_ref_fasta_index = mt_shifted_fasta_index,
      shifted_ref_dict = mt_shifted_dict
  }

  output {
    File subset_bam = SubsetBam.output_bam
    File subset_bai = SubsetBam.output_bai
    File mt_aligned_bam = AlignToMt.mt_aligned_bam
    File mt_aligned_bai = AlignToMt.mt_aligned_bai
    File out_vcf = LiftoverAndCombineVcfs.final_vcf
    File out_vcf_index = LiftoverAndCombineVcfs.final_vcf_index
    File duplicate_metrics = AlignToMt.duplicate_metrics
    File coverage_metrics = CollectWgsMetrics.metrics
    File theoretical_sensitivity_metrics = CollectWgsMetrics.theoretical_sensitivity
    File contamination_metrics = GetContamination.contamination_file
    File base_level_coverage_metrics = CoverageAtEveryBase.table
    Int mean_coverage = CollectWgsMetrics.mean_coverage
    String major_haplogroup = GetContamination.major_hg
    Float contamination = GetContamination.minor_level
  }
}

task SubsetBam {
  File input_bam
  String cram_basename = basename(input_bam, ".cram")
  String basename = basename(cram_basename, ".bam")
  File interval_list
  File? ref_fasta
  File? ref_fasta_index
  File? ref_dict

  # runtime
  Int? preemptible_tries
  Float bam_size = size(input_bam, "GB")
  Float ref_size = if defined(ref_fasta) then size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB") else 0
  Int disk_size = ceil(bam_size * 2 + ref_size) + 20
  Int final_preemptible_tries = if bam_size > 110.0 then 0 else select_first([preemptible_tries, 5])

  meta {
    description: "Subsets a whole genome bam to just NuMT and Mitochondria reads and their mates"
  }
  parameter_meta {
    interval_list: "Interval list with NuMTs and chrM"
    ref_fasta: "Reference is only required for cram input. If it is provided ref_fasta_index and ref_dict are also required."
  }
  command <<<
    java -jar -Xmx2500m /usr/gitc/picard.jar FilterSamReads \
      I=${input_bam} \
      ${"R=" + ref_fasta} \
      O=${basename}.bam \
      FILTER=includePairedIntervals \
      INTERVAL_LIST=${interval_list} \
      CREATE_INDEX=true
  >>>
  runtime {
    memory: "3 GB"
    disks: "local-disk " + disk_size + " HDD"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.1-1540490856"
    preemptible: final_preemptible_tries
  }
  output {
    File output_bam = "${basename}.bam"
    File output_bai = "${basename}.bai"
  }
}

task AddOriginalAlignmentTags {
  File input_bam
  File input_bam_index
  String basename = basename(input_bam, ".bam")
  File? gatk_override

  # runtime
  Int? preemptible_tries
  Int disk_size = ceil(size(input_bam, "GB")) + 20

  meta {
    description: "Adds OriginalAlignment tag (OA) and an Original Mate Contig (XM) tag."
  }
  command <<<
    set -e

    export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

    gatk --java-options "-Xmx2500m" AddOriginalAlignmentTags \
      -I ${input_bam} \
      -O ${basename}.bam
  >>>
  runtime {
      memory: "3 GB"
      disks: "local-disk " + disk_size + " HDD"
      docker: "us.gcr.io/broad-gatk/gatk:4.1.0.0"
      preemptible: select_first([preemptible_tries, 5])
  }
  output {
    File output_bam = "${basename}.bam"
    File output_bai = "${basename}.bai"
  }
}

task RevertSam {
  File input_bam
  String basename = basename(input_bam, ".bam")

  # runtime
  Int? preemptible_tries
  Int disk_size = ceil(size(input_bam, "GB") * 2.5) + 20

  meta {
    description: "Removes alignment information while retaining recalibrated base qualities and original alignment tags"
  }
  parameter_meta {
    input_bam: "aligned bam"
  }
  command {
    java -Xmx1000m -jar /usr/gitc/picard.jar \
    RevertSam \
    INPUT=${input_bam} \
    OUTPUT_BY_READGROUP=false \
    OUTPUT=${basename}.bam \
    VALIDATION_STRINGENCY=LENIENT \
    ATTRIBUTE_TO_CLEAR=FT \
    ATTRIBUTE_TO_CLEAR=CO \
    SORT_ORDER=queryname \
    RESTORE_ORIGINAL_QUALITIES=false
  }
  runtime {
    disks: "local-disk " + disk_size + " HDD"
    memory: "2 GB"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.1-1540490856"
    preemptible: select_first([preemptible_tries, 5])
  }
  output {
    File unmapped_bam = "${basename}.bam"
  }
}

task LiftoverAndCombineVcfs {
  File shifted_vcf
  File vcf
  String basename = basename(shifted_vcf, ".vcf")

  File ref_fasta
  File ref_fasta_index
  File ref_dict

  File shift_back_chain

  # runtime
  Int? preemptible_tries
  Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB")
  Int disk_size = ceil(size(shifted_vcf, "GB") + ref_size) + 20

  meta {
    description: "Lifts over shifted vcf of control region and combines it with the rest of the chrM calls."
  }
  parameter_meta {
    shifted_vcf: "VCF of control region on shifted reference"
    vcf: "VCF of the rest of chrM on original reference"
    ref_fasta: "Original (not shifted) chrM reference"
    shift_back_chain: "Chain file to lift over from shifted reference to original chrM"
  }
  command<<<
    set -e

    java -jar /usr/gitc/picard.jar LiftoverVcf \
      I=${shifted_vcf} \
      O=${basename}.shifted_back.vcf \
      R=${ref_fasta} \
      CHAIN=${shift_back_chain} \
      REJECT=${basename}.rejected.vcf

    java -jar /usr/gitc/picard.jar MergeVcfs \
      I=${basename}.shifted_back.vcf \
      I=${vcf} \
      O=${basename}.final.vcf
    >>>
    runtime {
      disks: "local-disk " + disk_size + " HDD"
      memory: "1200 MB"
      docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.1-1540490856"
      preemptible: select_first([preemptible_tries, 5])
    }
    output{
        # rejected_vcf should always be empty
        File rejected_vcf = "${basename}.rejected.vcf"
        File final_vcf = "${basename}.final.vcf"
        File final_vcf_index = "${basename}.final.vcf.idx"
    }
}

task GetContamination {
  File input_bam
  File input_bam_index
  File ref_fasta
  File ref_fasta_index
  Int qual = 20
  Int map_qual = 30
  Float vaf = 0.01

  String basename = basename(input_bam)

  # runtime
  Int? preemptible_tries
  Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB")
  Int disk_size = ceil(size(input_bam, "GB") + ref_size) + 20

  meta {
    description: "Uses Haplochecker to estimate levels of contamination in mitochondria"
  }
  parameter_meta {
    input_bam: "Bam aligned to chrM"
    ref_fasta: "chrM reference"
  }
  command {
  set -e

  java -jar /usr/mtdnaserver/mitolib-0.1.0.jar haplochecker \
    --in ${input_bam} \
    --ref ${ref_fasta} \
    --out haplochecker_out \
    --QUAL ${qual} \
    --MAPQ ${map_qual} \
    --VAF ${vaf}

python3 <<CODE

import csv

with open("haplochecker_out/${basename}.contamination.txt") as output:
    reader = csv.DictReader(output, delimiter='\t')
    for row in reader:
        print(row["MajorHG"], file=open("major_hg.txt", 'w'))
        print(row["MajorLevel"], file=open("major_level.txt", 'w'))
        print(row["MinorHG"], file=open("minor_hg.txt", 'w'))
        print(row["MinorLevel"], file=open("minor_level.txt", 'w'))
CODE
  }
  runtime {
    preemptible: select_first([preemptible_tries, 5])
    memory: "3 GB"
    disks: "local-disk " + disk_size + " HDD"
    docker: "gatkworkflows/mtdnaserver:1.0"
  }
  output {
    File contamination_file = "haplochecker_out/${basename}.contamination.txt"
    String major_hg = read_string("major_hg.txt")
    Float major_level = read_float("major_level.txt")
    String minor_hg = read_string("minor_hg.txt")
    Float minor_level = read_float("minor_level.txt")
  }
}

task CollectWgsMetrics {
  File input_bam
  File input_bam_index
  File ref_fasta
  File ref_fasta_index
  Int? read_length
  Int read_length_for_optimization = select_first([read_length, 151])
  Int? coverage_cap

  Int? preemptible_tries
  Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB")
  Int disk_size = ceil(size(input_bam, "GB") + ref_size) + 20

  meta {
    description: "Collect coverage metrics"
  }
  parameter_meta {
    read_length: "Read length used for optimization only. If this is too small CollectWgsMetrics might fail. Default is 151."
  }

  command <<<
    set -e

    java -Xms2000m -jar /usr/gitc/picard.jar \
      CollectWgsMetrics \
      INPUT=${input_bam} \
      VALIDATION_STRINGENCY=SILENT \
      REFERENCE_SEQUENCE=${ref_fasta} \
      OUTPUT=metrics.txt \
      USE_FAST_ALGORITHM=true \
      READ_LENGTH=${read_length_for_optimization} \
      ${"COVERAGE_CAP=" + coverage_cap} \
      INCLUDE_BQ_HISTOGRAM=true \
      THEORETICAL_SENSITIVITY_OUTPUT=theoretical_sensitivity.txt

    R --vanilla <<CODE
      df = read.table("metrics.txt",skip=6,header=TRUE,stringsAsFactors=FALSE,sep='\t',nrows=1)
      write.table(floor(df[,"MEAN_COVERAGE"]), "mean_coverage.txt", quote=F, col.names=F, row.names=F)
    CODE
  >>>
  runtime {
    preemptible: select_first([preemptible_tries, 5])
    memory: "3 GB"
    disks: "local-disk " + disk_size + " HDD"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.1-1540490856"
  }
  output {
    File metrics = "metrics.txt"
    File theoretical_sensitivity = "theoretical_sensitivity.txt"
    Int mean_coverage = read_int("mean_coverage.txt")
  }
}

task CoverageAtEveryBase {
  File input_bam_regular_ref
  File input_bam_regular_ref_index
  File input_bam_shifted_ref
  File input_bam_shifted_ref_index
  File shift_back_chain
  File control_region_shifted_reference_interval_list
  File non_control_region_interval_list
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File shifted_ref_fasta
  File shifted_ref_fasta_index
  File shifted_ref_dict

  Int? preemptible_tries
  Int disk_size = ceil(size(input_bam_regular_ref, "GB") + size(input_bam_shifted_ref, "GB") + size(ref_fasta, "GB") * 2) + 20

  meta {
    description: "Remove this hack once there's a GVCF solution."
  }
  command <<<
    set -e

    java -jar /usr/gitc/picard.jar CollectHsMetrics \
      I=${input_bam_regular_ref} \
      R=${ref_fasta} \
      PER_BASE_COVERAGE=non_control_region.tsv \
      O=non_control_region.metrics \
      TI=${non_control_region_interval_list} \
      BI=${non_control_region_interval_list} \
      COVMAX=20000 \
      SAMPLE_SIZE=1

    java -jar /usr/gitc/picard.jar CollectHsMetrics \
      I=${input_bam_shifted_ref} \
      R=${shifted_ref_fasta} \
      PER_BASE_COVERAGE=control_region_shifted.tsv \
      O=control_region_shifted.metrics \
      TI=${control_region_shifted_reference_interval_list} \
      BI=${control_region_shifted_reference_interval_list} \
      COVMAX=20000 \
      SAMPLE_SIZE=1

    R --vanilla <<CODE
      shift_back = function(x) {
        if (x < 8570) {
          return(x + 8000)
        } else {
          return (x - 8569)
        }
      }

      control_region_shifted = read.table("control_region_shifted.tsv", header=T)
      shifted_back = sapply(control_region_shifted[,"pos"], shift_back)
      control_region_shifted[,"pos"] = shifted_back

      beginning = subset(control_region_shifted, control_region_shifted[,'pos']<8000)
      end = subset(control_region_shifted, control_region_shifted[,'pos']>8000)

      non_control_region = read.table("non_control_region.tsv", header=T)
      combined_table = rbind(beginning, non_control_region, end)
      write.table(combined_table, "per_base_coverage.tsv", row.names=F, col.names=T, quote=F, sep="\t")

    CODE
  >>>
  runtime {
    disks: "local-disk " + disk_size + " HDD"
    memory: "1200 MB"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.1-1540490856"
    preemptible: select_first([preemptible_tries, 5])
  }
  output {
    File table = "per_base_coverage.tsv"
  }
}
