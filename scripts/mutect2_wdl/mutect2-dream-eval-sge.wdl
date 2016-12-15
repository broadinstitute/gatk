# how to run
# 1) start a cromwell server:
# java -Dconfig.file=sge.conf -jar $cromwell server
# 2) open swagger in a web browser. URL: <server>:<port>
# 3) load the wdl and json and 'try it out'.

task SplitIntervals {
  Int scatterCount
  File intervals

  command <<<
    mkdir intervalFileDir

    java -jar /seq/software/picard/current/bin/picard-private.jar IntervalListTools \
    I=${intervals} \
    O=intervalFileDir \
    SCATTER_COUNT=${scatterCount}
  >>>

  output {
    Array[File] intervalFiles = glob("intervalFileDir/temp_*/scattered.interval_list")
  }

}

task M2 {
  File gatk
  File tumor
  File tumorIndex
  String tumorSampleName
  File normal
  File normalIndex
  String normalSampleName
  String vcfBasename
  File intervals

  File pon = "/xchip/cga/reference/hg19/wgs_hg19_125_cancer_blood_normal_panel.vcf"
  File ponIndex = "/xchip/cga/reference/hg19//wgs_hg19_125_cancer_blood_normal_panel.vcf.idx"

  File ref_fasta = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta"
  File ref_fasta_index = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta.fai"
  File ref_dict = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.dict"
  
  File dbsnp = "/humgen/gsa-hpprojects/GATK/bundle/current/b37/dbsnp_138.b37.vcf"
  File dbsnpIndex = "/humgen/gsa-hpprojects/GATK/bundle/current/b37/dbsnp_138.b37.vcf.idx"
  File cosmic = "/xchip/cga/reference/hg19/hg19_cosmic_v54_120711.vcf"
  File cosmicIndex = "/xchip/cga/reference/hg19/hg19_cosmic_v54_120711.vcf.idx"
  
  command <<<
  	java -Xmx4g -jar ${gatk} Mutect2 \
    -R ${ref_fasta} \
    -I ${tumor} \
    -I ${normal} \
    -tumor ${tumorSampleName} \
    -normal ${normalSampleName} \
    --dbsnp ${dbsnp} \
    --cosmic ${cosmic} \
    --normal_panel ${pon} \
    -L ${intervals} \
    -O ${vcfBasename}.vcf
  >>>

  output {
    File outputVcf = "${vcfBasename}.vcf"
  }
}

task EvaluateVcfs {
  File mutect2Vcf
  Int dreamVersion

  File evalScript = "/dsde/working/mutect/dream_smc/dream_eval.pl"
  File snp21="dream${dreamVersion}-chr21-indel.txt"
  File snpwgs="dream${dreamVersion}-wgs-snp.txt"
  File indel21="dream${dreamVersion}-chr21-indel.txt"
  File indelwgs="dream${dreamVersion}-wgs-indel.txt"

  command <<<
    source /humgen/gsa-hpprojects/dev/tsato/wdl/pyvenv/bin/activate

    ${evalScript} ${dreamVersion} 21 SNV ${mutect2Vcf} > ${snp21}
    ${evalScript} ${dreamVersion} 21 INDEL ${mutect2Vcf} > ${indel21}
    ${evalScript} ${dreamVersion} wgs SNV ${mutect2Vcf} > ${snpwgs}
    ${evalScript} ${dreamVersion} wgs INDEL ${mutect2Vcf} > ${indelwgs}

    echo sensitivity,specificity,balanced-accuracy,dream,interval,masked,type > "dream${dreamVersion}.csv"
    grep sensitivity ${snp21} | head -n 1 | cut -d ':' -f 2 | awk -v d=${dreamVersion} '{ print $1 $2 "," d ",chr21," "masked" "," "snp"}' >> "dream${dreamVersion}.csv"
    grep sensitivity ${snpwgs} | tail -n 1 | cut -d ':' -f 2 | awk -v d=${dreamVersion} '{ print $1 $2 "," d ",wgs," "unmasked" "," "snp" }' >> "dream${dreamVersion}.csv"
    grep sensitivity ${indel21} | head -n 1 | cut -d ':' -f 2 | awk -v d=${dreamVersion} '{ print $1 $2 "," d ",chr21," "masked" "," "indel"}' >> "dream${dreamVersion}.csv"
    grep sensitivity ${indelwgs} | tail -n 1 | cut -d ':' -f 2 | awk -v d=${dreamVersion} '{ print $1 $2 "," d ",wgs," "unmasked" "," "indel" }' >> "dream${dreamVersion}.csv"
  >>>

  output {
    File evaluationSummary = "dream${dreamVersion}.csv"
  }
}

# copied from https://github.com/broadinstitute/dsde-pipelines/blob/master/genomes_in_the_cloud/single_sample/PairedSingleSampleWf.wdl
# combines multiple vcfs from scattered Mutect2 runs
task GatherVCFs {
  Array[File] input_vcfs
  String output_vcf_name

  # using MergeVcfs instead of GatherVcfs so we can create indices
  # WARNING 2015-10-28 15:01:48 GatherVcfs  Index creation not currently supported when gathering block compressed VCFs.
  command <<<
    java -Xmx2g -jar /seq/software/picard/current/bin/picard-private.jar \
    MergeVcfs \
    INPUT=${sep=' INPUT=' input_vcfs} \
    OUTPUT=${output_vcf_name}
  >>>

  output {
    File output_vcf = "${output_vcf_name}"
  }
}

workflow M2DreamChallenge  {
  File gatk
  Int scatterCount

  File intervals = "/dsde/working/mutect/dream_smc/bams/wgs_calling_regions.v1.interval_list"
  # note these filepaths are symlinks
  File d3tumor = "/dsde/working/mutect/dream_smc/bams/synthetic.challenge.set3.tumor.bam"
  File d3tumorIndex = "/dsde/working/mutect/dream_smc/bams/synthetic.challenge.set3.tumor.bai"
  File d3normal = "/dsde/working/mutect/dream_smc/bams/synthetic.challenge.set3.normal.bam"
  File d3normalIndex = "/dsde/working/mutect/dream_smc/bams/synthetic.challenge.set3.normal.bai"
  File d4tumor = "/dsde/working/mutect/dream_smc/bams/synthetic.challenge.set4.tumor.bam"
  File d4tumorIndex = "/dsde/working/mutect/dream_smc/bams/synthetic.challenge.set4.tumor.bai"
  File d4normal = "/dsde/working/mutect/dream_smc/bams/synthetic.challenge.set4.normal.bam"
  File d4normalIndex = "/dsde/working/mutect/dream_smc/bams/synthetic.challenge.set4.normal.bai"

  call SplitIntervals {
    input:
      scatterCount = scatterCount,
      intervals = intervals
  }

  scatter (subIntervals in SplitIntervals.intervalFiles ) {
    call M2 as Dream3 {
      input: 
        gatk = gatk,
        tumor = d3tumor,
        tumorIndex = d3tumorIndex,
        tumorSampleName = "IS3.snv.indel.sv",
        normal = d3normal,
        normalIndex = d3normalIndex,
        normalSampleName = "G15512.prenormal.sorted",
        vcfBasename = "d3",
        intervals = subIntervals
    } 
    
    call M2 as Dream4 {
      input: 
        gatk = gatk,
        tumor = d4tumor,
        tumorIndex = d4tumorIndex,
        tumorSampleName = "synthetic.challenge.set4.tumour",
        normal = d4normal,
        normalIndex = d4normalIndex,
        normalSampleName = "synthetic.challenge.set4.normal",
        vcfBasename = "d4",
        intervals = subIntervals
    } 
  }

  call GatherVCFs as Gather3 {
    input:
      input_vcfs = Dream3.outputVcf,
      output_vcf_name = "dream3.vcf"
  }

  call GatherVCFs as Gather4 {
    input:
      input_vcfs = Dream4.outputVcf,
      output_vcf_name = "dream4.vcf"
  }

  call EvaluateVcfs as EvaluateD3 {
    input:
      mutect2Vcf = Gather3.output_vcf,
      dreamVersion = 3
  }

  call EvaluateVcfs as EvaluateD4 {
    input:
      mutect2Vcf = Gather4.output_vcf,
      dreamVersion = 4
  }
}