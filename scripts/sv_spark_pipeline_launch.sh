#!/usr/bin/env bash

if [[ -z ${1+x} ]]; then
    echo "please provide address for output directory you have write access to"
    exit 1
else
    OUTPUT=$1
    echo "output directory is set by user to '$1'"
fi

# STAGE 1 args
ONE_MEGA_BYTE=1000000
TEN_MEGA_BYTES=10000000
EXCLUSIONINTERVALS=hdfs://svdev-1-m:8020/reference/GRCh37.kill.intervals
KMERSTOIGNORE=hdfs://svdev-1-m:8020/reference/Homo_sapiens_assembly38.k51.dups
FASTQDIR=$OUTPUT/fastq
EVIDENCE=$OUTPUT/evidence
INTERVAL=$OUTPUT/intervals
KMERINTERVALS=$OUTPUT/kmerIntervals
QNAMEINTERVALSMAPPED=$OUTPUT/qnameIntervalsMapped
QNAMEINTERVALSFORASSEMBLY=$OUTPUT/qnameIntervalsForAssembly
MAXFASTQSIZE=$TEN_MEGA_BYTES

# STAGE 2 args
SGA=/usr/local/bin/sga
ASSEMBLY=$OUTPUT/assembly
ASSEMBLY_SUCC=$ASSEMBLY"_0"
# STAGE 3 and 4 args
ALIGNMENTS=$OUTPUT/aligned_assemblies
SITEVCFDIR=$OUTPUT/variants
SITEVCF=$SITEVCFDIR"/inversions.vcf"

# STAGE 5 args
GENOTYPEDVCF=$OUTPUT"/genotype/inversions.vcf"
GENOTYPEDEBUG=$OUTPUT"/genotype/rlldebug"

# SHARED args
REFERENCE=hdfs://svdev-1-m:8020/reference/Homo_sapiens_assembly19.fasta
TWOBITREF=hdfs://svdev-1-m:8020/reference/Homo_sapiens_assembly19.2bit
BAM=hdfs://svdev-1-m:8020/data/NA12878_PCR-_30X.bam

SPARKARGS_HIGH="--sparkRunner GCS --cluster svdev-1 --driver-memory 30G --num-executors 20 --executor-memory 30G --conf spark.yarn.executor.memoryOverhead=5000"
SPARKARGS_LOW="--sparkRunner GCS --cluster svdev-1 --driver-memory 30G --executor-memory 8G --conf spark.yarn.executor.memoryOverhead=5000"

./gatk-launch FindBreakpointEvidenceSpark \
    -R $REFERENCE \
    -I $BAM \
    -O $FASTQDIR \
    --exclusionIntervals $EXCLUSIONINTERVALS \
    --kmersToIgnore $KMERSTOIGNORE \
    --kmerIntervals KMERINTERVALS \
    --breakpointIntervals $INTERVAL \
    --breakpointEvidenceDir $EVIDENCE \
    --qnameIntervalsMapped $QNAMEINTERVALSMAPPED \
    --qnameIntervalsForAssembly $QNAMEINTERVALSFORASSEMBLY \
    --maxFASTQSize $MAXFASTQSIZE \
    -- \
    $SPARKARGS_HIGH

./gatk-launch RunSGAViaProcessBuilderOnSpark \
    --fullPathToSGA $SGA \
    --inDir $FASTQDIR \
    --outDirPrefix $ASSEMBLY \
    -- \
    $SPARKARGS_LOW

./gatk-launch AlignAssembledContigsSpark \
    --inputFile $ASSEMBLY_SUCC \
    -R $REFERENCE \
    -O $ALIGNMENTS \
    -- \
    --sparkRunner GCS \
    --cluster svdev-1 \
    --driver-memory 30G \
    --executor-memory 30G \
    --num-executors 20 \
    --conf spark.yarn.executor.memoryOverhead=16000

./gatk-launch CallVariantsFromAlignedContigsSpark \
    -R $TWOBITREF \
    --fastaReference $REFERENCE \
    --inputAssemblies $ASSEMBLY_SUCC \
    --inputAlignments $ALIGNMENTS \
    --outputPath $SITEVCFDIR \
    -- \
    $SPARKARGS_LOW

./gatk-launch SingleDiploidSampleBiallelicInversionGenotyperSpark \
    -I $BAM \
    -R $REFERENCE \
    -fastaReference $REFERENCE \
    -V $SITEVCF \
    -fastq $FASTQDIR \
    -assembly $ASSEMBLY_SUCC \
    -aln $ALIGNMENTS \
    -O $GENOTYPEDVCF \
    -dso $GENOTYPEDEBUG \
    -- \
    --sparkRunner GCS \
    --cluster svdev-1 \
    --driver-memory 30G \
    --num-executors 20 --executor-cores 8 --executor-memory 30G \
    --conf spark.yarn.executor.memoryOverhead=1000