#!/usr/bin/env bash

# example to run this script
# ./sv_spark_pipeline_launch.sh svdev-1 8020 /user/svshare/NA12878_PCR-_30X/experiments/pipeline

if [[ -z ${1+x} ]]; then
    echo "please provide cluster name to run these jobs on..."
    exit 1
else
    CLUSTER=$1
fi

if [[ -z ${2+x} ]]; then
    echo "please provide port number to connect to the cluster"
    exit 1
else
    MASTER=$CLUSTER"-m:"$2
fi

if [[ -z ${3+x} ]]; then
    echo "please provide address for output directory you have write access to"
    exit 1
else
    OUTPUT="hdfs://"$MASTER"/"$3
    echo "output directory is set by user to "
    echo $OUTPUT
fi

# SHARED args
REFERENCE="hdfs://"$MASTER"/reference/Homo_sapiens_assembly19.fasta"
TWOBITREF="hdfs://"$MASTER"/reference/Homo_sapiens_assembly19.2bit"
BAM="hdfs://"$MASTER"/data/NA12878_PCR-_30X.bam"

# STAGE 1 args
EXCLUSIONINTERVALS="hdfs://"$MASTER"/reference/GRCh37.kill.intervals"
KMERSTOIGNORE="hdfs://"$MASTER"/reference/Homo_sapiens_assembly38.k51.dups"
FASTQDIR=$OUTPUT"/fastq"
EVIDENCE=$OUTPUT"/evidence"
INTERVAL=$OUTPUT"/intervals"
KMERINTERVALS=$OUTPUT"/kmerIntervals"
QNAMEINTERVALSMAPPED=$OUTPUT"/qnameIntervalsMapped"
QNAMEINTERVALSFORASSEMBLY=$OUTPUT"/qnameIntervalsForAssembly"
MAXFASTQSIZE=10000000

# STAGE 2 args
SGA=/usr/local/bin/sga
ASSEMBLY=$OUTPUT"/assembly"
ASSEMBLY_SUCC=$ASSEMBLY"_0"
# STAGE 3 and 4 args
ALIGNMENTS=$OUTPUT"/aligned_assemblies"
SITEVCFDIR=$OUTPUT"/variants"
SITEVCF=$SITEVCFDIR"/inversions.vcf"

# Spark args
SPARKARGS_COMM="--sparkRunner GCS --cluster "$CLUSTER" --driver-memory 30G"
SPARKARGS_LOW=$SPARKARGS_COMM" --executor-memory 8G --conf spark.yarn.executor.memoryOverhead=5000"
SPARKARGS_HIGH=$SPARKARGS_COMM" --num-executors 20 --executor-memory 30G --conf spark.yarn.executor.memoryOverhead=5000"
SPARKARGS_MEMHOG=$SPARKARGS_COMM" --num-executors 20 --executor-memory 30G --conf spark.yarn.executor.memoryOverhead=16000"

 cd ..

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

# arg "subStringToStrip" will be removed once the tool is cleaned up in a PR
./gatk-launch RunSGAViaProcessBuilderOnSpark \
    --fullPathToSGA $SGA \
    --inDir $FASTQDIR \
    --outDirPrefix $ASSEMBLY \
    --subStringToStrip assembly \
    -- \
    $SPARKARGS_LOW

./gatk-launch AlignAssembledContigsSpark \
    --inputFile $ASSEMBLY_SUCC \
    -R $REFERENCE \
    -O $ALIGNMENTS \
    -- \
    $SPARKARGS_MEMHOG

./gatk-launch CallVariantsFromAlignedContigsSpark \
    -R $TWOBITREF \
    --fastaReference $REFERENCE \
    --inputAssemblies $ASSEMBLY_SUCC \
    --inputAlignments $ALIGNMENTS \
    --outputPath $SITEVCFDIR \
    -- \
    $SPARKARGS_LOW