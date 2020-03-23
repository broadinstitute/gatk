#!/usr/bin/env bash

# This script takes a bam (plus context), rolls back the effects of MarkDuplicates and BQSR, and runs the BAM
# through Picard/GATK 3 and GATK 4 and compares the legacy results with GATK 4.
#
# This requires that PICARD_JAR and GATK3_JAR be set.
#
# Example usage
# ./validate-reads-spark-pipeline.sh \
#   src/test/resources/org/broadinstitute/hellbender/tools/validation/CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.bam \
#   src/test/resources/large/human_g1k_v37.20.21.fasta \
#   src/test/resources/large/human_g1k_v37.20.21.2bit \
#   src/test/resources/large/dbsnp_138.b37.20.21.vcf \
#   [optionally, throwOnDiff]
#   [optionally, saveIntermediateFiles]

set -ex

# Find the path to the binary
if [ -z "$GATK_HOME" ]; then
  GATK_HOME=`pwd`
fi

if [ -z "$PICARD_JAR" ]; then
  echo "ERROR: PICARD_JAR must be set to the Picard jar, you can get it here: http://broadinstitute.github.io/picard/"  1>&2
  exit 1
fi

if [ -z "$GATK3_JAR" ]; then
  echo "ERROR: GATK3_JAR must be set to the GATK 3 jar, you can get it here: https://www.broadinstitute.org/gatk/download/"  1>&2
  exit 1
fi

if [ "$#" -lt 4 ]; then
    echo "Wrong number of parameters, needs: <your.bam> <reference.fasta> <reference.2bit> <knownSites.vcf> [throwOnDiff] [saveIntermediateBams]"
fi

CRASH_ON_DIFF="false"
if [ "$#" -ge 5 ] && [ "$5" == "throwOnDiff" ]; then
    CRASH_ON_DIFF="true"
fi
if [ "$#" -ge 6 ] && [ "$6" == "throwOnDiff" ]; then
    CRASH_ON_DIFF="true"
fi

SAVE_INTERMEDIATE_FILES="false"
if [ "$#" -ge 5 ] && [ "$5" == "saveIntermediateFiles" ]; then
    SAVE_INTERMEDIATE_FILES="true"
fi
if [ "$#" -ge 6 ] && [ "$6" == "saveIntermediateFiles" ]; then
    SAVE_INTERMEDIATE_FILES="true"
fi

GATK4_JAR="${GATK_HOME}/build/install/gatk/bin/gatk"

TMP_FILE=${1%.bam}
FILE_WO_EXTENSION=${TMP_FILE##*/}


# Copy the bam
TEMP_BAM=tmp.${FILE_WO_EXTENSION}.bam
UNMARKED_BAM=tmp.unmarked.${FILE_WO_EXTENSION}.bam
CLEAN_BAM=tmp.clean.${FILE_WO_EXTENSION}.bam
MARKED_BAM=tmp.marked.${FILE_WO_EXTENSION}.bam
GATK4_BAM=tmp.gatk4.${FILE_WO_EXTENSION}.bam
RECAL_TABLE=tmp.${FILE_WO_EXTENSION}.recal.table

PICARD_MARKED_BAM=tmp.picard.marked.${FILE_WO_EXTENSION}.bam
LEGACY_RECAL_TABLE=tmp.legacy.${FILE_WO_EXTENSION}.recal.table
LEGACY_BAM=tmp.legacy.${FILE_WO_EXTENSION}.bam

cp $1 ${TEMP_BAM}
echo ${TEMP_BAM}



echo "################################################"
echo "###### Undoing MarkDuplicates and BQSR #########"
echo "################################################"
# Undo MarkDuplicates
${GATK4_JAR} UnmarkDuplicates -I ${TEMP_BAM} -O ${UNMARKED_BAM}

# Undo BQSR
${GATK4_JAR} RevertBaseQualityScores -I ${UNMARKED_BAM} -O ${CLEAN_BAM}


echo "################################################"
echo "############## Run GATK 4 tools ################"
echo "################################################"
# Run MarkDuplicates
${GATK4_JAR} MarkDuplicatesSpark -I ${CLEAN_BAM} -O ${MARKED_BAM} --metrics-file=tmp.metrics -DS SUM_OF_BASE_QUALITIES

# Run BQSR
${GATK4_JAR} BaseRecalibratorSpark -I ${MARKED_BAM} -R $3 -knownSites $4 -O ${RECAL_TABLE}
${GATK4_JAR} ApplyBQSRSpark -bqsr_recal_file ${RECAL_TABLE} -I ${MARKED_BAM} -O ${GATK4_BAM}



echo "################################################"
echo "########### Running legacy tools ###############"
echo "################################################"
# Run Picard MarkDuplicates
java -jar ${PICARD_JAR} MarkDuplicates I=${CLEAN_BAM} O=${PICARD_MARKED_BAM} METRICS_FILE=tmp.picard.metrics VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true PROGRAM_RECORD_ID=null DS=SUM_OF_BASE_QUALITIES

# Run GATK 3 BQSR
java -jar ${GATK3_JAR} -T BaseRecalibrator -I ${PICARD_MARKED_BAM} -R $2 -knownSites $4 -o ${LEGACY_RECAL_TABLE}
java -jar ${GATK3_JAR} -T PrintReads -R $2 -BQSR ${LEGACY_RECAL_TABLE} -I ${PICARD_MARKED_BAM} -o ${LEGACY_BAM} -DIQ

echo "################################################"
echo "############### Run diff tools #################"
echo "################################################"
${GATK4_JAR} CompareDuplicatesSpark -I ${LEGACY_BAM} -I2 ${GATK4_BAM} -cd ${CRASH_ON_DIFF}

${GATK4_JAR} CompareBaseQualities -cd ${CRASH_ON_DIFF} --VALIDATION_STRINGENCY SILENT ${LEGACY_BAM} ${GATK4_BAM} 

if [ "$SAVE_INTERMEDIATE_FILES" == "false" ]; then
    echo "deleting temp results"
    rm tmp.*.bam tmp.*.bai tmp.metrics tmp.picard.metrics .tmp*
fi
