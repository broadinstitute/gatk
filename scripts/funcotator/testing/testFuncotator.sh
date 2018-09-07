#!/bin/bash

################################################################################

#Setup variables for the script:
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SCRIPTNAME=$( echo $0 | sed 's#.*/##g' )
MINARGS=0
MAXARGS=1

GATKDIR=${SCRIPTDIR}/../../../

################################################################################

doUnitTests=false
doManualRun=false

doClean=false

REF_VER=hg19
OUT_FORMAT=VCF
useAOUDataSources=false

DATA_SOURCES_PATH=/Users/jonn/Development/funcotator_dataSources_latest
HG19=/Users/jonn/Development/references/Homo_sapiens_assembly19.fasta
HG38=/Users/jonn/Development/references/Homo_sapiens_assembly38.fasta

################################################################################

#Read args:
while [ $# -gt 0 ] ; do

  case "$1" in
		-c)
		doClean=true
		;;
		-u)
		doUnitTests=true
		;;
		-19)
		REF_VER=hg19
		;;
		-38)
		REF_VER=hg38
		;;
    -MAF)
		OUT_FORMAT=MAF
		;;
		-AOU)
		useAOUDataSources=true
		;;
		-t)
		doManualRun=true
		;;
		-r)
			shift
			REF_VER="${1}"
		;;
    *)
    ;; 
  esac
  
  #Get next argument in $1:
  shift
done

################################################################################

r=1
if ${doClean} ; then
	./gradlew clean compileJava compileTestJava installDist
	r=$?
else
	./gradlew compileJava compileTestJava installDist
	r=$?
fi 

if [[ $r -eq 0 ]] && ${doUnitTests} ; then
	echo "################################################################################"
	echo "## Running Unit Tests... "
	./gradlew test --tests org.broadinstitute.hellbender.tools.funcotator* --stacktrace  
	r=$?
fi

################################################################################

if [[ $r -eq 0 ]] && ${doManualRun} ; then
	
	echo "################################################################################"
	echo "## Running Large Tests... "
	echo
	echo
	echo "########################################"
	echo "## Using Reference: ${REF_VER}              ##"
	echo "########################################"

	if [[ "${REF_VER}" == "hg19" ]] ; then
		INPUT=/Users/jonn/Development/NON_PUBLIC/0816201804HC0_R01C01.vcf
		INPUT=/Users/jonn/Development/gatk/src/test/resources/large/funcotator/regressionTestVariantSet1.vcf
		#INPUT=/Users/jonn/Development/gatk/src/test/resources/large/funcotator/regressionTestVariantSet2.vcf
		#INPUT=/Users/jonn/Development/gatk/hg38_trio_liftoverb37.vcf
		#INPUT=/Users/jonn/Development/gatk/tmp.vcf
		#INPUT=/Users/jonn/Development/data_to_run/problem_samples/splice_site_should_not_be_splice_site/error_case.vcf
		
		#HG19=/Users/jonn/Development/references/ucsc.hg19.fasta
		#HG19=/Users/jonn/Development/references/ucsc.hg19.fasta
		#HG19=/Users/jonn/Development/references/GRCh37.p13.genome.fasta
		REF=$HG19
	else
		INPUT=/Users/jonn/Development/FUNCOTATOR_LARGE_TEST_INPUTS/hg38_trio.vcf
		#INPUT=/Users/jonn/Development/gatk/src/test/resources/large/funcotator/regressionTestVariantSetHG38.vcf
		REF=$HG38
	fi

	# Use the AOU data sources if we need them:
	$useAOUDataSources && DATA_SOURCES_PATH=/Users/jonn/Development/funcotator_dataSources.vAoU3

	OUT_FORMAT_LOWER=$( echo "${OUT_FORMAT}" | tr 'A-Z' 'a-z' )
	OUT_FILE_NAME=FUNCOTATOR_OUT.${OUT_FORMAT_LOWER}

	./gatk Funcotator \
		-V ${INPUT} \
		-O ${OUT_FILE_NAME} \
		-R ${REF} \
		--verbosity DEBUG \
		--data-sources-path ${DATA_SOURCES_PATH} \
		--ref-version ${REF_VER} \
		--output-file-format ${OUT_FORMAT} -- --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true'
	
	r=$?
fi

exit $r

