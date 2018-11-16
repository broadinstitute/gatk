#!/usr/bin/env bash

################################################################################
#
# WARNING: THIS SCRIPT IS UNSUPPORTED!
# USE AT YOUR OWN RISK
#
# DESCRIPTION:
#
# This script runs Oncotator and will process the output files so that they can
# be ingested by `makeComparisonForOncotator.sh` to compare with Funcotator 
# results.
# It must be internally configured to point at data sources for Oncotator and 
# an input file.  The Oncotator data sources must be analogous to the data 
# sources with which Funcotator will be run.  It must also be told where
# Oncotator lives.
# 
# This will not work for you out-of-the-box unless you are Jonn Smith.
#
# EXAMPLE:
#     ./run_oncotator_VCF_in.sh
#
# AUTHOR: Jonn Smith
#
################################################################################

# Change all thhese paths to point to your oncotator area and the data source
# locations as indicated by the variable names:

TEST_DIR="/Users/jonn/Development/oncotator_testing"

ALL_DB="/Users/jonn/Development/oncotator_v1_ds_April052016/"
QUICK_DB="~/Development/oncotator_testing/dbdir/"
FUNCOTATOR_EQUIVALENT_DB="~/Development/oncotator_testing/funcotator_dbdir/"

VCF_IN="/Users/jonn/Development/gatk/src/test/resources/org/broadinstitute/hellbender/tools/funcotator/validationTestData/regressionTestVariantSet1.vcf"
#VCF_IN="/Users/jonn/Development/gatk/src/test/resources/org/broadinstitute/hellbender/tools/funcotator/validationTestData/regressionTestVariantSet2.vcf"
#VCF_IN="tmp.vcf"
#VCF_IN="/Users/jonn/Development/gatk/src/test/resources/org/broadinstitute/hellbender/tools/funcotator/validationTestData/regressionTestHg19Large.vcf"

source ~/Development/oncotator_venv/bin/activate

doAll=false
doFunk=false
for arg in "${@}" ; do
	if [[ "${arg}" == "-all" ]] ; then
		doAll=true
		echo "USING ALL DATA SOURCES: ${ALL_DB}"
		break
	fi
	if [[ "${arg}" == "-funk" ]] ; then
		doFunk=true
		echo "USING ONLY THE MOST FUNKY OF DATA SOURCES: ${FUNCOTATOR_EQUIVALENT_DB}"
		break
	fi
done

if $doAll ; then
	oncotator -i VCF ${VCF_IN} ${TEST_DIR}/test.maf.raw.tsv hg19 --db-dir ${ALL_DB} 
elif $doFunk ; then
	oncotator -i VCF ${VCF_IN} ${TEST_DIR}/test.maf.raw.tsv hg19 --db-dir ${FUNCOTATOR_EQUIVALENT_DB} 
else
	oncotator -i VCF ${VCF_IN} ${TEST_DIR}/test.maf.raw.tsv hg19 --db-dir ${QUICK_DB} 
fi

#cut -d$'\t' -f 1,5,6,7,9,10,11,12,13,35,36,37,40,41,42,43 ${TEST_DIR}/test.maf.raw.tsv > ${TEST_DIR}/test.maf.tsv 
awk -F$'\t' '{ print $1,$5,$6,$7,$9,$10,$11,$12,$13,$35,$36,$37,$40,$41,$42,$66,$43 }' OFS='\t' ${TEST_DIR}/test.maf.raw.tsv > ${TEST_DIR}/test.maf.tsv 

sed -e '1,4d' ${TEST_DIR}/test.maf.tsv | awk '{print "new Object[] { \""$1"\",", $2 ",", $3 ",", $4 ",", "GencodeFuncotation.VariantClassification." toupper($5) ",", "GencodeFuncotation.VariantType." toupper($6) ",", "\"" toupper($7) "\",", "\"" toupper($8) "\",", "\"" toupper($9) "\",", "\""$10"\",", "\""$11"\",", "\""$12"\",", "\""$13"\"," ,"\""$14"\"", "},"}' | sed -e 's#""#null#g' > test.maf.java 

sed -e "s#De_novo_Start_OutOfFrame#DE_NOVO_START_OUT_FRAME#g" \
    -e "s#Nonsense_Mutation#NONSENSE#g"\
    -e "s#Nonstop_Mutation#NONSTOP#g"\
    -e "s#Missense_Mutation#MISSENSE#g"\
    -e "s#De_novo_Start_InFrame#DE_NOVO_START_IN_FRAME#g"\
    -e "s#In_Frame_Del#IN_FRAME_DEL#g"\
    -e "s#In_Frame_Ins#IN_FRAME_INS#g"\
    -e "s#Frame_Shift_Del#FRAME_SHIFT_DEL#g"\
    -e "s#Frame_Shift_Ins#FRAME_SHIFT_INS#g"\
    -e "s#Start_Codon_SNP#START_CODON_SNP#g"\
    -e "s#Start_Codon_Del#START_CODON_DEL#g"\
    -e "s#Start_Codon_Ins#START_CODON_INS#g"\
    -e "s#Stop_Codon_SNP#STOP_CODON_SNP#g"\
    -e "s#Stop_Codon_Del#STOP_CODON_DEL#g"\
    -e "s#Stop_Codon_Ins#STOP_CODON_INS#g"\
    -e "s#Splice_Site#SPLICE_SITE#g"\
    -e "s#lincRNA#LINCRNA#g"\
    -e "s#Silent#SILENT#g"\
    -e "s#3'UTR#THREE_PRIME_UTR#g"\
    -e "s#5'UTR#FIVE_PRIME_UTR#g"\
    -e "s#Intron#INTRON#g"\
    -e "s#5'Flank#FIVE_PRIME_FLANK#g"\
    -e "s#3'Flank#THREE_PRIME_FLANK#g" ${TEST_DIR}/test.maf.tsv > ${TEST_DIR}/test.maf.FUNCOTATOR_VCLASS.tsv


