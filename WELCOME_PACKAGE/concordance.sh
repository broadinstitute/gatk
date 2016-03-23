#!/bin/bash

GATK3_VCF="$1"
GATK4_VCF="$2"
HELLBENDER_LARGE="/Users/droazen/src/hellbender/src/test/resources/large/"

rm hc_concordance*
rm hc_3discordance4*
rm hc_4discordance3*

java -jar GATK3.4-46.jar -T SelectVariants -R "${HELLBENDER_LARGE}/human_g1k_v37.20.21.fasta" -V "${GATK3_VCF}" --concordance "${GATK4_VCF}" -o hc_concordance.vcf

java -jar GATK3.4-46.jar -T SelectVariants -R "${HELLBENDER_LARGE}/human_g1k_v37.20.21.fasta" -V "${GATK3_VCF}" --discordance "${GATK4_VCF}" -o hc_3discordance4.vcf

java -jar GATK3.4-46.jar -T SelectVariants -R "${HELLBENDER_LARGE}/human_g1k_v37.20.21.fasta" -V "${GATK4_VCF}" --discordance "${GATK3_VCF}" -o hc_4discordance3.vcf

GATK3_COUNT=`cat "${GATK3_VCF}" | grep -v '^#' | wc -l`
GATK4_COUNT=`cat "${GATK4_VCF}" | grep -v '^#' | wc -l`

CONC_COUNT=`cat "hc_concordance.vcf" | grep -v '^#' | wc -l`
DIS3_4_COUNT=`cat "hc_3discordance4.vcf" | grep -v '^#' | wc -l`
DIS4_3_COUNT=`cat "hc_4discordance3.vcf" | grep -v '^#' | wc -l`

echo "GATK3 total: ${GATK3_COUNT}"
echo "GATK4 total: ${GATK4_COUNT}"

echo "Concordant: ${CONC_COUNT}"
echo "3 Discordant 4: ${DIS3_4_COUNT}"
echo "4 Discordant 3: ${DIS4_3_COUNT}"

exit 0
