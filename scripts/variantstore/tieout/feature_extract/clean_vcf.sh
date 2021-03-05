REFERENCE=$1
INPUT_VCF=$2
OUTPUT_PREFIX=$3


bcftools view -G $INPUT_VCF > /tmp/${OUTPUT_PREFIX}.so.vcf
bcftools norm -N -m -any /tmp/${OUTPUT_PREFIX}.so.vcf | awk '{ if ($5 != "*") print $0 }' > /tmp/${OUTPUT_PREFIX}.so.split.vcf 

gatk LeftAlignAndTrimVariants \
   -R $REFERENCE \
	 --max-indel-length 1000 \
   -V /tmp/${OUTPUT_PREFIX}.so.split.vcf \
   -O /tmp/${OUTPUT_PREFIX}.so.split.reduced.vcf

bcftools norm -f $REFERENCE \
/tmp/${OUTPUT_PREFIX}.so.split.reduced.vcf > /tmp/${OUTPUT_PREFIX}.so.split.reduced.again.vcf

bcftools sort -O z /tmp/${OUTPUT_PREFIX}.so.split.reduced.again.vcf > ${OUTPUT_PREFIX}.clean.vcf.gz
