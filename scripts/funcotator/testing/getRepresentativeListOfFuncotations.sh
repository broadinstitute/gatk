#!/usr/bin/env bash

################################################################################
#
# WARNING: THIS SCRIPT IS UNSUPPORTED!
# USE AT YOUR OWN RISK
#
# DESCRIPTION:
#
# This script creates a new VCF file by extracting a set of N variants from a 
# given Funcotator-annotated VCF file of each type of VARIANT_CLASSIFICATION
# (N is internally configured). 
#
# EXAMPLE:
#     ./getRepresentativeListOfFuncotations.sh VCF_FILE.vcf
#
# AUTHOR: Jonn Smith
#
###############################################################################

inputFile=$1

[[ ! -f $inputFile ]] && echo "Error: Input file does not exist: ${inputFile}" 1>&2 && exit 1

maxResults=100
variantClassifications="COULD_NOT_DETERMINE FIVE_PRIME_FLANK START_CODON_INS DE_NOVO_START_IN_FRAME DE_NOVO_START_OUT_FRAME FIVE_PRIME_UTR FRAME_SHIFT_DEL FRAME_SHIFT_INS IGR INTRON IN_FRAME_DEL IN_FRAME_INS LINCRNA MISSENSE NONSENSE NONSTOP RNA SILENT SPLICE_SITE START_CODON_DEL START_CODON_SNP THREE_PRIME_UTR" 

let numVcs=$( echo "${variantClassifications}" | tr ' ' '\n' | wc -l | sed 's#^[\t ]*##g' )

rawHeaderFile=$( mktemp )
tmpHeaderFile=$( mktemp )
tmpCountFile=$( mktemp )
tmpVcfFile=$( mktemp )
tmpSortedVcFile=$( mktemp )
tmpFile2=$( mktemp )

#rawHeaderFile=rawHeaderFile.txt
#tmpHeaderFile=tmpHeaderFile.txt
#tmpCountFile=tmpCountFile.txt
#tmpVcfFile=tmpVcfFile.txt
#tmpSortedVcFile=tmpSortedVcfFile.txt
#tmpFile2=tmpFile2.txt

#rm -f ${tmpHeaderFile} ${tmpVcfFile} ${rawHeaderFile} ${tmpFile2} ${tmpCountFile}

echo "Beginning extraction of up to ${maxResults} of each variant classification (total #VCs: ${numVcs})..." 1>&2

# Get the VCF header:
echo -n '    Copying Header Core ... ' 1>&2
grep '^#' ${inputFile} > ${rawHeaderFile}
grep -v '#contig=<' ${rawHeaderFile} | grep -v '^#CHROM' > $tmpHeaderFile
echo 'DONE!' 1>&2

# Get the variants:
i=1
for variantClassification in ${variantClassifications} ; do
	printf "    Processing %-23s (%2d of %d) ... " "${variantClassification}" "$i" "$numVcs" 1>&2
	grep "|${variantClassification}|" ${inputFile} | grep -v "^#" | head -n ${maxResults} > $tmpCountFile
	numFound=$( wc -l ${tmpCountFile} | awk '{print $1}' )
	cat ${tmpCountFile} >> ${tmpVcfFile} 
	printf 'DONE%s (Found: %2d)\n' '!' "${numFound}" 1>&2
	let i=$i+1
done

# Sort VCF records:
sort -nk2 ${tmpVcfFile} > ${tmpSortedVcFile}

# Get the contig names from the variants: 
grep -v "^#" ${tmpVcfFile} | awk '{print $1}' | sort | uniq > ${tmpFile2}

## check for chr in contig names in the sequence dictionary:
#grep '^#*contig=<ID' ${rawHeaderFile} | grep 'ID=chr' > /dev/null
#r=$?
#seqDictHasChr=false
#[ $r -eq 0 ] && seqDictHasChr=true
#
## Check for chr in contig names in the variants:
#grep -v '^#' ${tmpFile2} | head -n 1 | awk '{print $1}' | grep '^chr'  > /dev/null
#r=$?
#variantHasChr=false
#[ $r -eq 0 ] && variantHasChr=true
#
## Now generate the contigs in order:
## Get new contig values if the contig names don't match:
#if $seqDictHasChr ^ $variantHasChr; then
#	grep -v "^#" ${tmpVcfFile} | awk '{print $1}' | sort -tr -nk2 | uniq > ${tmpFile2}
#fi
#
#while read contig ; do 
#	if $seqDictHasChr && ! $variantHasChr ; then
#		grep "#*contig=<ID=chr${contig}," ${rawHeaderFile}
#	elif ! $seqDictHasChr && $variantHasChr ; then
#		contig=$(echo ${contig} | sed 's#chr##g' )
#		grep "#*contig=<ID=${contig}," ${rawHeaderFile}
#	else
#		grep "#*contig=<ID=${contig}," ${rawHeaderFile}
#	fi
#done < ${tmpFile2} >> ${tmpHeaderFile} 

# Add the sequence dictionary to the new file:
grep "#*contig=<ID=" ${rawHeaderFile} >> ${tmpHeaderFile}

# Assemble our output:
cat ${tmpHeaderFile} 
echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' 
while read contig ; do 
	grep "^${contig}\t" ${tmpSortedVcFile}
done < ${tmpFile2} | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t.\tPASS\t."}' 

# Clean our garbage - we are not barbarians.
rm -f ${tmpHeaderFile} ${tmpVcfFile} ${rawHeaderFile} ${tmpFile2} ${tmpCountFile}
echo 'Extraction Complete!' 1>&2

