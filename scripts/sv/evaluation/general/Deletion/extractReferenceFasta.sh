#!/bin/bash

echo -e "#################################################"
echo -e "  extracting reference context"
echo -e "#################################################"
echo -e "  The FASTA file is a dictionary file where the sequence is the reference sequence being deleted and"
echo -e "  the sequence name is the corresponding variant ID"

REFERECE_FASTA="/Users/shuang/igv/genomes/hg38.fa"

echo -e "  We depend on bedtools for doing it, and this step might take sometime."
bedtools --version

echo -e "\n########## GATK\n"

bedtools getfasta -fi "$REFERECE_FASTA" -bed gatk.cleanDel.bed -fo temp.fasta
# extract ref seq on even lines only
sed '1d; n; d' temp.fasta > temp.refSeq.txt
grep -v '^#' GATK_primaryContigs_cleanDel.vcf | awk '{print $3}' > temp.id.txt
paste -d '\n' temp.id.txt temp.refSeq.txt > temp.gatk.idRef.dictionary.fasta
sed 's/^/>/;n' < temp.gatk.idRef.dictionary.fasta> gatk.idRef.dictionary.fasta
rm temp* 
echo -e "\nDone for GATK."


echo -e "\n########## Manta\n"

# dictionary for mapping from variant ID to reference sequence
grep -v '^#' Manta103_PASS_PRECISE_nonBND_primaryContigs_cleanDel.vcf | \
    grep -v '	<DEL>	' | \
    awk 'BEGIN {OFS="\n"}; {print $3, $4}' > \
    temp.manta.idRef.nonSymb.fasta

# extract reference sequence from bedtools and create FASTA dictionary
grep -v '^#' Manta103_PASS_PRECISE_nonBND_primaryContigs_cleanDel.vcf | \
    grep '	<DEL>	' | >
    awk 'BEGIN {OFS="	"}; {print $1, $2, $8}' | >
    grep -Eo '^([0-9]{1,2}|X|Y)	[0-9]+	END=[0-9]+' | \
    sed 's/END=//' > \
    temp.forBedtools.bed
bedtools getfasta -fi "$REFERECE_FASTA" -bed temp.forBedtools.bed -fo temp.manta.symbDel.refseq.fasta
# extract ref seq on even lines only
sed '1d; n; d' temp.manta.symbDel.refseq.fasta > temp.manta.symbDel.refSeq.txt
grep -v '^#' Manta103_PASS_PRECISE_nonBND_primaryContigs_cleanDel.vcf | \
    grep '	<DEL>	' | awk '{print $3}' > temp.manta.symbDel.ids.txt
paste -d '\n' temp.manta.symbDel.ids.txt temp.manta.symbDel.refSeq.txt \
    > temp.manta.idRef.symb.fasta
########## merge the bed files, sort, merge the FASTA dictionaries and clean up
echo -e "\nNow merge results for Manta symbolic and non-symbolic clean deletions."

cat temp.manta.idRef.nonSymb.fasta temp.manta.idRef.symb.fasta > temp.manta.idRef.dictionary.fasta
sed 's/^/>/;n' < temp.manta.idRef.dictionary.fasta > manta.idRef.dictionary.fasta

rm temp*

echo -e "\nDone for Manta."

echo -e "\n########## PacBio\n"
