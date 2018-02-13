#!/bin/bash

## contains misc helper functions (or templates on how to make a customized one) for reviewing variants

########## for use in reviewing specific assembly contigs that did or did not induce variants
function review() {
    if [[ "$#" -eq 2 ]]; then
        TIG_NAME=$(printf "asm%06u:tig%05u" "$1" "$2")
    else
        TIG_NAME=$1
    fi
    ASSEMBLY_ALIGNMENTS=""
    MASTER_VCF="/Users/shuang/Desktop/callsets/gatk/master.vcf"
    FEATURE_VCF="/Users/shuang/Desktop/callsets/gatk/experimentalVariantInterpretations/nonComplex.vcf"
    AMBIGUOUS_SAM="/Users/shuang/Desktop/callsets/feature/Ambiguous.sam"
    INCOMPLETE_SAM="/Users/shuang/Desktop/callsets/feature/Incomplete.sam"
    grep -F "$TIG_NAME" "$ASSEMBLY_ALIGNMENTS"
    echo;
    echo -e '\033[0;35mmaster\033[0m' && grep -F "$TIG_NAME" "$MASTER_VCF"
    echo -e '\033[0;35mfeature\033[0m' && grep -F "$TIG_NAME" "$FEATURE_VCF"
    echo;
    echo -e '\033[0;35mambiguous\033[0m' && grep -F "$TIG_NAME" "$AMBIGUOUS_SAM"
    echo -e '\033[0;35mincomplete\033[0m' && grep -F "$TIG_NAME" "$INCOMPLETE_SAM"
}

########## error and clean up
on_error() {
    directory=$(pwd)
    echo "ERROR: RUN INTO ERRORS. QUIT"
    rm -rf "$directory/temp*"
    exit 1
}

########## get the correct "sort" and "uniq" function (GNU version "gsort" and "gnuiq" on darwin systems)
SORT="sort"
UNIQ="uniq"
if [[ "$OSTYPE" == "darwin"* ]]; then
        command -v gsort >/dev/null 2>&1 || { echo >&2 "I require gsort for sorting output bed file but it's not installed. Aborting."; exit 1; }
        command -v guniq >/dev/null 2>&1 || { echo >&2 "I require guniq for sorting output bed file but it's not installed. Aborting."; exit 1; }
        SORT="gsort"
        UNIQ="guniq"
elif [[ "$OSTYPE" == "linux-gnu" ]]; then
        SORT="sort"
        UNIQ="uniq"
else
        echo "unsupported OS, quit"
        exit 1
fi

########## get primary contigs pattern given requested assembly version (must be 19 or 38)
function get_primary_contigs_pattern() {
    PRIMARY_CONTIGS_PATTERN=""
    if [[ $1 == "19" ]]; then
        PRIMARY_CONTIGS_PATTERN="^([0-9]{1,2}|X|Y)  "
    elif [[ $1 == "38" ]]; then
        PRIMARY_CONTIGS_PATTERN="^chr([0-9]{1,2}|X|Y)   "
    else
        echo "reference version must be either 19 or 38"
        exit 1
    fi
}
