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
export -f review

########## error and clean up
on_error() {
    directory=$(pwd)
    echo "ERROR: RUN INTO ERRORS. QUIT"
    rm -rf "$directory/temp*"
    exit 1
}
export -f on_error

########## check for required programs
command -v bedtools >/dev/null 2>&1 || { echo >&2 "I require bedtools for interval intersection analysis but it's not installed. Aborting."; exit 1; }
command -v parallel >/dev/null 2>&1 || { echo >&2 "I require parallel but it's not installed. Aborting."; exit 1; }
command -v R >/dev/null 2>&1 || { echo >&2 "I require R but it's not installed. Aborting."; exit 1; }
temp=$(Rscript -e 'packages <-c("stringr", "plyr");  uninstalled <- packages[ which(packages %in% installed.packages()[, "Package"]==F) ]; uninstalled')
if [[ $temp != "character(0)" ]]; then echo "Required R package(s) $temp not installed; Aborting."; exit 1; fi


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
export -f get_primary_contigs_pattern

########## parallel grep for taking pattern from files (naive grep -f is too slow)
function parallel_grep() {

    MAX_MATCH_CNT=$1

    PATTERN_FILE=$2
    TARGET_FILE=$3
    if [[ "$#" -eq 3 ]]; then
        parallel -a "$PATTERN_FILE" "grep -m $MAX_MATCH_CNT -F {} $TARGET_FILE"
    elif [[ "$#" -ge 4 ]]; then
        OUTPUT=$4
        parallel -a "$PATTERN_FILE" "grep -m $MAX_MATCH_CNT -F {} $TARGET_FILE" > "$OUTPUT"
    fi
}
export -f parallel_grep

# function parallel_grep_v() { # this requires target file to be sorted

#     parallel_grep "$1" "$2" temp_parallel_grep.out

#     comm -23 "$2" temp_parallel_grep.out > temp_parallel_grep_v.out 
#     mv temp_parallel_grep_v.out "$3"
# }

function cite() {
    echo "We use a program parallel for fast greping."
    echo "@book{tange_ole_2018_1146014,
      author       = {Tange, Ole},
      title        = {GNU Parallel 2018},
      publisher    = {Ole Tange},
      year         = 2018,
      month        = apr,
      doi          = {10.5281/zenodo.1146014},
      url          = {https://doi.org/10.5281/zenodo.1146014}
    }"
}

export -f cite