SAMPLE=$1
POS=$2

GVCF=$(cat legacy_wdl/sample_set_membership_v6.tsv | grep $SAMPLE | cut -f2)
gsutil cat $GVCF | zgrep -C 10 $POS