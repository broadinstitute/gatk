#!/bin/bash
#

if [ $# -ne 1 -a $# -ne 2 ]
then
    echo "Usage: ./evoquer.sh interval [\"extra args\"]"
    exit 1
fi

./gatk Evoquer \
--project-id broad-dsp-spec-ops \
-L "$1" \
-R "src/test/resources/large/Homo_sapiens_assembly38.fasta.gz" \
--dataset-map evoquer_dataset_map \
-O evoquer.vcf \
$2

exit 0
