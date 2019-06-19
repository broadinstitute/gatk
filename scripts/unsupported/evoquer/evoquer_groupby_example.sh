#!/bin/bash
#
# Example of running Evoquer in "group by" mode
# Must be run from the root of a GATK clone

if [ $# -ne 2 -a $# -ne 3 ]
then
    echo "Usage: ./evoquer.sh dataset_map_file interval [\"extra args\"]"
    exit 1
fi

./gatk Evoquer \
--project-id broad-dsp-spec-ops \
-L "$2" \
-R "src/test/resources/large/Homo_sapiens_assembly38.fasta.gz" \
--dataset-map "$1" \
-O evoquer.vcf \
$3

exit 0
