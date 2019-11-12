#!/bin/sh
BAM=$1
REF=$3
OUT=/tmp/out.2.g.vcf.gz
DRAGSTR=$2
./gatk --java-options "-agentlib:jdwp=transport=dt_socket,server=y,suspend=n,address=5006 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Xmx6G -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
    HaplotypeCaller \
    -R ${REF} \
    -I ${BAM} \
    -O ${OUT} \
    --pcr-indel-model NONE \
    -contamination 0 -ERC GVCF



#    dd--dragstr-params-path ${DRAGSTR} \
