#!/bin/sh
BAM=$1
REF=$3
OUT=/tmp/out.3.pr.vcf.gz
DRAGSTR=$2
./gatk --java-options "-agentlib:jdwp=transport=dt_socket,server=y,suspend=n,address=5006 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Xmx6G -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
    HaplotypeCaller \
    -R ${REF} \
    -I ${BAM} \
    -O ${OUT} \
    -gam USE_POSTERIOR_PROBABILITIES \
    -L chr20:10000000-25000000 \
    --pcr-indel-model NONE \
    --dragstr-params-path ${DRAGSTR} \
    -contamination 0 #-ERC GVCF





