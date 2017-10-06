#!/bin/bash

set -eu

if [[ "$#" -lt 6 ]]; then
    echo -e "Please provide:"
    echo -e "  [1] local directory of GATK build (required)"
    echo -e "  [2] cluster name (required)"
    echo -e "  [3] absolute path to the output directory on the cluster (HDFS,required)"
    echo -e "  [4] absolute path to the BAM on the cluster (index assumed accompanying the bam) (HDFS,required)"
    echo -e "  [5] absolute path to the 2 bit reference on the cluster (skip list is assumed accompanying with same basename and extension \".kill.intervals\") (HDFS,required)"
    echo -e "  [6] absolute path to the reference index image on each worker node's local file system (required)"
    echo -e "  [*] extra command-line arguments to StructuralVariationDiscoveryPipelineSpark"
    echo -e "Example:"
    echo -e " bash svDiscover.sh \\"
    echo -e "      ~/GATK/gatk \\"
    echo -e "      my-test-cluster \\"
    echo -e "      /test-sample \\"
    echo -e "      /data/NA12878_test.bam \\"
    echo -e "      /reference/Homo_sapiens_assembly38.2bit"
    echo -e "      /mnt/1/reference/Homo_sapiens_assembly38.fasta.img (the /mnt/1/ prefix is the default mounting location of the ssd if the cluster is created with a local ssd for the worker nodes)"
    exit 1
fi


GATK_DIR="$1"
CLUSTER_NAME="$2"
MASTER_NODE="hdfs://${CLUSTER_NAME}-m:8020"
PROJECT_OUTPUT_DIR="${MASTER_NODE}$3"
INPUT_BAM="${MASTER_NODE}$4"
REF_TWOBIT="${MASTER_NODE}$5"
REF_INDEX_IMAGE="$6"
INTERVAL_KILL_LIST=$(echo "${REF_TWOBIT}" | sed 's/.2bit$/.kill.intervals/')
KMER_KILL_LIST=$(echo "${REF_TWOBIT}" | sed 's/.2bit$/.kill.kmers/')
ALTS_KILL_LIST=$(echo "${REF_TWOBIT}" | sed 's/.2bit$/.kill.alts/')

# extract any extra arguments to StructuralVariationDiscoveryPipelineSpark
shift $(($# < 6 ? $# : 6))
SV_ARGS=${*:-${SV_ARGS:-""}}
# expand any local variables passed as strings (e.g. PROJECT_OUTPUT_DIR)
eval "SV_ARGS=\"${SV_ARGS}\""

# Choose NUM_EXECUTORS = 2 * NUM_WORKERS
# NOTE: this would find preemptible workers, but it produces
# (erroneous?) deprecation warnings
#NUM_WORKERS=$(gcloud compute instances list --filter="name ~ ${CLUSTER_NAME}-[sw].*" | grep RUNNING | wc -l)
# this works but does not see preemptible workers
NUM_WORKERS=$(gcloud dataproc clusters list --filter "clusterName = ${CLUSTER_NAME}" | tail -n 1 | awk '{print $2}')
if [ -z "${NUM_WORKERS}" ]; then
    echo "Cluster \"${CLUSTER_NAME}\" not found"
    exit 1
fi
NUM_EXECUTORS=$((2 * ${NUM_WORKERS}))

"${GATK_DIR}/gatk-launch" StructuralVariationDiscoveryPipelineSpark \
    -I "${INPUT_BAM}" \
    -O "${PROJECT_OUTPUT_DIR}/variants/inv_del_ins.vcf" \
    -R "${REF_TWOBIT}" \
    --alignerIndexImage "${REF_INDEX_IMAGE}" \
    --exclusionIntervals "${INTERVAL_KILL_LIST}" \
    --kmersToIgnore "${KMER_KILL_LIST}" \
    --crossContigsToIgnore "${ALTS_KILL_LIST}" \
    --breakpointIntervals "${PROJECT_OUTPUT_DIR}/intervals" \
    --fastqDir "${PROJECT_OUTPUT_DIR}/fastq" \
    --contigSAMFile "${PROJECT_OUTPUT_DIR}/assemblies.sam" \
    --targetLinkFile ${PROJECT_OUTPUT_DIR}/target_links.bedpe \
    ${SV_ARGS} \
    -- \
    --sparkRunner GCS \
    --cluster "${CLUSTER_NAME}" \
    --num-executors ${NUM_EXECUTORS} \
    --driver-memory 30G \
    --executor-memory 30G \
    --conf spark.yarn.executor.memoryOverhead=5000 \
    --conf spark.network.timeout=600 \
    --conf spark.executor.heartbeatInterval=120
