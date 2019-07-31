#!/bin/bash

set -eu

if [[ "$#" -lt 6 ]]; then
    echo -e "Please provide:"
    echo -e "  [1] local directory of GATK build (required)"
    echo -e "  [2] project name (required)"
    echo -e "  [3] cluster name (required)"
    echo -e "  [4] absolute path to the output directory on the cluster (HDFS,required)"
    echo -e "  [5] absolute path to the BAM on the cluster (index assumed accompanying the bam) (HDFS,required)"
    echo -e "  [6] absolute path to the 2 bit reference on the cluster (skip list is assumed accompanying with same basename and extension \".kill.intervals\") (HDFS,required)"
    echo -e "  [7] absolute path to the reference index image on each worker node's local file system (required)"
    echo -e "  [*] extra command-line arguments to StructuralVariationDiscoveryPipelineSpark"
    echo -e "Example:"
    echo -e " bash run_whole_pipeline.sh \\"
    echo -e "      ~/GATK/gatk \\"
    echo -e "      my-project \\"
    echo -e "      my-test-cluster \\"
    echo -e "      /test-sample \\"
    echo -e "      /data/NA12878_test.bam \\"
    echo -e "      /reference/Homo_sapiens_assembly38.2bit"
    echo -e "      /reference/Homo_sapiens_assembly38.fasta.img"
    exit 1
fi


GATK_DIR="$1"
PROJECT_NAME="$2"
CLUSTER_NAME="$3"
MASTER_NODE="hdfs://${CLUSTER_NAME}-m:8020"
PROJECT_OUTPUT_DIR="${MASTER_NODE}$4"
INPUT_BAM="${MASTER_NODE}$5"
REF_TWOBIT="${MASTER_NODE}$6"
REF_INDEX_IMAGE="$7"
INTERVAL_KILL_LIST=$(echo "${REF_TWOBIT}" | sed 's/.2bit$/.kill.intervals/')
KMER_KILL_LIST=$(echo "${REF_TWOBIT}" | sed 's/.2bit$/.kill.kmers/')
ALTS_KILL_LIST=$(echo "${REF_TWOBIT}" | sed 's/.2bit$/.kill.alts/')

# extract any extra arguments to StructuralVariationDiscoveryPipelineSpark
shift $(($# < 7 ? $# : 7))
SV_ARGS=${*:-${SV_ARGS:-""}}
# expand any local variables passed as strings (e.g. PROJECT_OUTPUT_DIR)
eval "SV_ARGS=\"${SV_ARGS}\""

NUM_WORKERS=$(gcloud dataproc clusters list --project ${PROJECT_NAME} \
                    --filter "clusterName = ${CLUSTER_NAME}" \
                    | awk 'NR==1 { for (i=1; i<=NF; i++) { f[$i] = i } }
                           { print $(f["WORKER_COUNT"]) + $(f["PREEMPTIBLE_WORKER_COUNT"]) }' \
                    | tail -1)
if [ -z "${NUM_WORKERS}" ]; then
    echo "Cluster \"${CLUSTER_NAME}\" not found"
    exit 1
fi
NUM_EXECUTORS=$((2 * ${NUM_WORKERS}))

GATK_SV_TOOL=${GATK_SV_TOOL:-"StructuralVariationDiscoveryPipelineSpark"}

case ${GATK_SV_TOOL} in
    "StructuralVariationDiscoveryPipelineSpark")
        TOOL_OPTIONS="\
            -I ${INPUT_BAM} \
            -O ${PROJECT_OUTPUT_DIR}/variants/ \
            -R ${REF_TWOBIT} \
            --aligner-index-image ${REF_INDEX_IMAGE} \
            --exclusion-intervals ${INTERVAL_KILL_LIST} \
            --kmers-to-ignore ${KMER_KILL_LIST} \
            --cross-contigs-to-ignore ${ALTS_KILL_LIST} \
            --breakpoint-intervals ${PROJECT_OUTPUT_DIR}/intervals \
            --high-coverage-intervals "${PROJECT_OUTPUT_DIR}/highCoverageIntervals.bed" \
            --fastq-dir ${PROJECT_OUTPUT_DIR}/fastq \
	    --read-metadata ${PROJECT_OUTPUT_DIR}/read-metadata.txt \
            --contig-sam-file ${PROJECT_OUTPUT_DIR}/assemblies.bam \
            --target-link-file ${PROJECT_OUTPUT_DIR}/target_links.bedpe \
            --exp-interpret"
        ;;
    "ExtractSVEvidenceSpark")
        TOOL_OPTIONS="\
            -I ${INPUT_BAM} \
            -O ${PROJECT_OUTPUT_DIR}/evidence \
            -R ${REF_TWOBIT} \
            --aligner-index-image ${REF_INDEX_IMAGE} \
            --exclusion-intervals ${INTERVAL_KILL_LIST} \
            --kmers-to-ignore ${KMER_KILL_LIST} \
            --cross-contigs-to-ignore ${ALTS_KILL_LIST} \
            --breakpoint-intervals ${PROJECT_OUTPUT_DIR}/intervals \
            --high-coverage-intervals "${PROJECT_OUTPUT_DIR}/highCoverageIntervals.bed" \
            --fastq-dir ${PROJECT_OUTPUT_DIR}/fastq \
            --target-link-file ${PROJECT_OUTPUT_DIR}/target_links.bedpe"
        ;;
    *)
        echo "Unknown tool: ${GATK_SV_TOOL}" 1>&2
        exit 1
        ;;
esac

"${GATK_DIR}/gatk" ${GATK_SV_TOOL} \
    ${TOOL_OPTIONS} ${SV_ARGS} \
    -- \
    --spark-runner GCS \
    --cluster "${CLUSTER_NAME}" \
    --project "${PROJECT_NAME}" \
    --num-executors ${NUM_EXECUTORS} \
    --driver-memory 30G \
    --executor-memory 30G \
    --conf spark.executor.memoryOverhead=5000 \
    --conf spark.network.timeout=600 \
    --conf spark.executor.heartbeatInterval=120
