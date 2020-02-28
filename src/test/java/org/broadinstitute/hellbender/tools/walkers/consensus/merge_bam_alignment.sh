if [ "$#" -ne 3 ]; then
  echo "this script must 3 arguments: <unaligned_bam> <aligned_bam> <base_name>"
  exit 1
fi

export unaligned_bam=$1
export aligned_bam=$2
export base_name=$3


export picard_jar="/Users/tsato/workspace/picard/build/libs/picard.jar"
export reference="/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta"
export bwa_version="0.7.17-r1188"
export bwa_commandline="bwa mem -K 100000000 -t 16 ${reference}"

java -jar -Xmx6144m ${picard_jar} MergeBamAlignment \
REFERENCE_SEQUENCE=${reference} \
UNMAPPED_BAM="${unaligned_bam}" \
ALIGNED_BAM=${aligned_bam} \
OUTPUT="${base_name}".merged.bam \
PROGRAM_RECORD_ID="bwamem" \
PROGRAM_GROUP_VERSION="${bwa_version}" \
PROGRAM_GROUP_COMMAND_LINE="${bwa_commandline}" \
PROGRAM_GROUP_NAME="bwamem" \
CREATE_INDEX=true \
TMP_DIR=.