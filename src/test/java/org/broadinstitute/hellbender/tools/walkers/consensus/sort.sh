if [ "$#" -ne 3 ]; then
  echo "this script must 2 arguments: <unsorted_consensus_bam> <out_dir> <prefix>"
  exit 1
fi

export unsorted_consensus_bam=$1
export out_dir=$2
export prefix=$3


samtools sort $unsorted_consensus_bam -o ${out_dir}.${prefix}.sorted.bam
samtools index ${out_dir}.${prefix}.sorted.bam ${out_dir}.${prefix}.sorted.bai