if [ "$#" -lt 3 ]; then
  echo "this script must 3 arguments: <grouped_bam> <base_name> <out_dir>"
  exit 1
fi

export grouped_bam=$1
export base_name=$2
export out_dir=$3
export debug=$4

export home_dir="/dsde/working/tsato/consensus/"
# export fgbio_jar="/dsde/working/tsato/consensus/fgbio.jar"
export fgbio_jar="/Users/tsato/Downloads/fgbio-1.1.0.jar"
# export tp53_dir="/dsde/working/tsato/consensus/tp53/"
# export bam="/dsde/working/tsato/consensus/tp53/Jonna_Grimsby_A04_denovo_bloodbiopsy_1pct_rep1.tp53.17_7578712.bam"
# export out="/dsde/working/tsato/consensus/tp53/Jonna_Grimsby_A04_denovo_bloodbiopsy_1pct_rep1.tp53.17_7578712.grouped.bam"

if [ "$debug" ]; then
  export java_command="java -jar -agentlib:jdwp=transport=dt_socket,server=y,suspend=n,address=49015 $fgbio_jar"
else
  export java_command="java -jar $fgbio_jar"
fi

$java_command \
CallDuplexConsensusReads \
-p $base_name \
-i $grouped_bam \
--min-reads 1 1 1 \
-o "$out_dir$base_name.consensus.bam"

samtools index "$out_dir$base_name.consensus.bam"