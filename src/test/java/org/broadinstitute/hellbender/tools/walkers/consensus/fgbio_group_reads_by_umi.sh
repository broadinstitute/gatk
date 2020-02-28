if [ "$#" -lt 3 ]; then
  echo "this script takes 3 arguments: <bam> <out> <out_dir> (<debug>)"
  exit 1
fi

# (2/27/20): If you get a segfault from groupByUMI--

export bam=$1
export out=$2
export out_dir=$3
export debug=$4
export min_mapq=28

echo $bam
echo $out


export fgbio_jar="/Users/tsato/Downloads/fgbio-1.1.0.jar"

#export home_dir="/dsde/working/tsato/consensus/"

#export tp53_dir="/dsde/working/tsato/consensus/tp53/"
# export bam="/dsde/working/tsato/consensus/tp53/Jonna_Grimsby_A04_denovo_bloodbiopsy_1pct_rep1.tp53.17_7578712.bam"
# export out="/dsde/working/tsato/consensus/tp53/Jonna_Grimsby_A04_denovo_bloodbiopsy_1pct_rep1.tp53.17_7578712.grouped.bam"

# export bam="/dsde/working/tsato/consensus/tp53/Jonna_Grimsby_A04_denovo_bloodbiopsy_1pct_rep1.tp53.tiny.bam"
# export bam="/dsde/working/tsato/consensus/tp53/Jonna_Grimsby_A04_denovo_bloodbiopsy_1pct_rep1.tp53.bam"
# export out="/dsde/working/tsato/consensus/tp53/Jonna_Grimsby_A04_denovo_bloodbiopsy_1pct_rep1.tp53.grouped.bam"


if [ "$debug" ]; then
  export java_command="java -jar -agentlib:jdwp=transport=dt_socket,server=y,suspend=n,address=49015 $fgbio_jar"
else
  export java_command="java -jar $fgbio_jar"
fi

$java_command \
GroupReadsByUmi \
-e 1 \
-i $bam \
-o $out \
--strategy=paired \
--min-map-q="$min_mapq" \
-f "$out_dir"/family_size.tsv

# create an indel only sam file
# it contains some 146M's but that's OK
# export out_indel_only="/dsde/working/tsato/consensus/tp53/Jonna_Grimsby_A04_denovo_bloodbiopsy_1pct_rep1.tp53.grouped.I.sam"
# samtools view -h $out | egrep "^@|[0-9]+I[0-9]+M" > $out_indel_only


# fgbio fails after print reads
# need to create a good indel set....
# but let's really study what is happening with fgbio - what it's getting rid of etc.