#
# This is not efficient at all... but gets the job done
#
echo -e "index\tcontig\tmin\tmax\ttotal_bases"
for f in ${1}/*.interval_list; do    
    SHARD=`basename $f | cut -d"-" -f1`
    START=`cat $f | grep -v "@" | head -1 | cut -f1-2`
    END=`cat $f| grep -v "@" | tail -1 | cut -f3`
    BASES=`cat $f | grep -v "@" | awk '{sum+=($3-$2);} END{print sum;}'`
    echo -e "$SHARD\t$START\t$END\t$BASES"
done