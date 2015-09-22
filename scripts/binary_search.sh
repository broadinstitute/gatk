#!/bin/sh


# 
# This script does a binary search to look for the smallest input that will
# break something.
#
# JPM

# The folder for us to store our temporary files. There may be many of them.
# This folder should exist before you start the script.
WORKDIR=tmp/search

# The code below assumes there's a "hb" script on the root that runs hellbender.
# Its contents would be:
#
# #!/bin/sh
# build/install/$HELLBENDER/bin/$HELLBENDER "$@"
#
# where you replace $HELLBENDER with the name of your hellbender folder
# (typically just "hellbender")

#
# Generates the input to be tested.
# It takes two arguments: the first and last element of the
# interval.
#
genInput () {
  INPUT=$WORKDIR/CEUTrio.test.ch20.$1-$2.bam
  STDINPUT=$WORKDIR/search_input.bam
  rm $STDINPUT
  # keep a copy of every input that we generated
  ./hb PrintReads -I src/test/resources/org/broadinstitute/hellbender/tools/BQSR/CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.bam -O $INPUT -L 20:$1-$2 2>&1 > $WORKDIR/inputgen.log || echo "input gen reports an error"
  # Standard filename for the next stage
  cp $INPUT $STDINPUT
  samtools index $STDINPUT
  CLOUDINPUT=gs://jpmartin/hellbender-test-inputs/search_input.bam
  gsutil cp $STDINPUT $CLOUDINPUT
  gsutil acl set public-read $CLOUDINPUT
  gsutil cp $STDINPUT.bai $CLOUDINPUT.bai
  gsutil acl set public-read $CLOUDINPUT.bai
}

#
# Generates the "known good" output.
#
genCorrectOutput () {
  STDINPUT=$WORKDIR/search_input.bam
  ./hb BaseRecalibrator -R ../../sample-data/GRCh37.ch20.fixup.fa -I $STDINPUT --knownSites tstBQSR/dbsnp_132.b37.excluding_sites_after_129.chr17_69k_70k.vcf -sortAllCols --RECAL_TABLE_FILE $WORKDIR/recaltable.correct 2>&1 > $WORKDIR/outputgen.log  || echo "correct output gen reports an error"
}

#
# Runs the computation under test.
#
genTestOutput () {
  CLOUDINPUT=gs://jpmartin/hellbender-test-inputs/search_input.bam
  ./hb BaseRecalibratorDataflow --runner=BLOCKING --numWorkers=1 -R gg://reference/EOSsjdnTicvzwAE -I $CLOUDINPUT --baseRecalibrationKnownVariants tstBQSR/dbsnp_132.b37.excluding_sites_after_129.chr17_69k_70k.vcf --apiKey=$GOOGLE_API_KEY --project=genomics-pipelines --staging=gs://jpmartin/staging-jpmartin -sortAllCols --RECAL_TABLE_FILE $WORKDIR/recaltable.cloud 2>&1 > $WORKDIR/testgen.log || echo "test output gen reports an error"
}

# This tries to shrink a failing interval by removing less than half.
# (It runs after the binary search, so we know that removing half is too much)
rogner() {
# scrunch the sides
# left side
while true; do
  MID=$(( $FIRST + ($LAST - $FIRST) / 5 ))
  if [ $FIRST -eq $MID ]; then
    echo "Done shrinking the left"
    break
  fi
  rm -f $WORKDIR/recaltable.correct $WORKDIR/recaltable.cloud
  genInput $MID $LAST
  genCorrectOutput
  genTestOutput
  if [ ! -f $WORKDIR/recaltable.correct ]; then
    echo "missing correct output for $MID-$LAST."
    break
  fi
  if [ ! -f $WORKDIR/recaltable.cloud ]; then
    echo "Output missing (test) for $MID-$LAST."
    FIRST=$MID;
    continue
  fi
  if ! `diff $WORKDIR/recaltable.correct $WORKDIR/recaltable.cloud`; then
    # they differ, explore this branch
    echo "Output is wrong for $MID-$LAST, going deeper"
    FIRST=$MID;
  else 
    # they are the same, try the other half
    echo "Output matches for $MID-$LAST, going to the right side"
    break
  fi
done
# right side
while true; do
  MID=$(( $LAST - ($LAST - $FIRST) / 5 ))
  if [ $LAST -eq $MID ]; then
    echo "Done shrinking the right"
    break
  fi
  rm -f $WORKDIR/recaltable.correct $WORKDIR/recaltable.cloud
  genInput $FIRST $MID
  genCorrectOutput
  genTestOutput
  if [ ! -f $WORKDIR/recaltable.correct ]; then
    echo "missing correct output for $FIRST-$MID."
    break
  fi
  if [ ! -f $WORKDIR/recaltable.cloud ]; then
    echo "Output missing (test) for $FIRST-$MID."
    LAST=$MID;
    continue
  fi
  if ! `diff $WORKDIR/recaltable.correct $WORKDIR/recaltable.cloud`; then
    # they differ, explore this branch
    echo "Output is wrong for $FIRST-$MID, going deeper"
    LAST=$MID;
  else 
    # they are the same, try the other half
    echo "Output matches for $FIRST-$MID, done"
    break
  fi
done

echo "Output is now $FIRST-$LAST"
}



binarySearch () {
while true; do

  MID=$(( ( $FIRST + $LAST ) / 2 ))
  if [ $FIRST -eq $MID ]; then
    echo "We have narrowed it down to just $FIRST-$LAST, we're done now."
    break
  fi
  rm -f $WORKDIR/recaltable.correct $WORKDIR/recaltable.cloud
  genInput $FIRST $MID
  genCorrectOutput
  genTestOutput
  if [ ! -f $WORKDIR/recaltable.correct ]; then
    echo "missing correct output for $FIRST-$MID."
    break
  fi
  if [ ! -f $WORKDIR/recaltable.cloud ]; then
    echo "missing test output for $FIRST-$MID."
  fi
  if ! `diff $WORKDIR/recaltable.correct $WORKDIR/recaltable.cloud`; then
    # they differ, explore this branch
    echo "Output is wrong for $FIRST-$MID, going deeper"
    LAST=$MID;
  else 
    # they are the same, try the other half
    echo "Output is correct for $FIRST-$MID, trying the latter half"
    rm -f $WORKDIR/recaltable.correct $WORKDIR/recaltable.cloud
    genInput $MID $LAST
    genCorrectOutput
    genTestOutput
    if [ ! -f $WORKDIR/recaltable.correct ]; then
      echo "missing correct output for $MID-$LAST."
      break
    fi
    if [ ! -f $WORKDIR/recaltable.cloud ]; then
      echo "missing test output for $MID-$LAST."
    fi
    if ! `diff $WORKDIR/recaltable.correct $WORKDIR/recaltable.cloud`; then
      # they differ, explore this branch
      echo "Output is wrong for $MID-$LAST, going deeper"
      FIRST=$MID;
    else
      echo "Output is correct for $MID-$LAST too!"
      # they are the same, too!
      break
    fi
  fi

done
}


#
# Set here the genomic interval you'd like to start from.
#
FIRST=1000250;
LAST=1000312;

binarySearch
rogner
