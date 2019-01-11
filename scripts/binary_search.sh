#!/bin/sh


# 
# This script does a binary search to look for the smallest input that will
# break something.
#

# The folder for us to store our temporary files. There may be many of them.
WORKDIR=tmp/search

# This setting determines where we save the input we generate.
# You should change it to something you have write access to.
CLOUDINPUT=gs://jpmartin/hellbender-test-inputs/search_input.bam

# file containing the known good output
GOODOUTPUT=$WORKDIR/recaltable.correct

# file containing the output under test.
TESTOUTPUT=$WORKDIR/recaltable.cloud

#
# Set here the genomic interval you'd like to start from.
#
FIRST=1000250;
LAST=1000312;


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
# The output will be saved to $CLOUDINPUT
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
  # copy to $CLOUDINPUT so the next phases can get to it
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
  # change this line to run the known good version of the computation you're interested in.
  ./hb BaseRecalibrator -R ../../sample-data/GRCh37.ch20.fixup.fa -I $STDINPUT --knownSites tstBQSR/dbsnp_132.b37.excluding_sites_after_129.chr17_69k_70k.vcf -sortAllCols --RECAL_TABLE_FILE $GOODOUTPUT 2>&1 > $WORKDIR/outputgen.log  || echo "correct output gen reports an error"
}

#
# Runs the computation under test.
#
genTestOutput () {
  # change this line to run the computation you're testing.
  ./hb BaseRecalibratorDataflow --runner=BLOCKING --numWorkers=1 -R gg://reference/EOSsjdnTicvzwAE -I $CLOUDINPUT --baseRecalibrationKnownVariants tstBQSR/dbsnp_132.b37.excluding_sites_after_129.chr17_69k_70k.vcf --project=genomics-pipelines --staging=gs://jpmartin/staging-jpmartin -sortAllCols --RECAL_TABLE_FILE $TESTOUTPUT 2>&1 > $WORKDIR/testgen.log || echo "test output gen reports an error"
}

# This tries to shrink a failing interval by removing less than half.
# (It runs after the binary search, so we know that removing half is too much)
rogner() {
# try to move the left side by 1/5th at a time
while true; do
  MID=$(( $FIRST + ($LAST - $FIRST) / 5 ))
  if [ $FIRST -eq $MID ]; then
    echo "Done shrinking the left"
    break
  fi
  rm -f $GOODOUTPUT $TESTOUTPUT
  genInput $MID $LAST
  genCorrectOutput
  genTestOutput
  if [ ! -f $GOODOUTPUT ]; then
    echo "missing correct output for $MID-$LAST."
    break
  fi
  if [ ! -f $TESTOUTPUT ]; then
    echo "Output missing (test) for $MID-$LAST."
    FIRST=$MID;
    continue
  fi
  if ! `diff $WORKDIR/recaltable.correct $TESTOUTPUT`; then
    # they differ, explore this branch
    echo "Output is wrong for $MID-$LAST, going deeper"
    FIRST=$MID;
  else 
    # they are the same, try the other half
    echo "Output matches for $MID-$LAST, going to the right side"
    break
  fi
done
# right side now, again 1/5th at a time
while true; do
  MID=$(( $LAST - ($LAST - $FIRST) / 5 ))
  if [ $LAST -eq $MID ]; then
    echo "Done shrinking the right"
    break
  fi
  rm -f $GOODOUTPUT $TESTOUTPUT
  genInput $FIRST $MID
  genCorrectOutput
  genTestOutput
  if [ ! -f $GOODOUTPUT ]; then
    echo "missing correct output for $FIRST-$MID."
    break
  fi
  if [ ! -f $TESTOUTPUT ]; then
    echo "Output missing (test) for $FIRST-$MID."
    LAST=$MID;
    continue
  fi
  if ! `diff $GOODOUTPUT $TESTOUTPUT`; then
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



#
# This runs a binary search to look for the smallest interval there the output
# of the test program differs from the output of the good program.
#
binarySearch () {
while true; do

  MID=$(( ( $FIRST + $LAST ) / 2 ))
  if [ $FIRST -eq $MID ]; then
    echo "We have narrowed it down to just $FIRST-$LAST, we're done now."
    break
  fi
  rm -f $GOODOUTPUT $TESTOUTPUT
  genInput $FIRST $MID
  genCorrectOutput
  genTestOutput
  if [ ! -f $GOODOUTPUT ]; then
    echo "missing correct output for $FIRST-$MID."
    break
  fi
  if [ ! -f $TESTOUTPUT ]; then
    echo "missing test output for $FIRST-$MID."
  fi
  if ! `diff $GOODOUTPUT $TESTOUTPUT`; then
    # they differ, explore this branch
    echo "Output is wrong for $FIRST-$MID, going deeper"
    LAST=$MID;
  else 
    # they are the same, try the other half
    echo "Output is correct for $FIRST-$MID, trying the latter half"
    rm -f $GOODOUTPUT $TESTOUTPUT
    genInput $MID $LAST
    genCorrectOutput
    genTestOutput
    if [ ! -f $GOODOUTPUT ]; then
      echo "missing correct output for $MID-$LAST."
      break
    fi
    if [ ! -f $TESTOUTPUT ]; then
      echo "missing test output for $MID-$LAST."
    fi
    if ! `diff $GOODOUTPUT $TESTOUTPUT`; then
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


# create the work directory if it's not already there
[ -e $WORKDIR ] || mkdir -p $WORKDIR

# shrink the range, halving every time
binarySearch
# if there is still something left, try removing 1/5th at a time
rogner
