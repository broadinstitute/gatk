package org.broadinstitute.hellbender.tools.spark.pathseq;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Created by markw on 5/9/17.
 */
public class PSPairedUnpairedSplitterSparkTest {


    @Test(groups = "spark")
    public void testReadPairing() {
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        final List<GATKRead> readList = new ArrayList<>(2);

        final SAMSequenceDictionary seq = new SAMSequenceDictionary();
        seq.addSequence(new SAMSequenceRecord("test_seq", 1000));
        final SAMFileHeader header = new SAMFileHeader(seq);

        final int numReadPairs = 300;
        final int numUnpairedReads = 205;
        final int numFormerlyPairedReads = 400;
        final List<String> pairedReadNames = new ArrayList<>(numReadPairs * 2);
        final List<String> unpairedReadNames = new ArrayList<>(numUnpairedReads + numFormerlyPairedReads);

        for (int i = 0; i < numReadPairs; i++) {
            final String readName = "paired_" + i;
            final List<GATKRead> readPair = ArtificialReadUtils.createPair(header, readName, 101, 1, 102, false, false);
            readList.addAll(readPair);
            pairedReadNames.add(readName);
            pairedReadNames.add(readName);
        }

        for (int i = 0; i < numUnpairedReads; i++) {
            final String readName = "unpaired_" + i;
            final GATKRead unpairedRead = ArtificialReadUtils.createRandomRead(101);
            unpairedRead.setName(readName);
            readList.add(unpairedRead);
            unpairedReadNames.add(readName);
        }

        for (int i = 0; i < numFormerlyPairedReads; i++) {
            final String readName = "formerly_paired_" + i;
            final List<GATKRead> readPair = ArtificialReadUtils.createPair(header, readName, 101, 1, 102, false, false);
            final GATKRead formerlyPairedRead = readPair.get(0);
            readList.add(formerlyPairedRead);
            unpairedReadNames.add(readName);
        }

        Collections.shuffle(readList);
        final int numPartitions = 3;
        final JavaRDD<GATKRead> readRDD = ctx.parallelize(readList, numPartitions);
        final int readsPerPartitionGuess = (numReadPairs + numUnpairedReads + numFormerlyPairedReads) / numPartitions;
        final PSPairedUnpairedSplitterSpark splitter = new PSPairedUnpairedSplitterSpark(readRDD, readsPerPartitionGuess, true);

        final List<GATKRead> resultPaired = splitter.getPairedReads().collect();
        final List<GATKRead> resultUnpaired = splitter.getUnpairedReads().collect();
        splitter.close();

        Assert.assertEquals(resultPaired.size(), pairedReadNames.size());
        for (final GATKRead read : resultPaired) {
            final String readName = read.getName();
            Assert.assertTrue(pairedReadNames.contains(readName));
            pairedReadNames.remove(readName);
        }

        Assert.assertEquals(resultUnpaired.size(), unpairedReadNames.size());
        for (final GATKRead read : resultUnpaired) {
            final String readName = read.getName();
            Assert.assertTrue(unpairedReadNames.contains(readName));
            unpairedReadNames.remove(readName);
        }
    }
}