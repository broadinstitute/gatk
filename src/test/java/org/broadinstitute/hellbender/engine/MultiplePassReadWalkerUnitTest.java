package org.broadinstitute.hellbender.engine;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.TestProgramGroup;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;
import org.testng.annotations.DataProvider;

public class MultiplePassReadWalkerUnitTest extends CommandLineProgramTest{

    final private static String numberOfPassesArg = "number-of-passes";

    @CommandLineProgramProperties(
            summary = "Dummy that reads a file and counts the reads container n times",
            oneLineSummary = "empty class",
            programGroup = TestProgramGroup.class
    )
    private static class dummyMultiplePassReadWalker extends MultiplePassReadWalker {
        public int[] readCounts;
        boolean betweenTraversals = false;

        @Argument(fullName=numberOfPassesArg)
        public int numberOfPasses;

        @Override
        public void onTraversalStart() {
            readCounts = new int[numberOfPasses];
        }

        @Override
        public void traverseReads() {
            for (int i = 0; i < numberOfPasses; i++) {
                final int passNumber = i;
                forEachRead((GATKRead read, ReferenceContext reference, FeatureContext features) -> readCounts[passNumber]++);
                betweenTraversals = true;
            }
        }
    }

    @Test(dataProvider = "unsortedFiles")
    public void testTwoPass(String file) {
        final MultiplePassReadWalkerUnitTest.dummyMultiplePassReadWalker tool = new MultiplePassReadWalkerUnitTest.dummyMultiplePassReadWalker();

        final String[] args = {
                "-I", getTestDataDir()+ file,
                "-R", getTestDataDir()+ "/count_reads.fasta",
                "--" + numberOfPassesArg, "2"
        };

        tool.instanceMain(args);

        Assert.assertEquals(tool.readCounts[0], 8);
        Assert.assertEquals(tool.readCounts[1], 8);
        Assert.assertTrue(tool.betweenTraversals);
    }

    @DataProvider(name="passCounts")
    public Object[] getPassCounts() {
        return new Object[] {
                0, 1, 2, 5, 25
        };
    }

    @Test(dataProvider="passCounts")
    public void testPassCounts(final Integer numberOfPasses) {
        final MultiplePassReadWalkerUnitTest.dummyMultiplePassReadWalker tool = new MultiplePassReadWalkerUnitTest.dummyMultiplePassReadWalker();

        final String[] args = {
                "-I", getTestDataDir()+ "/count_reads.bam",
                "--" + numberOfPassesArg, numberOfPasses.toString()
        };

        tool.instanceMain(args);

        for (int i = 0; i < numberOfPasses; i++) {
            Assert.assertEquals(tool.readCounts[i], 8);
        }

        Assert.assertEquals(numberOfPasses == 0 ? false : true, tool.betweenTraversals);
    }

    @DataProvider(name = "unsortedFiles")
    public Object[][] getUnsortedFiles() {
        return new Object[][] {
                {"/count_reads.bam"},
                {"/count_reads.sam"},
                {"/count_reads.cram"}
        };
    }

    @Test(dataProvider = "sortedFiles")
    public void testIntervalFiltering(String file) {
        final MultiplePassReadWalkerUnitTest.dummyMultiplePassReadWalker tool = new MultiplePassReadWalkerUnitTest.dummyMultiplePassReadWalker();

        final String[] args = {
                "-I", getTestDataDir()+ file,
                "-R", getTestDataDir()+ "/count_reads.fasta",
                "-L", "chr7:10-40",
                "--" + numberOfPassesArg, "2"
        };

        tool.instanceMain(args);

        Assert.assertEquals(tool.readCounts[0], 5);
        Assert.assertEquals(tool.readCounts[1], 5);
        Assert.assertTrue(tool.betweenTraversals);
    }

    @DataProvider(name = "sortedFiles")
    public Object[][] getSortedFiles() {
        return new Object[][] {
                {"/count_reads_sorted.bam"},
                {"/count_reads_sorted.cram"}};
    }

    @Test(dataProvider = "unsortedFiles", expectedExceptions = UserException.class)
    public void testRejectAttemptToIndexNonIndexableInput(String file) {
        final MultiplePassReadWalkerUnitTest.dummyMultiplePassReadWalker tool = new MultiplePassReadWalkerUnitTest.dummyMultiplePassReadWalker();

        final String[] args = {
                "-I", getTestDataDir()+ file,
                "-R", getTestDataDir()+ "/count_reads.fasta",
                "-L", "chr7:10-40",
                "--" + numberOfPassesArg, "2"
        };

        tool.instanceMain(args);
    }

}
