package org.broadinstitute.hellbender.engine;


import org.broadinstitute.barclay.argparser.CommandLineParser;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.TestProgramGroup;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;
import org.testng.annotations.DataProvider;

import java.io.IOException;

public class TwoPassReadsWalkerUnitTest extends CommandLineProgramTest{

    //public class
    @CommandLineProgramProperties(
            summary = "Dummy that reads file and counts it twice",
            oneLineSummary = "empty class",
            programGroup = TestProgramGroup.class
    )
    private static class dummyTwoPassReadsWalker extends TwoPassReadWalker {
        public int firstPass = 0;
        public int secondPass = 0;
        boolean betweenTraversals = false;
        @Override
        protected void firstPassApply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {
            firstPass++;
        }
        @Override
        protected void secondPassApply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {
            secondPass++;
        }
        @Override
        protected void afterFirstPass() {
            betweenTraversals = true;
        }
    }

    @Test(dataProvider = "unsortedFiles")
    public void testDifferentFormatEquivalentBehavior(String file) throws IOException {
        final TwoPassReadsWalkerUnitTest.dummyTwoPassReadsWalker tool = new TwoPassReadsWalkerUnitTest.dummyTwoPassReadsWalker();

        final String[] args = {
                "-I", getTestDataDir()+ file,
                "-R", getTestDataDir()+ "/count_reads.fasta"
        };

        tool.instanceMain(args);

        Assert.assertEquals(tool.firstPass, 8);
        Assert.assertEquals(tool.secondPass, 8);
        Assert.assertTrue(tool.betweenTraversals);
    }

    @DataProvider(name = "unsortedFiles")
    public Object[][] makeExtensions() {
        return new Object[][] {{"/count_reads.bam"}, {"/count_reads.sam"}, {"/count_reads.cram"}};
    }

    @Test(dataProvider = "sortedFiles")
    public void testIntervalFiltering(String file) {
        final TwoPassReadsWalkerUnitTest.dummyTwoPassReadsWalker tool = new TwoPassReadsWalkerUnitTest.dummyTwoPassReadsWalker();

        final String[] args = {
                "-I", getTestDataDir()+ file,
                "-R", getTestDataDir()+ "/count_reads.fasta",
                "-L", "chr7:10-40"
        };

        tool.instanceMain(args);

        Assert.assertEquals(tool.firstPass, 5);
        Assert.assertEquals(tool.secondPass, 5);
        Assert.assertTrue(tool.betweenTraversals);
    }

    @DataProvider(name = "sortedFiles")
    public Object[][] makeSortedExtensions() {
        return new Object[][] {{"/count_reads_sorted.bam"}, {"/count_reads_sorted.cram"}};
    }
}
