package org.broadinstitute.hellbender.engine;


import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.CommandLineParser;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.TestProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkToolIntegrationTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;
import org.testng.annotations.DataProvider;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class TwoPassReadsWalkerUnitTest extends CommandLineProgramTest{

    //public class
    @CommandLineProgramProperties(
            summary = "Dummy that reads file and counts it twice .",
            oneLineSummary = "empty class",
            programGroup = TestProgramGroup.class
    )
    private static class dummyTwoPassReadsWalker extends TwoPassReadWalker {
        int firstPass = 0;
        int secondPass = 0;
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

        public boolean everythingWorked() {
            if (firstPass != secondPass || !betweenTraversals) {
                return false;
            }
            return true;
        }
    }

    @Override
    public String getTestedClassName() {
        return TwoPassReadsWalkerUnitTest.dummyTwoPassReadsWalker.class.getSimpleName();
    }

    @Test(dataProvider = "extensions")
    public void testDifferentFormatEquivalentBehavior(String ext) throws IOException {
        final TwoPassReadsWalkerUnitTest.dummyTwoPassReadsWalker tool = new TwoPassReadsWalkerUnitTest.dummyTwoPassReadsWalker();
        final CommandLineParser clp = new CommandLineParser(tool);

        final String[] args = {
                "-I", getTestDataDir()+ "/count_reads" + ext,
                "-R", getTestDataDir()+ "/count_reads.fasta"
        };


        clp.parseArguments(System.out, args);
        tool.onStartup();
        tool.doWork();
        tool.onShutdown();

        Assert.assertTrue(tool.everythingWorked());
    }

    @DataProvider(name = "extensions")
    public Object[][] makeExtensions() {
        return new Object[][] {{".bam"}, {".sam"}, {".cram"}};
    }
}
