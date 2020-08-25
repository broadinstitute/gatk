package org.broadinstitute.hellbender.tools.spark.pipelines;

import htsjdk.samtools.util.Locatable;
import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;

public final class CountBasesSparkIntegrationTest extends CommandLineProgramTest {

    @Override
    public String getTestedClassName() {
        return CountBasesSpark.class.getSimpleName();
    }


    @DataProvider(name="countBases")
    public Object[][] filenames() {
        return new Object[][]{
                {"count_bases.sam", null, 808L},
                {"count_bases.bam", null, 808L},
                {"count_bases.cram", "count_bases.fasta", 808L}
        };
    }

    @Test(dataProvider = "countBases", groups = "spark")
    public void countBases(final String fileName, final String referenceFileName, final long expectedCount) throws Exception {
        final File unsortedBam = new File(getTestDataDir(), fileName);
        final File outputTxt = createTempFile("count_bases", ".txt");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addInput(unsortedBam);
        args.addOutput(outputTxt);
        if (null != referenceFileName) {
            final File ref = new File(getTestDataDir(), referenceFileName);
            args.addReference(ref);
        }
        this.runCommandLine(args.getArgsArray());

        final String readIn = FileUtils.readFileToString(outputTxt.getAbsoluteFile(), StandardCharsets.UTF_8);
        Assert.assertEquals((int) Integer.valueOf(readIn), expectedCount);
    }

    @DataProvider
    public Object[][] getIntervals(){
        return new Object[][]{
                // The first case triggered https://github.com/broadinstitute/gatk/issues/6319 before the fix
                {new SimpleInterval("chr8"), 0L},
                {new SimpleInterval("chr7"), 707L},
                {new SimpleInterval("chr7:1-2"), 303L}
        };
    }

    @Test(groups = "spark", dataProvider = "getIntervals")
    public void testIntervals(Locatable interval, long expectedCount) throws IOException {
        final File sortedBam = new File(getTestDataDir(), "count_bases_sorted.bam");
        final File outputTxt = createTempFile("count_bases", ".txt");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addInput(sortedBam);
        args.addOutput(outputTxt);
        args.addInterval(interval);
        this.runCommandLine(args);

        final String readIn = FileUtils.readFileToString(outputTxt.getAbsoluteFile(), StandardCharsets.UTF_8);
        Assert.assertEquals((int) Integer.valueOf(readIn), expectedCount);
    }

    @Test(groups = "spark")
    public void testNoNPRWhenOutputIsUnspecified(){
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addInput(new File(getTestDataDir(), "count_bases.bam"));
        this.runCommandLine(args.getArgsArray());
    }
}
