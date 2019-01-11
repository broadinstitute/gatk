package org.broadinstitute.hellbender.tools.spark.pipelines;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
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

    @Test(groups = "spark")
    public void testNoNPRWhenOutputIsUnspecified(){
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addInput(new File(getTestDataDir(), "count_bases.bam"));
        this.runCommandLine(args.getArgsArray());
    }
}
