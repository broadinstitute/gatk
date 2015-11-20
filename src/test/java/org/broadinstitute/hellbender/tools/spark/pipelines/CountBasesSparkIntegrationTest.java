package org.broadinstitute.hellbender.tools.spark.pipelines;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

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

    @Test(dataProvider = "countBases")
    public void countBases(final String fileName, final String referenceFileName, final long expectedCount) throws Exception {
        final File unsortedBam = new File(getTestDataDir(), fileName);
        final File outputTxt = createTempFile("count_bases", ".txt");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--" + StandardArgumentDefinitions.INPUT_LONG_NAME);
        args.add(unsortedBam.getCanonicalPath());
        args.add("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        args.add(outputTxt.getCanonicalPath());
        if (null != referenceFileName) {
            final File REF = new File(getTestDataDir(), referenceFileName);
            args.add("-R");
            args.add(REF.getAbsolutePath());
        }
        this.runCommandLine(args.getArgsArray());

        final String readIn = FileUtils.readFileToString(outputTxt.getAbsoluteFile());
        Assert.assertEquals((int) Integer.valueOf(readIn), expectedCount);
    }
}
