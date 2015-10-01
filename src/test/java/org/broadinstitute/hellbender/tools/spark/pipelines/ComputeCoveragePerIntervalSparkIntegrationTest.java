package org.broadinstitute.hellbender.tools.spark.pipelines;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.exome.ExomeReadCounts;
import org.broadinstitute.hellbender.utils.read.SamAssertionUtils;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.ExomeToolsTestUtils;
import org.junit.Assert;
import org.testng.annotations.Test;

import java.io.File;

public final class ComputeCoveragePerIntervalSparkIntegrationTest extends CommandLineProgramTest {

    private File testFile(final String fileName) {
        return new File(getToolTestDataDir() +"/"+ fileName);
    }

    private final File NA12878_BAM = testFile("exome-read-counts-NA12878.bam");

    private final File INTERVALS_BED = testFile("exome-read-counts-intervals.bed");

    private final File SAMPLE_COUNT_EXPECTED_OUTPUT = testFile("exome-read-counts-sample.output");

    @Test
    public void test() throws Exception {
        final File unsortedBam =  NA12878_BAM;
        final File intervals =  INTERVALS_BED;
        final File expectedTxt =  SAMPLE_COUNT_EXPECTED_OUTPUT;
        final File outputTxt = createTempFile("na12878.coverage", ".txt");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--"+ StandardArgumentDefinitions.INPUT_LONG_NAME); args.add(unsortedBam.getCanonicalPath());
        args.add("--"+ ExomeReadCounts.EXOME_FILE_FULL_NAME); args.add(intervals.getCanonicalPath());
        args.add("--"+StandardArgumentDefinitions.OUTPUT_LONG_NAME); args.add(outputTxt.getCanonicalPath());

        this.runCommandLine(args.getArgsArray());

        final String readInActual = FileUtils.readFileToString(outputTxt.getAbsoluteFile());
        final String readInExpected = FileUtils.readFileToString(expectedTxt.getAbsoluteFile());

        Assert.assertEquals(readInExpected, readInActual);
    }

}