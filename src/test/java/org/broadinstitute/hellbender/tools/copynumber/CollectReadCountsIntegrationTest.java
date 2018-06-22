package org.broadinstitute.hellbender.tools.copynumber;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SimpleCountCollection;
import org.broadinstitute.hellbender.utils.IntervalMergingRule;
import org.broadinstitute.hellbender.utils.IntervalSetRule;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

/**
 * Integration test for {@link CollectReadCounts}.
 *
 * @author Andrey Smirnov &lt;asmirnov@broadinstitute.org&gt;
 */
public final class CollectReadCountsIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_SUB_DIR = new File(toolsTestDir, "copynumber");

    private static final File NA12878_BAM = new File(TEST_SUB_DIR, "collect-read-counts-NA12878.bam");
    private static final File NA12878_READ_COUNTS_EXPECTED_OUTPUT = new File(TEST_SUB_DIR, "collect-read-counts-NA12878-expected.tsv");
    private static final File INTERVALS_FILE = new File(TEST_SUB_DIR, "collect-read-counts-test.interval_list");

    @DataProvider(name = "testData")
    public Object[][] testData() {
        return new Object[][] {
                {NA12878_BAM, NA12878_READ_COUNTS_EXPECTED_OUTPUT}
        };
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testIntervalSetRule() {
        final File resultOutputFile = createTempFile("collect-read-counts-test", ".tsv");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addInput(NA12878_BAM)
                .addArgument(StandardArgumentDefinitions.INTERVALS_LONG_NAME, INTERVALS_FILE.getAbsolutePath())
                .addArgument(IntervalArgumentCollection.INTERVAL_SET_RULE_LONG_NAME, IntervalSetRule.INTERSECTION.toString())
                .addArgument(CollectReadCounts.FORMAT_LONG_NAME, CollectReadCounts.Format.TSV.toString())
                .addOutput(resultOutputFile);
        runCommandLine(argsBuilder);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testIntervalExclusionPadding() {
        final File resultOutputFile = createTempFile("collect-read-counts-test", ".tsv");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addInput(NA12878_BAM)
                .addArgument(StandardArgumentDefinitions.INTERVALS_LONG_NAME, INTERVALS_FILE.getAbsolutePath())
                .addArgument(IntervalArgumentCollection.INTERVAL_EXCLUSION_PADDING_LONG_NAME, "1")
                .addArgument(CollectReadCounts.FORMAT_LONG_NAME, CollectReadCounts.Format.TSV.toString())
                .addOutput(resultOutputFile);
        runCommandLine(argsBuilder);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testIntervalPadding() {
        final File resultOutputFile = createTempFile("collect-read-counts-test", ".tsv");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addInput(NA12878_BAM)
                .addArgument(StandardArgumentDefinitions.INTERVALS_LONG_NAME, INTERVALS_FILE.getAbsolutePath())
                .addArgument(IntervalArgumentCollection.INTERVAL_PADDING_LONG_NAME, "1")
                .addArgument(CollectReadCounts.FORMAT_LONG_NAME, CollectReadCounts.Format.TSV.toString())
                .addOutput(resultOutputFile);
        runCommandLine(argsBuilder);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testIntervalMergingRule() {
        final File resultOutputFile = createTempFile("collect-read-counts-test", ".tsv");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addInput(NA12878_BAM)
                .addArgument(StandardArgumentDefinitions.INTERVALS_LONG_NAME, INTERVALS_FILE.getAbsolutePath())
                .addArgument(IntervalArgumentCollection.INTERVAL_MERGING_RULE_LONG_NAME, IntervalMergingRule.ALL.toString())
                .addArgument(CollectReadCounts.FORMAT_LONG_NAME, CollectReadCounts.Format.TSV.toString())
                .addOutput(resultOutputFile);
        runCommandLine(argsBuilder);
    }

    @Test(dataProvider = "testData")
    public void testTSVOutput(final File inputBAMFile, final File expectedOutputFile) {
        final File resultOutputFile = createTempFile("collect-read-counts-test", ".tsv");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addInput(inputBAMFile)
                .addArgument(StandardArgumentDefinitions.INTERVALS_LONG_NAME, INTERVALS_FILE.getAbsolutePath())
                .addArgument(IntervalArgumentCollection.INTERVAL_MERGING_RULE_LONG_NAME, IntervalMergingRule.OVERLAPPING_ONLY.toString())
                .addArgument(CollectReadCounts.FORMAT_LONG_NAME, CollectReadCounts.Format.TSV.toString())
                .addOutput(resultOutputFile);
        runCommandLine(argsBuilder);
        final SimpleCountCollection expectedCounts = SimpleCountCollection.read(expectedOutputFile);
        final SimpleCountCollection resultCounts = SimpleCountCollection.read(resultOutputFile);
        Assert.assertEquals(expectedCounts, resultCounts);
    }

    @Test(dataProvider = "testData")
    public void testHDF5Output(final File inputBAMFile, final File expectedOutputFile) {
        final File resultOutputFile = createTempFile("collect-read-counts-test", ".hdf5");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addInput(inputBAMFile)
                .addArgument(StandardArgumentDefinitions.INTERVALS_LONG_NAME, INTERVALS_FILE.getAbsolutePath())
                .addArgument(IntervalArgumentCollection.INTERVAL_MERGING_RULE_LONG_NAME, IntervalMergingRule.OVERLAPPING_ONLY.toString())
                .addArgument(CollectReadCounts.FORMAT_LONG_NAME, CollectReadCounts.Format.HDF5.toString())
                .addOutput(resultOutputFile);
        runCommandLine(argsBuilder);
        final SimpleCountCollection expectedCounts = SimpleCountCollection.read(expectedOutputFile);
        final SimpleCountCollection resultCounts = SimpleCountCollection.read(resultOutputFile);
        Assert.assertEquals(expectedCounts, resultCounts);
    }
}