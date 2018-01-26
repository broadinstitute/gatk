package org.broadinstitute.hellbender.tools.copynumber;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SimpleCountCollection;
import org.broadinstitute.hellbender.utils.IntervalMergingRule;
import org.broadinstitute.hellbender.utils.IntervalSetRule;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.List;

/**
 * Integration test for {@link CollectFragmentCounts}.
 *
 * @author Andrey Smirnov &lt;asmirnov@broadinstitute.org&gt;
 */
public final class CollectFragmentCountsIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_SUB_DIR = new File(toolsTestDir, "copynumber");

    private static final File NA12878_BAM = new File(TEST_SUB_DIR, "collect-fragment-counts-NA12878.bam");
    private static final File NA12878_FRAGMENT_COUNTS_EXPECTED_OUTPUT = new File(TEST_SUB_DIR, "collect-fragment-counts-NA12878-expected.tsv");
    private static final File INTERVALS_FILE = new File(TEST_SUB_DIR, "collect-fragment-counts-test.interval_list");

    @DataProvider(name = "testData")
    public Object[][] testData() {
        return new Object[][] {
                {NA12878_BAM, NA12878_FRAGMENT_COUNTS_EXPECTED_OUTPUT}
        };
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testIntervalSetRule() {
        final File resultOutputFile = createTempFile("collect-fragment-counts-test", ".tsv");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addInput(NA12878_BAM)
                .addArgument(StandardArgumentDefinitions.INTERVALS_LONG_NAME, INTERVALS_FILE.getAbsolutePath())
                .addArgument(IntervalArgumentCollection.INTERVAL_SET_RULE_LONG_NAME, IntervalSetRule.INTERSECTION.toString())
                .addArgument(CollectFragmentCounts.FORMAT_LONG_NAME, CollectFragmentCounts.Format.TSV.toString())
                .addOutput(resultOutputFile);
        runCommandLine(argsBuilder);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testIntervalExclusionPadding() {
        final File resultOutputFile = createTempFile("collect-fragment-counts-test", ".tsv");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addInput(NA12878_BAM)
                .addArgument(StandardArgumentDefinitions.INTERVALS_LONG_NAME, INTERVALS_FILE.getAbsolutePath())
                .addArgument(IntervalArgumentCollection.INTERVAL_EXCLUSION_PADDING_LONG_NAME, "1")
                .addArgument(CollectFragmentCounts.FORMAT_LONG_NAME, CollectFragmentCounts.Format.TSV.toString())
                .addOutput(resultOutputFile);
        runCommandLine(argsBuilder);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testIntervalPadding() {
        final File resultOutputFile = createTempFile("collect-fragment-counts-test", ".tsv");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addInput(NA12878_BAM)
                .addArgument(StandardArgumentDefinitions.INTERVALS_LONG_NAME, INTERVALS_FILE.getAbsolutePath())
                .addArgument(IntervalArgumentCollection.INTERVAL_PADDING_LONG_NAME, "1")
                .addArgument(CollectFragmentCounts.FORMAT_LONG_NAME, CollectFragmentCounts.Format.TSV.toString())
                .addOutput(resultOutputFile);
        runCommandLine(argsBuilder);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testIntervalMergingRule() {
        final File resultOutputFile = createTempFile("collect-fragment-counts-test", ".tsv");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addInput(NA12878_BAM)
                .addArgument(StandardArgumentDefinitions.INTERVALS_LONG_NAME, INTERVALS_FILE.getAbsolutePath())
                .addArgument(IntervalArgumentCollection.INTERVAL_MERGING_RULE_LONG_NAME, IntervalMergingRule.ALL.toString())
                .addArgument(CollectFragmentCounts.FORMAT_LONG_NAME, CollectFragmentCounts.Format.TSV.toString())
                .addOutput(resultOutputFile);
        runCommandLine(argsBuilder);
    }

    @Test(dataProvider = "testData")
    public void testTSVOutput(final File inputBAMFile, final File expectedOutputFile) {
        final File resultOutputFile = createTempFile("collect-fragment-counts-test", ".tsv");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addInput(inputBAMFile)
                .addArgument(StandardArgumentDefinitions.INTERVALS_LONG_NAME, INTERVALS_FILE.getAbsolutePath())
                .addArgument(IntervalArgumentCollection.INTERVAL_MERGING_RULE_LONG_NAME, IntervalMergingRule.OVERLAPPING_ONLY.toString())
                .addArgument(CollectFragmentCounts.FORMAT_LONG_NAME, CollectFragmentCounts.Format.TSV.toString())
                .addOutput(resultOutputFile);
        runCommandLine(argsBuilder);
        final SimpleCountCollection expectedCounts = SimpleCountCollection.read(expectedOutputFile);
        final SimpleCountCollection resultCounts = SimpleCountCollection.read(resultOutputFile);
        Assert.assertEquals(expectedCounts, resultCounts);
    }

    @Test(dataProvider = "testData")
    public void testHDF5Output(final File inputBAMFile, final File expectedOutputFile) {
        final File resultOutputFile = createTempFile("collect-fragment-counts-test", ".hdf5");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addInput(inputBAMFile)
                .addArgument(StandardArgumentDefinitions.INTERVALS_LONG_NAME, INTERVALS_FILE.getAbsolutePath())
                .addArgument(IntervalArgumentCollection.INTERVAL_MERGING_RULE_LONG_NAME, IntervalMergingRule.OVERLAPPING_ONLY.toString())
                .addArgument(CollectFragmentCounts.FORMAT_LONG_NAME, CollectFragmentCounts.Format.HDF5.toString())
                .addOutput(resultOutputFile);
        runCommandLine(argsBuilder);
        final SimpleCountCollection expectedCounts = SimpleCountCollection.read(expectedOutputFile);
        final SimpleCountCollection resultCounts = SimpleCountCollection.read(resultOutputFile);
        Assert.assertEquals(expectedCounts, resultCounts);
    }

    @DataProvider(name = "artificialReadsData")
    public Object[][] artificialReadsTestData() {
        final SAMFileHeader samHeader = ArtificialReadUtils.createArtificialSamHeader();
        final List<GATKRead> readPairFirstInPairIsLeft = ArtificialReadUtils.createPair(samHeader,
                "firstReadPair", 101, 1000, 1300, true, false);
        final GATKRead firstInPairPositive = readPairFirstInPairIsLeft.get(0);

        final List<GATKRead> readPairFirstInPairIsRight = ArtificialReadUtils.createPair(samHeader,
                "secondReadPair", 101, 800, 1000, false, true);
        final GATKRead firstInPairNegative = readPairFirstInPairIsRight.get(1);

        return new Object[][] {
                {firstInPairPositive, 1200},
                {firstInPairNegative, 950}
        };
    }

    @Test(dataProvider = "artificialReadsData")
    public void testFragmentCenterComputation(final GATKRead read, final int expectedFragmentCenterPosition) {
        final SimpleInterval expectedFragmentCenter = new SimpleInterval(read.getContig(), expectedFragmentCenterPosition, expectedFragmentCenterPosition);
        final SimpleInterval resultFragmentCenter = CollectFragmentCounts.ReadOrientation.getFragmentCenter(read);
        Assert.assertEquals(expectedFragmentCenter, resultFragmentCenter);
    }
}