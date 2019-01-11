package org.broadinstitute.hellbender.tools.copynumber;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.utils.IntervalMergingRule;
import org.broadinstitute.hellbender.utils.IntervalSetRule;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public final class PreprocessIntervalsIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_SUB_DIR = new File(toolsTestDir, "copynumber");
    private static final File INTERVAL_LIST_FILE = new File(TEST_SUB_DIR, "preprocess-intervals-test.interval_list");
    private static final File REFERENCE_FILE = new File(b37_reference_20_21);

    @DataProvider(name = "intervalInputsFromCommandLine")
    public Object[][] testData() {
        // Test for separate intervals
        final int binLengthSeparateIntervalTest = 10_000;
        final int paddingLengthSeparateIntervalTest = 0;
        final List<Interval> inputIntervalsSeparateIntervalTest = Arrays.asList(
                new Interval("20", 103_000, 120_000),
                new Interval("20", 100_200, 101_900)
        );
        final List<Interval> expectedBinsSeparateIntervalTest = Arrays.asList(
                new Interval("20", 100_200, 101_900),
                new Interval("20", 103_000, 112_999),
                new Interval("20", 113_000, 120_000)
        );

        // Test for no binning (specified by zero bin length)
        final int binLengthNoBinningTest = 0;
        final int paddingLengthNoBinningTest = 0;
        final List<Interval> inputIntervalsNoBinningTest = Arrays.asList(
                new Interval("20", 103_000, 120_000),
                new Interval("20", 100_200, 101_900)
        );
        final List<Interval> expectedBinsNoBinningTest = Arrays.asList(
                new Interval("20", 100_200, 101_900),
                new Interval("20", 103_000, 120_000)
        );

        // Test for overlapping intervals
        final int binLengthOverlappingIntervalTest = 10_000;
        final int paddingLengthOverlappingIntervalTest = 500;
        final List<Interval> inputIntervalsOverlappingIntervalTest = Arrays.asList(
                new Interval("20", 103_001, 120_000),
                new Interval("20", 100_201, 102_100)
        );
        final List<Interval> expectedBinsOverlappingIntervalTest = Arrays.asList(
                new Interval("20", 99_701, 102_550),
                new Interval("20", 102_551, 112_550),
                new Interval("20", 112_551, 120_500)
        );

        // Test for intervals reaching the ends of the contigs
        final int binLengthEdgeIntervalTest = 100_000;
        final int paddingLengthEdgeIntervalTest = 5_000;
        final List<Interval> inputIntervalsEdgeIntervalTest = Arrays.asList(
                new Interval("20", 3_000, 200_000),
                new Interval("20", 62_935_000, 63_025_520)
        );
        final List<Interval> expectedBinsEdgeIntervalTest = Arrays.asList(
                new Interval("20", 1, 100_000),
                new Interval("20", 100_001, 200_000),
                new Interval("20", 200_001, 205_000),
                new Interval("20", 62_930_000, 63_025_520)
        );

        // Test for dropping of intervals with all Ns (first 60_000 bases in test reference are Ns)
        final int binLengthDropNsIntervalTest = 60_000;
        final int paddingLengthDropNsIntervalTest = 0;
        final List<Interval> inputIntervalsDropNsIntervalTest = Collections.singletonList(
                new Interval("20", 1, 120_000));
        final List<Interval> expectedBinsDropNsIntervalTest = Collections.singletonList(
                new Interval("20", 60_001, 120_000));

        // Test for whole chromosome
        final int binLengthWholeChromosomeTest = 10_000_000;
        final int paddingLengthWholeChromosomeTest = 500;
        final List<Interval> inputIntervalsWholeChromosomeTest = Collections.singletonList(
                new Interval("20", 1, 63_025_520));
        final List<Interval> expectedBinsWholeChromosomeTest = Arrays.asList(
                new Interval("20", 1, 10_000_000),
                new Interval("20", 10_000_001, 20_000_000),
                new Interval("20", 20_000_001, 30_000_000),
                new Interval("20", 30_000_001, 40_000_000),
                new Interval("20", 40_000_001, 50_000_000),
                new Interval("20", 50_000_001, 60_000_000),
                new Interval("20", 60_000_001, 63_025_520)
        );

        // Test for whole genome -- when we don't give any intervals, then the tool assumes that the user wants to sequence the whole genome
        final int binLengthWholeGenomeTest = 10_000_000;
        final int paddingLengthWholeGenomeTest = 500;
        final List<Interval> inputIntervalsWholeGenomeTest = Collections.emptyList();
        final List<Interval> expectedBinsWholeGenomeTest = Arrays.asList(
                new Interval("20", 1, 10_000_000),
                new Interval("20", 10_000_001, 20_000_000),
                new Interval("20", 20_000_001, 30_000_000),
                new Interval("20", 30_000_001, 40_000_000),
                new Interval("20", 40_000_001, 50_000_000),
                new Interval("20", 50_000_001, 60_000_000),
                new Interval("20", 60_000_001, 63_025_520),
                new Interval("21", 1, 10_000_000),
                new Interval("21", 10_000_001, 20_000_000),
                new Interval("21", 20_000_001, 30_000_000),
                new Interval("21", 30_000_001, 40_000_000),
                new Interval("21", 40_000_001, 48_129_895)
        );

        // Return all test data
        return new Object[][]{
                {binLengthSeparateIntervalTest, paddingLengthSeparateIntervalTest, inputIntervalsSeparateIntervalTest, expectedBinsSeparateIntervalTest},
                {binLengthNoBinningTest, paddingLengthNoBinningTest, inputIntervalsNoBinningTest, expectedBinsNoBinningTest},
                {binLengthOverlappingIntervalTest, paddingLengthOverlappingIntervalTest, inputIntervalsOverlappingIntervalTest, expectedBinsOverlappingIntervalTest},
                {binLengthEdgeIntervalTest, paddingLengthEdgeIntervalTest, inputIntervalsEdgeIntervalTest, expectedBinsEdgeIntervalTest},
                {binLengthDropNsIntervalTest, paddingLengthDropNsIntervalTest, inputIntervalsDropNsIntervalTest, expectedBinsDropNsIntervalTest},
                {binLengthWholeChromosomeTest, paddingLengthWholeChromosomeTest, inputIntervalsWholeChromosomeTest, expectedBinsWholeChromosomeTest},
                {binLengthWholeGenomeTest, paddingLengthWholeGenomeTest, inputIntervalsWholeGenomeTest, expectedBinsWholeGenomeTest}
        };
    }

    // Test for interval inputs given as -L command line arguments
    @Test(dataProvider = "intervalInputsFromCommandLine")
    public void testCommandLine(final int binLength, final int paddingLength, final List<Interval> inputIntervals, final List<Interval> binsExpected) {
        final File outputFile = createTempFile("GATK-preprocess-intervals-test", ".interval_list");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addReference(REFERENCE_FILE)
                .addArgument(PreprocessIntervals.BIN_LENGTH_LONG_NAME, Integer.toString(binLength))
                .addArgument(PreprocessIntervals.PADDING_LONG_NAME, Integer.toString(paddingLength))
                .addArgument(IntervalArgumentCollection.INTERVAL_MERGING_RULE_LONG_NAME, IntervalMergingRule.OVERLAPPING_ONLY.toString())
                .addOutput(outputFile);
        inputIntervals.forEach(i -> argsBuilder.addArgument(StandardArgumentDefinitions.INTERVALS_LONG_NAME, i.getContig() + ":" + i.getStart() + "-" + i.getEnd()));
        runCommandLine(argsBuilder);

        final IntervalList binsResult = IntervalList.fromFile(outputFile);
        Assert.assertEquals(binsResult, binsExpected);
    }

    // Test for interval inputs read from a file
    @Test
    public void singleFileTest() {
        final int binLength = 10_000;
        final int paddingLength = 5_000;
        final File outputFile = createTempFile("preprocess-intervals-test", ".interval_list");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addReference(REFERENCE_FILE)
                .addArgument(PreprocessIntervals.BIN_LENGTH_LONG_NAME, Integer.toString(binLength))
                .addArgument(PreprocessIntervals.PADDING_LONG_NAME, Integer.toString(paddingLength))
                .addArgument(StandardArgumentDefinitions.INTERVALS_LONG_NAME, INTERVAL_LIST_FILE.getAbsolutePath())
                .addArgument(IntervalArgumentCollection.INTERVAL_MERGING_RULE_LONG_NAME, IntervalMergingRule.OVERLAPPING_ONLY.toString())
                .addOutput(outputFile);
        runCommandLine(argsBuilder);
        final IntervalList binsResult = IntervalList.fromFile(outputFile);

        final IntervalList binsExpected = new IntervalList(binsResult.getHeader().getSequenceDictionary());
        binsExpected.add(new Interval("20", 95_001, 105_000));
        binsExpected.add(new Interval("20", 105_001, 115_000));
        binsExpected.add(new Interval("20", 115_001, 125_000));

        Assert.assertEquals(binsResult, binsExpected);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testIntervalSetRule() {
        final File resultOutputFile = createTempFile("preprocess-intervals-test", ".interval_list");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addReference(REFERENCE_FILE)
                .addArgument(StandardArgumentDefinitions.INTERVALS_LONG_NAME, INTERVAL_LIST_FILE.getAbsolutePath())
                .addArgument(IntervalArgumentCollection.INTERVAL_SET_RULE_LONG_NAME, IntervalSetRule.INTERSECTION.toString())
                .addOutput(resultOutputFile);
        runCommandLine(argsBuilder);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testIntervalExclusionPadding() {
        final File resultOutputFile = createTempFile("preprocess-intervals-test", ".interval_list");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addReference(REFERENCE_FILE)
                .addArgument(StandardArgumentDefinitions.INTERVALS_LONG_NAME, INTERVAL_LIST_FILE.getAbsolutePath())
                .addArgument(IntervalArgumentCollection.INTERVAL_EXCLUSION_PADDING_LONG_NAME,"1")
                .addOutput(resultOutputFile);
        runCommandLine(argsBuilder);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testIntervalPadding() {
        final File resultOutputFile = createTempFile("preprocess-intervals-test", ".interval_list");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addReference(REFERENCE_FILE)
                .addArgument(StandardArgumentDefinitions.INTERVALS_LONG_NAME, INTERVAL_LIST_FILE.getAbsolutePath())
                .addArgument(IntervalArgumentCollection.INTERVAL_PADDING_LONG_NAME,"1")
                .addOutput(resultOutputFile);
        runCommandLine(argsBuilder);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testIntervalMergingRule() {
        final File resultOutputFile = createTempFile("preprocess-intervals-test", ".interval_list");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addReference(REFERENCE_FILE)
                .addArgument(StandardArgumentDefinitions.INTERVALS_LONG_NAME, INTERVAL_LIST_FILE.getAbsolutePath())
                .addArgument(IntervalArgumentCollection.INTERVAL_MERGING_RULE_LONG_NAME, IntervalMergingRule.ALL.toString())
                .addOutput(resultOutputFile);
        runCommandLine(argsBuilder);
    }
}