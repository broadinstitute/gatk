package org.broadinstitute.hellbender.tools.copynumber;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

public final class PreprocessIntervalsIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_SUB_DIR = toolsTestDir + "copynumber/";
    private static final File INTERVAL_LIST_FILE = new File(TEST_SUB_DIR, "preprocess-intervals-test.interval_list");
    private static final File REFERENCE_FILE = new File(b37_reference_20_21);

    @DataProvider(name = "intervalInputsFromCommandLine")
    public Object[][] testData() {
        // Test for separate intervals
        final int binLengthSeparateIntervalTest = 10_000;
        final int paddingLengthSeparateIntervalTest = 0;
        final List<Interval> inputIntervalsSeparateIntervalTest = Arrays.asList(
                new Interval("20", 3_000, 20_000),
                new Interval("20", 200, 1_900)
        );
        final List<Interval> expectedBinsSeparateIntervalTest = Arrays.asList(
                new Interval("20", 200, 1_900),
                new Interval("20", 3_000, 12_999),
                new Interval("20", 13_000, 20_000)
        );

        // Test for overlapping intervals
        final int binLengthOverlappingIntervalTest = 10_000;
        final int paddingLengthOverlappingIntervalTest = 500;
        final List<Interval> inputIntervalsOverlappingIntervalTest = Arrays.asList(
                new Interval("20", 3_000, 20_000),
                new Interval("20", 200, 2_100)
        );
        final List<Interval> expectedBinsOverlappingIntervalTest = Arrays.asList(
                new Interval("20", 1, 10_000),
                new Interval("20", 10_001, 20_000),
                new Interval("20", 20_001, 20_500)
        );

        // Test for intervals reaching the ends of the contigs
        final int binLengthEdgeIntervalTest = 10_000;
        final int paddingLengthEdgeIntervalTest = 500;
        final List<Interval> inputIntervalsEdgeIntervalTest = Arrays.asList(
                new Interval("20", 3_000, 20_000),
                new Interval("20", 63_025_220, 63_025_520)
        );
        final List<Interval> expectedBinsEdgeIntervalTest = Arrays.asList(
                new Interval("20", 2_500, 12_499),
                new Interval("20", 12_500, 20_500),
                new Interval("20", 63_024_720, 63_025_520)
        );

        // Test for whole chromosome
        final int binLengthWholeChromosomeTest = 10_000_000;
        final int paddingLengthWholeChromosomeTest = 500;
        final List<Interval> inputIntervalsWholeChromosomeTest = Arrays.asList(new Interval("20", 1, 63_025_520));
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
        final List<Interval> inputIntervalsWholeGenomeTest = Arrays.asList();
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
                {binLengthOverlappingIntervalTest, paddingLengthOverlappingIntervalTest, inputIntervalsOverlappingIntervalTest, expectedBinsOverlappingIntervalTest},
                {binLengthEdgeIntervalTest, paddingLengthEdgeIntervalTest, inputIntervalsEdgeIntervalTest, expectedBinsEdgeIntervalTest},
                {binLengthWholeChromosomeTest, paddingLengthWholeChromosomeTest, inputIntervalsWholeChromosomeTest, expectedBinsWholeChromosomeTest},
                {binLengthWholeGenomeTest, paddingLengthWholeGenomeTest, inputIntervalsWholeGenomeTest, expectedBinsWholeGenomeTest}
        };
    }

    // Test for interval inputs given as -L command line arguments
    @Test(dataProvider = "intervalInputsFromCommandLine")
    public void testCommandLine(final int binLength, final int paddingLength, final List<Interval> inputIntervals, final List<Interval> binsExpected) {
        final String[] outputFileName =  {"GATK-preprocess-intervals-test", ".tmp"};
        final File outputFile = createTempFile(outputFileName[0], outputFileName[1]);
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addReference(REFERENCE_FILE)
                .addArgument(PreprocessIntervals.BIN_LENGTH_SHORT_NAME, Integer.toString(binLength))
                .addArgument(PreprocessIntervals.PADDING_SHORT_NAME, Integer.toString(paddingLength))
                .addOutput(outputFile);
        inputIntervals.forEach(i -> argsBuilder.addArgument("L", i.getContig() + ":" + i.getStart() + "-" + i.getEnd()));
        runCommandLine(argsBuilder);

        final IntervalList binsResult = IntervalList.fromFile(outputFile);
        Assert.assertEquals(binsResult, binsExpected);
    }

    // Test for interval inputs read from a file
    @Test
    public void singleFileTest() {
        final int binLength = 10_000;
        final int paddingLength = 5_000;
        final File outputFile = createTempFile("preprocess-intervals-test", ".tsv");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addReference(REFERENCE_FILE)
                .addArgument(PreprocessIntervals.BIN_LENGTH_SHORT_NAME, Integer.toString(binLength))
                .addArgument(PreprocessIntervals.PADDING_SHORT_NAME, Integer.toString(paddingLength))
                .addArgument("L",  INTERVAL_LIST_FILE.getAbsolutePath())
                .addOutput(outputFile);
        runCommandLine(argsBuilder);
        final IntervalList binsResult = IntervalList.fromFile(outputFile);

        final IntervalList binsExpected = new IntervalList(binsResult.getHeader().getSequenceDictionary());
        binsExpected.add(new Interval("20", 1, 10_000));
        binsExpected.add(new Interval("20", 10_001, 16_000));
        binsExpected.add(new Interval("21", 1, 5_100));
        binsExpected.add(new Interval("21", 15_000, 24_999));
        binsExpected.add(new Interval("21", 25_000, 27_000));

        Assert.assertEquals(binsResult, binsExpected);
    }
}