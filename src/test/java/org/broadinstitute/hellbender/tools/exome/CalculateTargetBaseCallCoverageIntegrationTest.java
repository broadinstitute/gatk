package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.*;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Integration tests for {@link CalculateTargetBaseCallCoverage}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class CalculateTargetBaseCallCoverageIntegrationTest extends CommandLineProgramTest {

    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "exome");

    private static final File TEST_BAM_NA12872 = new File(TEST_DATA_DIR, "exome-read-counts-NA12872.bam");
    private static final File TEST_BAM_NA12778 = new File(TEST_DATA_DIR, "exome-read-counts-NA12778.bam");
    private static final File TEST_BAM_NA12878 = new File(TEST_DATA_DIR, "exome-read-counts-NA12878.bam");
    private static final File TARGETS_FILE = new File(TEST_DATA_DIR, "exome-read-counts-test-targets.tsv");
    private static final File EXPECTED_COUNTS_FILE = new File(TEST_DATA_DIR, "exome-read-counts.output");
    private static final File EXPECTED_COUNTS_MAX_OF_9_FILE = new File(TEST_DATA_DIR, "exome-read-counts-max-of-9.output");
    private static final File EXPECTED_COUNTS_MIN_MQ_30_FILE = new File(TEST_DATA_DIR, "exome-read-counts-min-MQ-30.output");
    private static final File EXPECTED_AVERAGE_DEPTH_COUNTS_FILE = new File(TEST_DATA_DIR, "exome-average-depth.output");
    private static final File TARGETS_FILE_WITHOUT_COORDINATES = new File(TEST_DATA_DIR, "exome-read-counts-test-targets-wo-coords.tsv");
    private static final File INEXISTENT_TARGETS_FILE = new File(TEST_DATA_DIR, "fantasy-exome-read-counts-test-targets.tsv");

    // Meta-parameters
    private static final List<Integer> TEST_MIN_MQ = Collections.unmodifiableList(Arrays.asList(0, 2, 10, 20, 30, 40, 60, 99, 255, 1000));
    private static final List<Integer> TEST_MIN_BQ = Collections.unmodifiableList(Arrays.asList(0, 2, 10, 20, 40, 60, 99));

    @Override
    public String getTestedClassName() {
        return CalculateTargetBaseCallCoverage.class.getSimpleName();
    }

    @Test(expectedExceptions = UserException.BadArgumentValue.class)
    public void testNegativeMinimumMQ() throws IOException {
        final File outputFile = createTempFile("ctc-test-", ".tsv");
        runCommandLine(
                new String[]{
                        "-" + TargetArgumentCollection.TARGET_FILE_SHORT_NAME,
                        TARGETS_FILE.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12878.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12778.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12872.toString(),
                        "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
                        outputFile.toString(),
                        "-" + CalculateTargetBaseCallCoverage.COVERAGE_UNIT_SHORT_NAME,
                        CoverageUnit.OVERLAPPING_READ.toString(),
                        "-" + CalculateTargetBaseCallCoverage.MINIMUM_MAPPING_QUALITY_SHORT_NAME,
                        "-10"
                });
    }

    @Test(expectedExceptions = UserException.BadArgumentValue.class)
    public void testNegativeMinimumBQ() throws IOException {
        final File outputFile = createTempFile("ctc-test-", ".tsv");
        runCommandLine(
                new String[]{
                        "-" + TargetArgumentCollection.TARGET_FILE_SHORT_NAME,
                        TARGETS_FILE.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12878.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12778.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12872.toString(),
                        "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
                        outputFile.toString(),
                        "-" + CalculateTargetBaseCallCoverage.COVERAGE_UNIT_SHORT_NAME,
                        CoverageUnit.OVERLAPPING_READ.toString(),
                        "-" + CalculateTargetBaseCallCoverage.MINIMUM_BASE_QUALITY_SHORT_NAME,
                        "-10"
                });
    }

    @Test(expectedExceptions = UserException.BadArgumentValue.class)
    public void testNegativeMaximumCoverage() throws IOException {
        final File outputFile = createTempFile("ctc-test-", ".tsv");
        runCommandLine(
                new String[]{
                        "-" + TargetArgumentCollection.TARGET_FILE_SHORT_NAME,
                        TARGETS_FILE.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12878.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12778.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12872.toString(),
                        "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
                        outputFile.toString(),
                        "-" + CalculateTargetBaseCallCoverage.COVERAGE_UNIT_SHORT_NAME,
                        CoverageUnit.OVERLAPPING_READ.toString(),
                        "-" + CalculateTargetBaseCallCoverage.MAXIMUM_COVERAGE_SHORT_NAME,
                        "-10"
                });
    }

    @Test(expectedExceptions = UserException.BadArgumentValue.class)
    public void testNaNMaximumCoverage() throws IOException {
        final File outputFile = createTempFile("ctc-test-", ".tsv");
        runCommandLine(
                new String[]{
                        "-" + TargetArgumentCollection.TARGET_FILE_SHORT_NAME,
                        TARGETS_FILE.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12878.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12778.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12872.toString(),
                        "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
                        outputFile.toString(),
                        "-" + CalculateTargetBaseCallCoverage.COVERAGE_UNIT_SHORT_NAME,
                        CoverageUnit.OVERLAPPING_READ.toString(),
                        "-" + CalculateTargetBaseCallCoverage.MAXIMUM_COVERAGE_SHORT_NAME,
                        "nan"
                });
    }

    @Test
    public void testSimpleRunReadsCoverage() throws IOException {
        final File outputFile = createTempFile("ctc-test-",".tsv");
        runCommandLine(
                new String[]{
                        "-" + TargetArgumentCollection.TARGET_FILE_SHORT_NAME,
                        TARGETS_FILE.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12878.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12778.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12872.toString(),
                        "-" + CalculateTargetBaseCallCoverage.COVERAGE_UNIT_SHORT_NAME,
                        CoverageUnit.OVERLAPPING_READ.toString(),
                        "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
                        outputFile.toString()
                });
        Assert.assertTrue(outputFile.canRead());
        final ReadCountCollection outputCounts = ReadCountCollectionUtils.parse(outputFile);
        final ReadCountCollection expectedCounts = ReadCountCollectionUtils.parse(EXPECTED_COUNTS_FILE);

        Assert.assertEquals(outputCounts.columnNames(), expectedCounts.columnNames());
        Assert.assertEquals(outputCounts.counts().getData(), expectedCounts.counts().getData());
    }

    @Test
    public void testSimpleRunWithMaximumCoverage() throws IOException {
        final File outputFile = createTempFile("ctc-test-",".tsv");
        runCommandLine(
                new String[]{
                        "-" + TargetArgumentCollection.TARGET_FILE_SHORT_NAME,
                        TARGETS_FILE.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12878.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12778.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12872.toString(),
                        "-" + CalculateTargetBaseCallCoverage.MAXIMUM_COVERAGE_SHORT_NAME,
                        "9",
                        "-" + CalculateTargetBaseCallCoverage.MINIMUM_BASE_QUALITY_SHORT_NAME,
                        "" + (CalculateTargetBaseCallCoverage.MINIMUM_BASE_QUALITY_DEFAULT + 30),
                        "-" + CalculateTargetBaseCallCoverage.COVERAGE_UNIT_SHORT_NAME,
                        CoverageUnit.OVERLAPPING_READ.toString(),
                        "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
                        outputFile.toString()
                });
        Assert.assertTrue(outputFile.canRead());
        final ReadCountCollection outputCounts = ReadCountCollectionUtils.parse(outputFile);
        final ReadCountCollection expectedCounts = ReadCountCollectionUtils.parse(EXPECTED_COUNTS_MAX_OF_9_FILE);

        Assert.assertEquals(outputCounts.columnNames(), expectedCounts.columnNames());
        Assert.assertEquals(outputCounts.counts().getData(), expectedCounts.counts().getData(), Arrays.toString(outputCounts.counts().getData()[0]));
    }

    @Test
    public void testSimpleRunWithMinimumMQ() throws IOException {
        final File outputFile = createTempFile("ctc-test-",".tsv");
        runCommandLine(
                new String[]{
                        "-" + TargetArgumentCollection.TARGET_FILE_SHORT_NAME,
                        TARGETS_FILE.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12878.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12778.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12872.toString(),
                        "-" + CalculateTargetBaseCallCoverage.MINIMUM_MAPPING_QUALITY_SHORT_NAME,
                        "30",
                        "-" + CalculateTargetBaseCallCoverage.COVERAGE_UNIT_SHORT_NAME,
                        CoverageUnit.OVERLAPPING_READ.toString(),
                        "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
                        outputFile.toString()
                });
        Assert.assertTrue(outputFile.canRead());
        final ReadCountCollection outputCounts = ReadCountCollectionUtils.parse(outputFile);
        final ReadCountCollection expectedCounts = ReadCountCollectionUtils.parse(EXPECTED_COUNTS_MIN_MQ_30_FILE);

        Assert.assertEquals(outputCounts.columnNames(), expectedCounts.columnNames());
        Assert.assertEquals(outputCounts.counts().getData(), expectedCounts.counts().getData(), Arrays.toString(outputCounts.counts().getData()[0]));
    }

    @Test
    public void testSingleSampleRunOverlappingReads() throws IOException {
        final List<SimpleInterval> intervals;
        try (final TargetTableReader targetReader = new TargetTableReader(TARGETS_FILE)) {
            intervals = targetReader.stream().map(Target::getInterval).collect(Collectors.toList());
        }
        final List<long[]> expectedCounts = calculateOverlappingReads(TEST_BAM_NA12878, "NA12878", intervals, TEST_MIN_MQ);
        for (int i = 0; i < TEST_MIN_MQ.size(); i++) {
            final int minMQ = TEST_MIN_MQ.get(i);
            final File outputFile = createTempFile("ctc-test-", ".tsv");
            runCommandLine(
                    new String[]{
                            "-" + TargetArgumentCollection.TARGET_FILE_SHORT_NAME,
                            TARGETS_FILE.toString(),
                            "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                            TEST_BAM_NA12878.toString(),
                            "-" + CalculateTargetBaseCallCoverage.COVERAGE_UNIT_SHORT_NAME,
                            CoverageUnit.OVERLAPPING_READ.toString(),
                            "-" + CalculateTargetBaseCallCoverage.MINIMUM_MAPPING_QUALITY_SHORT_NAME,
                            String.valueOf(minMQ),
                            "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
                            outputFile.toString()
                    });
            Assert.assertTrue(outputFile.canRead());
            final ReadCountCollection outputCounts = ReadCountCollectionUtils.parse(outputFile);
            final double[] actualSampleCounts = outputCounts.counts().getColumn(0);
            for (int j = 0; j < outputCounts.targets().size(); j++) {
                Assert.assertEquals(actualSampleCounts[j], (double) expectedCounts.get(i)[j], 0.00001, "" + j + " intervals " + intervals.get(j) + " and " + outputCounts.targets().get(j).getInterval());
            }
            outputFile.delete();
        }
    }

    @Test
    public void testSingleSampleRunAverageDepth() throws IOException {
        final List<SimpleInterval> intervals;
        try (final TargetTableReader targetReader = new TargetTableReader(TARGETS_FILE)) {
            intervals = targetReader.stream().map(Target::getInterval).collect(Collectors.toList());
        }
        for (int i = 0; i < TEST_MIN_MQ.size(); i++) {
            final int minMQ = TEST_MIN_MQ.get(i);
            final List<long[]> expectedCounts = calculateOveralappingBaseCalls(TEST_BAM_NA12878, "NA12878", intervals, minMQ, TEST_MIN_BQ);
            for (int j = 0; j < TEST_MIN_BQ.size(); j++) {
                final int minBQ = TEST_MIN_BQ.get(j);
                final File outputFile = createTempFile("ctc-test-", ".tsv");
                runCommandLine(
                        new String[]{
                                "-" + TargetArgumentCollection.TARGET_FILE_SHORT_NAME,
                                TARGETS_FILE.toString(),
                                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                                TEST_BAM_NA12878.toString(),
                                "-" + CalculateTargetBaseCallCoverage.COVERAGE_UNIT_SHORT_NAME,
                                CoverageUnit.AVERAGE_DEPTH.toString(),
                                "-" + CalculateTargetBaseCallCoverage.MINIMUM_MAPPING_QUALITY_SHORT_NAME,
                                String.valueOf(minMQ),
                                "-" + CalculateTargetBaseCallCoverage.MINIMUM_BASE_QUALITY_SHORT_NAME,
                                String.valueOf(minBQ),
                                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
                                outputFile.toString()
                        });
                Assert.assertTrue(outputFile.canRead());
                final ReadCountCollection outputCounts = ReadCountCollectionUtils.parse(outputFile);
                final double[] actualSampleCounts = outputCounts.counts().getColumn(0);
                for (int k = 0; k < outputCounts.targets().size(); k++) {
                    final double intervalWidth = intervals.get(k).getEnd() - intervals.get(k).getStart() + 1;
                    Assert.assertEquals(actualSampleCounts[k], (double) expectedCounts.get(j)[k] / intervalWidth, 0.00001, "" + k + " intervals " + intervals.get(k) + " and " + outputCounts.targets().get(k).getInterval());
                }
                outputFile.delete();
            }
        }
    }

    @Test
    public void testSimpleRunBaseCallCoverage() throws IOException {
        final File outputFile = createTempFile("ctc-test-", ".tsv");
        runCommandLine(
                new String[]{
                        "-" + TargetArgumentCollection.TARGET_FILE_SHORT_NAME,
                        TARGETS_FILE.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12878.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12778.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12872.toString(),
                        "-" + CalculateTargetBaseCallCoverage.COVERAGE_UNIT_SHORT_NAME,
                        CoverageUnit.AVERAGE_DEPTH.toString(),
                        "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
                        outputFile.toString()
                });
        Assert.assertTrue(outputFile.canRead());
        final ReadCountCollection outputCounts = ReadCountCollectionUtils.parse(outputFile);

        //BAMIO.openBAM()
        final ReadCountCollection expectedCounts = ReadCountCollectionUtils.parse(EXPECTED_AVERAGE_DEPTH_COUNTS_FILE);

        Assert.assertEquals(outputCounts.columnNames(), expectedCounts.columnNames());
        assertEqualDataMatrix(outputCounts.counts().getData(), expectedCounts.counts().getData());
    }

    private void assertEqualDataMatrix(final double[][] m1, final double[][] m2) {
        Assert.assertEquals(m1.length, m2.length);
        for (int i = 0; i < m1.length; i++) {
            Assert.assertEquals(m1[i].length, m2[i].length);
            for (int j = 0; j < m1[i].length; j++) {
                Assert.assertEquals(m1[i][j], m2[i][j], 0.00001, "i,j = " + i + "," + j);
            }
        }
    }

    private static CountingReadFilter makeBasicReadFilter(final SAMFileHeader header) {
        return new CountingReadFilter("Wellformed", new WellformedReadFilter(header))
                .and(new CountingReadFilter("Mapped", ReadFilterLibrary.MAPPED))
                .and(new CountingReadFilter("Not_Duplicate", ReadFilterLibrary.NOT_DUPLICATE))
                .and(new CountingReadFilter("Non_Zero_Reference_Length", ReadFilterLibrary.NON_ZERO_REFERENCE_LENGTH_ALIGNMENT));
    }

    /**
     * Slow but simple overlapping read count calculators.
     * @param bamFile
     * @param minMQ list of qualifying min MQ to consider.
     * @return never {@code }
     */
    private List<long[]> calculateOveralappingBaseCalls(final File bamFile, final String sample, final List<SimpleInterval> intervals, final int minMQ, final List<Integer> minBQ) throws IOException {
        final List<long[]> result = new ArrayList<>(minBQ.size());
        for (int i = 0; i < minBQ.size(); i++) {
            result.add(new long[intervals.size()]);
        }
        final SamReader reader = SamReaderFactory.makeDefault().open(bamFile);
        final CountingReadFilter baseReadFilter = makeBasicReadFilter(reader.getFileHeader());
        final SAMRecordIterator iterator = reader.iterator();
        iterator.forEachRemaining(read -> {
            if (!baseReadFilter.test(new SAMRecordToGATKReadAdapter(read))
                    || read.getMappingQuality() < minMQ || !read.getReadGroup().getSample().equals(sample)) {
                return;
            }
           final byte[] baseQualities = read.getBaseQualities();
            for (int i = 0; i < intervals.size(); i++) {
                final SimpleInterval simpleInterval = intervals.get(i);
                if (read.getReferenceName().equals(simpleInterval.getContig())) {
                    for (final AlignmentBlock block : read.getAlignmentBlocks()) {
                        final int blockStart = block.getReferenceStart();
                        final int blockEnd = blockStart + block.getLength() - 1;
                        if ((blockStart <= simpleInterval.getEnd() && blockEnd >= simpleInterval.getStart())) {
                            final int overlapStart = Math.max(blockStart, simpleInterval.getStart());
                            final int overlapEnd = Math.min(blockEnd, simpleInterval.getEnd());
                            final int readStart = block.getReadStart() + (overlapStart - blockStart);
                            final int readEnd = readStart + overlapEnd - overlapStart;
                            for (int j = readStart; j <= readEnd; j++) {

                                final int qual = baseQualities[j - 1];
                                for (int k = 0; k < minBQ.size(); k++) {
                                    if (minBQ.get(k) <= qual) {
                                        result.get(k)[i]++;
                                    }
                                }
                            }
                        }
                    }

                }
            }
        });
        reader.close();
        return result;
    }


    /**
     * Slow but simple overlapping read count calculators.
     * @param bamFile
     * @param minMQ list of qualifying min MQ to consider.
     * @return never {@code }
     */
    private List<long[]> calculateOverlappingReads(final File bamFile, final String sample, final List<SimpleInterval> intervals, final List<Integer> minMQ) throws IOException {
        final List<long[]> result = new ArrayList<>(minMQ.size());
        for (int i = 0; i < minMQ.size(); i++) {
            result.add(new long[intervals.size()]);
        }
        final SamReader reader = SamReaderFactory.makeDefault().open(bamFile);
        final CountingReadFilter baseReadFilter = makeBasicReadFilter(reader.getFileHeader());
        final SAMRecordIterator iterator = reader.iterator();
        iterator.forEachRemaining(read -> {
            if (!baseReadFilter.test(new SAMRecordToGATKReadAdapter(read))) {
                return;
            }
            if (!read.getReadGroup().getSample().equals(sample)) {
                return;
            }
            for (int i = 0; i < intervals.size(); i++) {
                final SimpleInterval simpleInterval = intervals.get(i);
                if (read.getReferenceName().equals(simpleInterval.getContig())) {
                    boolean foundOverlap = false;
                    for (final AlignmentBlock block : read.getAlignmentBlocks()) {
                        final int blockStart = block.getReferenceStart();
                        final int blockEnd = blockStart + block.getLength() - 1;
                        if ((foundOverlap = blockStart <= simpleInterval.getEnd() && blockEnd >= simpleInterval.getStart())) {
                            break;

                        }
                    }
                    if (foundOverlap) {
                        for (int j = 0; j < minMQ.size(); j++) {
                            if (minMQ.get(j) <= read.getMappingQuality()) {
                                result.get(j)[i]++;
                            }
                        }
                    }

                }
            }
        });
        reader.close();
        return result;
    }

    @Test(expectedExceptions = UserException.class)
    public void testWithoutTargetFile() throws IOException {
        final File outputFile = createTempFile("ctc-test-", ".tsv");
        runCommandLine(
                new String[]{
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12878.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12778.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12872.toString(),
                        "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
                        outputFile.toString()
                });
   }

    @Test(expectedExceptions = UserException.class)
    public void testWithMissingTargetFile() throws IOException {
        final File outputFile = createTempFile("ctc-test-", ".tsv");
        runCommandLine(
                new String[]{
                        "-" + TargetArgumentCollection.TARGET_FILE_SHORT_NAME,
                        INEXISTENT_TARGETS_FILE.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12878.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12778.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12872.toString(),
                        "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
                        outputFile.toString()
                });
    }

    @Test(expectedExceptions = UserException.class)
    public void testWithTargetFileWithoutCoordinates() throws IOException {
        final File outputFile = createTempFile("ctc-test-", ".tsv");
        runCommandLine(
                new String[]{
                        "-" + TargetArgumentCollection.TARGET_FILE_SHORT_NAME,
                        TARGETS_FILE_WITHOUT_COORDINATES.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12878.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12778.toString(),
                        "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME,
                        TEST_BAM_NA12872.toString(),
                        "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
                        outputFile.toString()
                });
    }
}
