package org.broadinstitute.hellbender.tools.copynumber.utils;

import com.google.common.collect.ImmutableSortedMap;
import com.google.common.collect.Sets;
import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedInterval;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedIntervalCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.stream.Stream;

public final class CombineSegmentBreakpointsIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_SUB_DIR = toolsTestDir + "/copynumber/utils/";
    private static final String INPUT_SEGMENTS_FILE = TEST_SUB_DIR + "combine-segment-breakpoints.seg";
    private static final String INPUT_SEGMENTS_FILE_NO_SAMHEADER = TEST_SUB_DIR + "combine-segment-breakpoints-no-samheader.seg";
    private static final String INPUT_SEGMENTS_FILE_ALT_SAMHEADER = TEST_SUB_DIR + "combine-segment-breakpoints-alt-samheader.seg";
    private static final String GROUND_TRUTH_SEGMENTS_FILE = TEST_SUB_DIR + "combine-segment-breakpoints-comparison.seg";
    private static final String GROUND_TRUTH_SEGMENTS_FILE_NO_SAMHEADER = TEST_SUB_DIR + "combine-segment-breakpoints-comparison-no-samheader.seg";
    private static final String SEGMENTS_FILE_RECAPSEG = TEST_SUB_DIR + "combine-segment-breakpoints-comparison-recapseg.seg";
    private static final String SEGMENTS_FILE_JABBA = TEST_SUB_DIR + "combine-segment-breakpoints-comparison-jabba.seg";
    private static final String SEGMENTS_FILE_JABBA_EXPECTED = TEST_SUB_DIR + "combine-segment-breakpoints-comparison-jabba-expected.seg";
    private static final String SEGMENTS_FILE_PCAWG_CONSENSUS = TEST_SUB_DIR + "combine-segment-breakpoints-comparison-pcawg-consensus.seg";
    private static final String SEGMENTS_FILE_PCAWG_CONSENSUS_EXPECTED = TEST_SUB_DIR + "combine-segment-breakpoints-comparison-pcawg-consensus-expected.seg";
    private static final String SEGMENTS_FILE_UCSC_TRACK = TEST_SUB_DIR + "combine-segment-breakpoints-comparison-ucsc-track.seg";
    private static final String SEGMENTS_FILE_UCSC_TRACK_EXPECTED = TEST_SUB_DIR + "combine-segment-breakpoints-comparison-ucsc-track-expected.seg";
    private static final String SEGMENTS_FILE_ONCOTATOR_GENE_LIST = TEST_SUB_DIR + "combine-segment-breakpoints-comparison-oncotator-gene-list.seg";
    private static final String INPUT_SEGMENTS_FILE_WITH_DIFFERENT_HEADERS = TEST_SUB_DIR + "combine-segment-breakpoints-different-annotation-headers.seg";
    private static final String REFERENCE_FILE = hg19_chr1_1M_Reference;
    private static final String SEGMENT_CALL_1 = "CALL";
    private static final String SEGMENT_MEAN_1 = "MEAN_LOG2_COPY_RATIO";
    private static final String SEGMENT_MEAN_2 = "Segment_Mean";
    private static final String SEGMENT_CALL_2 = "Segment_Call";

    @Test
    public void testRunWithExactSegments() throws IOException {
        // Segment intervals are the same in the input files.  Therefore, the union should only generate more columns.
        final File outputFile = File.createTempFile("combineseg_", ".seg");
        final Set<String> columnSet = Sets.newHashSet("MEAN_LOG2_COPY_RATIO", "CALL", "Segment_Mean", "Segment_Call");
        final List<String> arguments = new ArrayList<>();
        arguments.add("--" + CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME);
        arguments.add(INPUT_SEGMENTS_FILE);
        arguments.add("--" + CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME);
        arguments.add(INPUT_SEGMENTS_FILE_WITH_DIFFERENT_HEADERS);
        arguments.add("-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME);
        arguments.add(REFERENCE_FILE);

        columnSet.forEach(s -> {
            arguments.add("--" + CombineSegmentBreakpoints.COLUMNS_OF_INTEREST_LONG_NAME);
            arguments.add(s);
        });

        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        runCommandLine(arguments);

        Assert.assertTrue(outputFile.exists());

        final AnnotatedIntervalCollection regions = AnnotatedIntervalCollection.create(outputFile.toPath(), Sets.newHashSet("MEAN_LOG2_COPY_RATIO", "CALL", "Segment_Mean", "Segment_Call"));
        Assert.assertEquals(regions.size(), 4);
        Assert.assertTrue(regions.getRecords().stream().allMatch(r -> r.getAnnotations().size() == columnSet.size()));
        Assert.assertTrue(regions.getRecords().stream().allMatch(r -> r.getAnnotations().keySet().containsAll(columnSet)));
        Assert.assertTrue(regions.getRecords().stream().noneMatch(r -> r.getAnnotations().values().contains("")));
    }

    @Test
    public void testRunWithExactSameFiles() throws IOException {
        // Input files are exactly the same.  Therefore, the union should only generate more columns.
        final File outputFile = File.createTempFile("combineseg_", ".seg");
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME);
        arguments.add(REFERENCE_FILE);
        arguments.add("--" + CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME);
        arguments.add(INPUT_SEGMENTS_FILE);
        arguments.add("--" + CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME);
        arguments.add(INPUT_SEGMENTS_FILE);
        arguments.add("--" + CombineSegmentBreakpoints.COLUMNS_OF_INTEREST_LONG_NAME);
        arguments.add("MEAN_LOG2_COPY_RATIO");
        arguments.add("--" + CombineSegmentBreakpoints.COLUMNS_OF_INTEREST_LONG_NAME);
        arguments.add("CALL");
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        runCommandLine(arguments);

        Assert.assertTrue(outputFile.exists());

        final AnnotatedIntervalCollection regions = AnnotatedIntervalCollection.create(outputFile.toPath(), Sets.newHashSet("MEAN_LOG2_COPY_RATIO_1", "CALL_1", "MEAN_LOG2_COPY_RATIO_2", "CALL_2"));
        final Set<String> gtColumnSet = Sets.newHashSet("MEAN_LOG2_COPY_RATIO_1", "CALL_1", "MEAN_LOG2_COPY_RATIO_2", "CALL_2");
        Assert.assertEquals(regions.size(), 4);
        Assert.assertTrue(regions.getRecords().stream().allMatch(r -> r.getAnnotations().size() == gtColumnSet.size()));
        Assert.assertTrue(regions.getRecords().stream().allMatch(r -> r.getAnnotations().keySet().containsAll(gtColumnSet)));
        Assert.assertTrue(regions.getRecords().stream().noneMatch(r -> r.getAnnotations().values().contains("")));
        Assert.assertTrue(regions.getRecords().stream().allMatch(r -> r.getAnnotations().get("MEAN_LOG2_COPY_RATIO_1").equals(r.getAnnotations().get("MEAN_LOG2_COPY_RATIO_2"))));
        Assert.assertTrue(regions.getRecords().stream().allMatch(r -> r.getAnnotations().get("CALL_1").equals(r.getAnnotations().get("CALL_2"))));
    }
    @Test
    public void testRunWithoutExactOverlaps() throws IOException {
        // This test is a bit more like the real world
        final File outputFile = File.createTempFile("combineseg_", ".seg");
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME);
        arguments.add(REFERENCE_FILE);
        arguments.add("--" + CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME);
        arguments.add(INPUT_SEGMENTS_FILE_WITH_DIFFERENT_HEADERS);
        arguments.add("--" + CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME);
        arguments.add(GROUND_TRUTH_SEGMENTS_FILE);
        arguments.add("--" + CombineSegmentBreakpoints.COLUMNS_OF_INTEREST_LONG_NAME);
        arguments.add("Segment_Mean");
        arguments.add("--" + CombineSegmentBreakpoints.COLUMNS_OF_INTEREST_LONG_NAME);
        arguments.add("Segment_Call");
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        runCommandLine(arguments);

        Assert.assertTrue(outputFile.exists());
        final String SEGMENT_CALL_1 = "Segment_Call_1";
        final String SEGMENT_MEAN_1 = "Segment_Mean_1";
        final String SEGMENT_MEAN_2 = "Segment_Mean_2";
        final String SEGMENT_CALL_2 = "Segment_Call_2";
        final AnnotatedIntervalCollection regions = AnnotatedIntervalCollection.create(outputFile.toPath(), Sets.newHashSet(SEGMENT_MEAN_1, SEGMENT_CALL_1, SEGMENT_MEAN_2, SEGMENT_CALL_2));
        Assert.assertEquals(regions.size(), 13);
        Assert.assertTrue(regions.getRecords().stream().allMatch(r -> r.getAnnotations().size() == 4));
        assertUnionedSegFiles(SEGMENT_CALL_1, SEGMENT_MEAN_1, SEGMENT_MEAN_2, SEGMENT_CALL_2, regions.getRecords());
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testFailureIfNotAllColsOfInterestExist() throws IOException {
        // This test is a bit more like the real world
        final File outputFile = File.createTempFile("combineseg_", ".seg");
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME);
        arguments.add(REFERENCE_FILE);
        arguments.add("--" + CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME);
        arguments.add(INPUT_SEGMENTS_FILE_WITH_DIFFERENT_HEADERS);
        arguments.add("--" + CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME);
        arguments.add(GROUND_TRUTH_SEGMENTS_FILE);
        arguments.add("--" + CombineSegmentBreakpoints.COLUMNS_OF_INTEREST_LONG_NAME);
        arguments.add("Segment_Mean_THAT_DOES_NOT_EXIST");
        arguments.add("--" + CombineSegmentBreakpoints.COLUMNS_OF_INTEREST_LONG_NAME);
        arguments.add("Segment_Call");
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        runCommandLine(arguments);

    }

    private void assertUnionedSegFiles(final String segmentCall1, final String segmentMean1, final String segmentMean2,
                                       final String segmentCall2, final List<AnnotatedInterval> regions) {
        // Painstakingly made by hand
        Assert.assertEquals(regions.get(0), new AnnotatedInterval(new SimpleInterval("1", 5000, 10000),
                ImmutableSortedMap.of(segmentCall1, "", segmentCall2, "0", segmentMean1, "", segmentMean2, "-0.04")));
        Assert.assertEquals(regions.get(1), new AnnotatedInterval(new SimpleInterval("1", 10001, 10500),
                ImmutableSortedMap.of(segmentCall1, "0", segmentCall2, "0", segmentMean1, "-0.03", segmentMean2, "-0.04")));
        Assert.assertEquals(regions.get(2), new AnnotatedInterval(new SimpleInterval("1", 10501, 52499),
                ImmutableSortedMap.of(segmentCall1, "", segmentCall2, "0", segmentMean1, "", segmentMean2, "-0.04")));
        Assert.assertEquals(regions.get(3), new AnnotatedInterval(new SimpleInterval("1", 52500, 60000),
                ImmutableSortedMap.of(segmentCall1, "-", segmentCall2, "0", segmentMean1, "-0.76", segmentMean2, "-0.04")));
        Assert.assertEquals(regions.get(4), new AnnotatedInterval(new SimpleInterval("1", 60001, 69999),
                ImmutableSortedMap.of(segmentCall1, "-", segmentCall2, "", segmentMean1, "-0.76", segmentMean2, "")));
        Assert.assertEquals(regions.get(5), new AnnotatedInterval(new SimpleInterval("1", 70000, 100000),
                ImmutableSortedMap.of(segmentCall1, "-", segmentCall2, "-", segmentMean1, "-0.76", segmentMean2, "-0.8")));
        Assert.assertEquals(regions.get(6), new AnnotatedInterval(new SimpleInterval("1", 100001, 109750),
                ImmutableSortedMap.of(segmentCall1, "-", segmentCall2, "", segmentMean1, "-0.76", segmentMean2, "")));
        Assert.assertEquals(regions.get(7), new AnnotatedInterval(new SimpleInterval("1", 109751, 119999),
                ImmutableSortedMap.of(segmentCall1, "0", segmentCall2, "", segmentMean1, "-0.10", segmentMean2, "")));
        Assert.assertEquals(regions.get(8), new AnnotatedInterval(new SimpleInterval("1", 120000, 220000),
                ImmutableSortedMap.of(segmentCall1, "0", segmentCall2, "0", segmentMean1, "-0.10", segmentMean2, "-0.1")));
        Assert.assertEquals(regions.get(9), new AnnotatedInterval(new SimpleInterval("1", 220001, 229999),
                ImmutableSortedMap.of(segmentCall1, "0", segmentCall2, "", segmentMean1, "-0.10", segmentMean2, "")));
        Assert.assertEquals(regions.get(10), new AnnotatedInterval(new SimpleInterval("1", 230000, 230500),
                ImmutableSortedMap.of(segmentCall1, "0", segmentCall2, "-", segmentMean1, "-0.10", segmentMean2, "-0.8")));
        Assert.assertEquals(regions.get(11), new AnnotatedInterval(new SimpleInterval("1", 230501, 258500),
                ImmutableSortedMap.of(segmentCall1, "-", segmentCall2, "-", segmentMean1, "-0.60", segmentMean2, "-0.8")));
        Assert.assertEquals(regions.get(12), new AnnotatedInterval(new SimpleInterval("1", 258501, 300000),
                ImmutableSortedMap.of(segmentCall1, "", segmentCall2, "-", segmentMean1, "", segmentMean2, "-0.8")));
    }

    @Test
    public void testRunWithoutExactOverlapsAndLabels() throws IOException {
        final String TEST = "test";
        final String GT = "gt";
        // This test is a bit more like the real world
        final File outputFile = File.createTempFile("combineseg_", ".seg");
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME);
        arguments.add(REFERENCE_FILE);
        arguments.add("--" + CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME);
        arguments.add(INPUT_SEGMENTS_FILE_WITH_DIFFERENT_HEADERS);
        arguments.add("--" + CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME);
        arguments.add(GROUND_TRUTH_SEGMENTS_FILE);
        arguments.add("--" + CombineSegmentBreakpoints.LABELS_LONG_NAME);
        arguments.add(TEST);
        arguments.add("--" + CombineSegmentBreakpoints.LABELS_LONG_NAME);
        arguments.add(GT);
        arguments.add("--" + CombineSegmentBreakpoints.COLUMNS_OF_INTEREST_LONG_NAME);
        arguments.add("Segment_Mean");
        arguments.add("--" + CombineSegmentBreakpoints.COLUMNS_OF_INTEREST_LONG_NAME);
        arguments.add("Segment_Call");
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        runCommandLine(arguments);

        Assert.assertTrue(outputFile.exists());

        final String SEGMENT_CALL_1 = "Segment_Call_" + TEST;
        final String SEGMENT_MEAN_1 = "Segment_Mean_" + TEST;
        final String SEGMENT_MEAN_2 = "Segment_Mean_" + GT;
        final String SEGMENT_CALL_2 = "Segment_Call_" + GT;
        final AnnotatedIntervalCollection regions = AnnotatedIntervalCollection.create(outputFile.toPath(), Sets.newHashSet(SEGMENT_MEAN_1, SEGMENT_CALL_1, SEGMENT_MEAN_2, SEGMENT_CALL_2));
        Assert.assertEquals(regions.size(), 13);
        Assert.assertTrue(regions.getRecords().stream().allMatch(r -> r.getAnnotations().size() == 4));
        assertUnionedSegFiles(SEGMENT_CALL_1, SEGMENT_MEAN_1, SEGMENT_MEAN_2, SEGMENT_CALL_2, regions.getRecords());
    }

    @Test(expectedExceptions = UserException.class)
    public void testErrorWithoutSamHeadersAndNoReference() throws IOException {

        final File outputFile = File.createTempFile("combineseg_", ".seg");
        runCombineSegmentBreakpoints(GROUND_TRUTH_SEGMENTS_FILE_NO_SAMHEADER, INPUT_SEGMENTS_FILE_NO_SAMHEADER,
                outputFile, null);
    }

    @Test
    public void testWithoutSamHeaderAndWithReference() throws IOException {

        // This test is a bit more like the real world
        final File outputFile = File.createTempFile("combineseg_", ".seg");
        runCombineSegmentBreakpoints(GROUND_TRUTH_SEGMENTS_FILE_NO_SAMHEADER, INPUT_SEGMENTS_FILE_NO_SAMHEADER, outputFile, REFERENCE_FILE);

        Assert.assertTrue(outputFile.exists());

        final AnnotatedIntervalCollection regions = AnnotatedIntervalCollection.create(outputFile.toPath(), Sets.newHashSet(SEGMENT_MEAN_1, SEGMENT_CALL_1, SEGMENT_MEAN_2, SEGMENT_CALL_2));
        Assert.assertEquals(regions.size(), 13);
        Assert.assertTrue(regions.getRecords().stream().allMatch(r -> r.getAnnotations().size() == 4));
        assertUnionedSegFiles(SEGMENT_CALL_1, SEGMENT_MEAN_1, SEGMENT_MEAN_2, SEGMENT_CALL_2, regions.getRecords());
        Assert.assertEquals(regions.getComments().size(), 3+3);
        Assert.assertEquals(regions.getComments().get(0), " This is a comment");
        Assert.assertEquals(regions.getComments().get(1), " This is another comment");
        Assert.assertEquals(regions.getComments().get(2), " This is yet another comment");
    }

    @Test
    public void testRunAndOneWithoutSamHeaderAndNoReference() throws IOException {

        // This test is a bit more like the real world
        final File outputFile = File.createTempFile("combineseg_", ".seg");
        runCombineSegmentBreakpoints(INPUT_SEGMENTS_FILE, GROUND_TRUTH_SEGMENTS_FILE_NO_SAMHEADER, outputFile, null);

        Assert.assertTrue(outputFile.exists());

        final AnnotatedIntervalCollection regions = AnnotatedIntervalCollection.create(outputFile.toPath(), Sets.newHashSet(SEGMENT_MEAN_1, SEGMENT_CALL_1, SEGMENT_MEAN_2, SEGMENT_CALL_2));
        Assert.assertEquals(regions.size(), 13);
        // Reminder:  Three comments are added to all outputs.
        Assert.assertTrue(regions.getRecords().stream().allMatch(r -> r.getAnnotations().size() == (1+3)));
        assertUnionedSegFiles(SEGMENT_CALL_1, SEGMENT_MEAN_1, SEGMENT_MEAN_2, SEGMENT_CALL_2, regions.getRecords());
        Assert.assertEquals(regions.getComments().size(), (1+3));
        Assert.assertEquals(regions.getComments().get(0), " This is a comment");

        // Can't simply test if the sam header is equal to the input, since the sort and group order gets overridden.
        final SAMFileHeader samFileHeader1 = AnnotatedIntervalCollection.create(new File(INPUT_SEGMENTS_FILE).toPath(), null).getSamFileHeader();
        Assert.assertEquals(samFileHeader1.getReadGroups(), regions.getSamFileHeader().getReadGroups());
        Assert.assertEquals(regions.getSamFileHeader().getSortOrder(), SAMFileHeader.SortOrder.coordinate);
        Assert.assertEquals(samFileHeader1.getSequenceDictionary(), regions.getSamFileHeader().getSequenceDictionary());
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testDisagreeingSamHeadersAndNoReference() throws IOException {
        final File outputFile = File.createTempFile("combineseg_", ".seg");
        runCombineSegmentBreakpoints(INPUT_SEGMENTS_FILE_ALT_SAMHEADER, GROUND_TRUTH_SEGMENTS_FILE, outputFile, null);
    }

    /** Reference is ignored if SAM Headers are specified. */
    @Test(expectedExceptions = UserException.BadInput.class)
    public void testDisagreeingSamHeadersAndReference() throws IOException {
        final File outputFile = File.createTempFile("combineseg_", ".seg");
        runCombineSegmentBreakpoints(INPUT_SEGMENTS_FILE_ALT_SAMHEADER, GROUND_TRUTH_SEGMENTS_FILE, outputFile, REFERENCE_FILE);
    }

    private void runCombineSegmentBreakpoints(final String file1, final String file2, final File outputFile, final String refFile) {

        runCombineSegmentBreakpoints(file1, file2, outputFile, refFile, SEGMENT_CALL_1, SEGMENT_MEAN_1, SEGMENT_CALL_2,
                SEGMENT_MEAN_2);
    }

    private void runCombineSegmentBreakpoints(final String file1, final String file2, final File outputFile, final String refFile,
                                              final String ... columnsOfInterest) {

        final List<String> arguments = new ArrayList<>();
        arguments.add("--" + CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME);
        arguments.add(file1);
        arguments.add("--" + CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME);
        arguments.add(file2);
        Stream.of(columnsOfInterest).forEach(ci -> {
            arguments.add("--" + CombineSegmentBreakpoints.COLUMNS_OF_INTEREST_LONG_NAME);
            arguments.add(ci);
        });
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        if (refFile != null) {
            arguments.add("-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME);
            arguments.add(refFile);
        }

        runCommandLine(arguments);
    }

    @DataProvider(name="oldFormats")
    public Object[][] provideOldFormats() {
        return new Object[][]{
                {INPUT_SEGMENTS_FILE, SEGMENTS_FILE_RECAPSEG, "Segment_Mean", "Segment_Call", 4},
                {INPUT_SEGMENTS_FILE, SEGMENTS_FILE_ONCOTATOR_GENE_LIST, "segment_mean", "segment_call", 4}
        };
    }

    @Test(dataProvider = "oldFormats")
    public void testOldFormats(final String inputSegFile, final String oldFormatSegFile, final String meanColumn, final String callColumn, int numAnnotations) throws IOException {
        final File outputFile = File.createTempFile("combineseg_old_formats_", ".seg");
        runCombineSegmentBreakpoints(inputSegFile, oldFormatSegFile, outputFile, REFERENCE_FILE, "MEAN_LOG2_COPY_RATIO", "CALL", meanColumn, callColumn);
        final AnnotatedIntervalCollection regions = AnnotatedIntervalCollection.create(outputFile.toPath(),
                Sets.newHashSet("MEAN_LOG2_COPY_RATIO", "CALL", meanColumn, callColumn));
        Assert.assertEquals(regions.size(), 13);
        Assert.assertTrue(regions.getRecords().stream().allMatch(r -> r.getAnnotations().size() == numAnnotations));
        assertUnionedSegFiles("CALL", "MEAN_LOG2_COPY_RATIO", meanColumn, callColumn, regions.getRecords());
    }

    @Test(dataProvider = "realTsvFormats")
    public void testRealFormats(final String segFile, final String externalSegFile, final String expectedResult, final String ... colsOfInterest) throws IOException {

        final File outputFile = File.createTempFile("combineseg_realformats_gt_", ".seg");
        runCombineSegmentBreakpoints(segFile, externalSegFile, outputFile, REFERENCE_FILE,
                colsOfInterest);

        final AnnotatedIntervalCollection regions = AnnotatedIntervalCollection.create(outputFile.toPath(),
                Sets.newHashSet(colsOfInterest));
        final AnnotatedIntervalCollection expectedRegions = AnnotatedIntervalCollection
                .create(new File(expectedResult).toPath(), null);
        Assert.assertEquals(regions, expectedRegions);

    }

    @DataProvider(name = "realTsvFormats")
    public Object[][] provideRealTsvFormats() {
        return new Object[][]{
                {INPUT_SEGMENTS_FILE, SEGMENTS_FILE_UCSC_TRACK, SEGMENTS_FILE_UCSC_TRACK_EXPECTED, "MEAN_LOG2_COPY_RATIO", "CALL", "type"},
                {INPUT_SEGMENTS_FILE, SEGMENTS_FILE_JABBA, SEGMENTS_FILE_JABBA_EXPECTED, "MEAN_LOG2_COPY_RATIO", "CALL", "cn"},
                {INPUT_SEGMENTS_FILE, SEGMENTS_FILE_PCAWG_CONSENSUS, SEGMENTS_FILE_PCAWG_CONSENSUS_EXPECTED, "MEAN_LOG2_COPY_RATIO", "CALL", "consensus_major_cn", "consensus_minor_cn", "final_total_cn", "star"}
        };
    }
}
