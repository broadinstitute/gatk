package org.broadinstitute.hellbender.tools.copynumber.utils;

import com.google.common.collect.ImmutableSortedMap;
import com.google.common.collect.Sets;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.copynumber.formats.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion.SimpleAnnotatedGenomicRegion;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

public class CombineSegmentBreakpointsIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/copynumber/utils/";
    private static final String INPUT_SEGMENTS_FILE = TEST_SUB_DIR + "combine-segment-breakpoints-with-legacy-header.tsv";
    private static final String GROUND_TRUTH_SEGMENTS_FILE = TEST_SUB_DIR + "combine-segment-breakpoints-with-legacy-header-ground-truth.tsv";
    private static final String INPUT_SEGMENTS_FILE_WITH_DIFFERENT_HEADERS = TEST_SUB_DIR + "combine-segment-breakpoints-different-annotation-headers-with-legacy-header.tsv";
    public static final String REFERENCE_FILE = hg19_chr1_1M_Reference;

    @Test
    public void testRunWithExactSegments() throws IOException {
        // Segment intervals are the same in the input files.  Therefore, the union should only generate more columns.
        final File outputFile = File.createTempFile("combineseg_", ".tsv");
        final Set<String> columnSet = Sets.newHashSet("MEAN_LOG2_COPY_RATIO", "CALL", "Segment_Mean", "Segment_Call");
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + CopyNumberStandardArgument.SEGMENTS_FILE_SHORT_NAME);
        arguments.add(INPUT_SEGMENTS_FILE);
        arguments.add("-" + CopyNumberStandardArgument.SEGMENTS_FILE_SHORT_NAME);
        arguments.add(INPUT_SEGMENTS_FILE_WITH_DIFFERENT_HEADERS);
        arguments.add("-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME);
        arguments.add(REFERENCE_FILE);

        columnSet.forEach(s -> {
            arguments.add("-" + CombineSegmentBreakpoints.COLUMNS_OF_INTEREST_SHORT_NAME);
            arguments.add(s);
        });

        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        runCommandLine(arguments);

        Assert.assertTrue(outputFile.exists());

        final List<SimpleAnnotatedGenomicRegion> regions = SimpleAnnotatedGenomicRegion.readAnnotatedRegions(outputFile, Sets.newHashSet("MEAN_LOG2_COPY_RATIO", "CALL", "Segment_Mean", "Segment_Call"));
        Assert.assertEquals(regions.size(), 4);
        Assert.assertTrue(regions.stream().allMatch(r -> r.getAnnotations().size() == columnSet.size()));
        Assert.assertTrue(regions.stream().allMatch(r -> r.getAnnotations().keySet().containsAll(columnSet)));
        Assert.assertTrue(regions.stream().noneMatch(r -> r.getAnnotations().values().contains("")));
    }

    @Test
    public void testRunWithExactSameFiles() throws IOException {
        // Input files are exactly the same.  Therefore, the union should only generate more columns.
        final File outputFile = File.createTempFile("combineseg_", ".tsv");
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME);
        arguments.add(REFERENCE_FILE);
        arguments.add("-" + CopyNumberStandardArgument.SEGMENTS_FILE_SHORT_NAME);
        arguments.add(INPUT_SEGMENTS_FILE);
        arguments.add("-" + CopyNumberStandardArgument.SEGMENTS_FILE_SHORT_NAME);
        arguments.add(INPUT_SEGMENTS_FILE);
        arguments.add("-" + CombineSegmentBreakpoints.COLUMNS_OF_INTEREST_SHORT_NAME);
        arguments.add("MEAN_LOG2_COPY_RATIO");
        arguments.add("-" + CombineSegmentBreakpoints.COLUMNS_OF_INTEREST_SHORT_NAME);
        arguments.add("CALL");
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        runCommandLine(arguments);

        Assert.assertTrue(outputFile.exists());

        final List<SimpleAnnotatedGenomicRegion> regions = SimpleAnnotatedGenomicRegion.readAnnotatedRegions(outputFile, Sets.newHashSet("MEAN_LOG2_COPY_RATIO_1", "CALL_1", "MEAN_LOG2_COPY_RATIO_2", "CALL_2"));

        final Set<String> gtColumnSet = Sets.newHashSet("MEAN_LOG2_COPY_RATIO_1", "CALL_1", "MEAN_LOG2_COPY_RATIO_2", "CALL_2");
        Assert.assertEquals(regions.size(), 4);
        Assert.assertTrue(regions.stream().allMatch(r -> r.getAnnotations().size() == gtColumnSet.size()));
        Assert.assertTrue(regions.stream().allMatch(r -> r.getAnnotations().keySet().containsAll(gtColumnSet)));
        Assert.assertTrue(regions.stream().noneMatch(r -> r.getAnnotations().values().contains("")));
        Assert.assertTrue(regions.stream().allMatch(r -> r.getAnnotations().get("MEAN_LOG2_COPY_RATIO_1").equals(r.getAnnotations().get("MEAN_LOG2_COPY_RATIO_2"))));
        Assert.assertTrue(regions.stream().allMatch(r -> r.getAnnotations().get("CALL_1").equals(r.getAnnotations().get("CALL_2"))));
    }
    @Test
    public void testRunWithNotExactOverlaps() throws IOException {
        // This test is a bit more like the real world
        final File outputFile = File.createTempFile("combineseg_", ".tsv");
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME);
        arguments.add(REFERENCE_FILE);
        arguments.add("-" + CopyNumberStandardArgument.SEGMENTS_FILE_SHORT_NAME);
        arguments.add(INPUT_SEGMENTS_FILE_WITH_DIFFERENT_HEADERS);
        arguments.add("-" + CopyNumberStandardArgument.SEGMENTS_FILE_SHORT_NAME);
        arguments.add(GROUND_TRUTH_SEGMENTS_FILE);
        arguments.add("-" + CombineSegmentBreakpoints.COLUMNS_OF_INTEREST_SHORT_NAME);
        arguments.add("Segment_Mean");
        arguments.add("-" + CombineSegmentBreakpoints.COLUMNS_OF_INTEREST_SHORT_NAME);
        arguments.add("Segment_Call");
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        runCommandLine(arguments);

        Assert.assertTrue(outputFile.exists());
        final String SEGMENT_CALL_1 = "Segment_Call_1";
        final String SEGMENT_MEAN_1 = "Segment_Mean_1";
        final String SEGMENT_MEAN_2 = "Segment_Mean_2";
        final String SEGMENT_CALL_2 = "Segment_Call_2";
        final List<SimpleAnnotatedGenomicRegion> regions = SimpleAnnotatedGenomicRegion.readAnnotatedRegions(outputFile, Sets.newHashSet(SEGMENT_MEAN_1, SEGMENT_CALL_1, SEGMENT_MEAN_2, SEGMENT_CALL_2));
        Assert.assertEquals(regions.size(), 13);
        Assert.assertTrue(regions.stream().allMatch(r -> r.getAnnotations().size() == 4));
        assertUnionedSegFiles(SEGMENT_CALL_1, SEGMENT_MEAN_1, SEGMENT_MEAN_2, SEGMENT_CALL_2, regions);
    }

    private void assertUnionedSegFiles(final String segmentCall1, final String segmentMean1, final String segmentMean2,
                                       final String segmentCall2, final List<SimpleAnnotatedGenomicRegion> regions) {
        // Painstakingly made by hand
        Assert.assertEquals(regions.get(0), new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 5000, 10000),
                ImmutableSortedMap.of(segmentCall1, "", segmentCall2, "0", segmentMean1, "", segmentMean2, "-0.04")));
        Assert.assertEquals(regions.get(1), new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 10001, 10500),
                ImmutableSortedMap.of(segmentCall1, "0", segmentCall2, "0", segmentMean1, "-0.03", segmentMean2, "-0.04")));
        Assert.assertEquals(regions.get(2), new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 10501, 52499),
                ImmutableSortedMap.of(segmentCall1, "", segmentCall2, "0", segmentMean1, "", segmentMean2, "-0.04")));
        Assert.assertEquals(regions.get(3), new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 52500, 60000),
                ImmutableSortedMap.of(segmentCall1, "-", segmentCall2, "0", segmentMean1, "-0.76", segmentMean2, "-0.04")));
        Assert.assertEquals(regions.get(4), new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 60001, 69999),
                ImmutableSortedMap.of(segmentCall1, "-", segmentCall2, "", segmentMean1, "-0.76", segmentMean2, "")));
        Assert.assertEquals(regions.get(5), new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 70000, 100000),
                ImmutableSortedMap.of(segmentCall1, "-", segmentCall2, "-", segmentMean1, "-0.76", segmentMean2, "-0.8")));
        Assert.assertEquals(regions.get(6), new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 100001, 109750),
                ImmutableSortedMap.of(segmentCall1, "-", segmentCall2, "", segmentMean1, "-0.76", segmentMean2, "")));
        Assert.assertEquals(regions.get(7), new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 109751, 119999),
                ImmutableSortedMap.of(segmentCall1, "0", segmentCall2, "", segmentMean1, "-0.10", segmentMean2, "")));
        Assert.assertEquals(regions.get(8), new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 120000, 220000),
                ImmutableSortedMap.of(segmentCall1, "0", segmentCall2, "0", segmentMean1, "-0.10", segmentMean2, "-0.1")));
        Assert.assertEquals(regions.get(9), new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 220001, 229999),
                ImmutableSortedMap.of(segmentCall1, "0", segmentCall2, "", segmentMean1, "-0.10", segmentMean2, "")));
        Assert.assertEquals(regions.get(10), new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 230000, 230500),
                ImmutableSortedMap.of(segmentCall1, "0", segmentCall2, "-", segmentMean1, "-0.10", segmentMean2, "-0.8")));
        Assert.assertEquals(regions.get(11), new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 230501, 258500),
                ImmutableSortedMap.of(segmentCall1, "-", segmentCall2, "-", segmentMean1, "-0.60", segmentMean2, "-0.8")));
        Assert.assertEquals(regions.get(12), new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 258501, 300000),
                ImmutableSortedMap.of(segmentCall1, "", segmentCall2, "-", segmentMean1, "", segmentMean2, "-0.8")));
    }

    @Test
    public void testRunWithNotExactOverlapsAndLabels() throws IOException {
        final String TEST = "test";
        final String GT = "gt";
        // This test is a bit more like the real world
        final File outputFile = File.createTempFile("combineseg_", ".tsv");
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME);
        arguments.add(REFERENCE_FILE);
        arguments.add("-" + CopyNumberStandardArgument.SEGMENTS_FILE_SHORT_NAME);
        arguments.add(INPUT_SEGMENTS_FILE_WITH_DIFFERENT_HEADERS);
        arguments.add("-" + CopyNumberStandardArgument.SEGMENTS_FILE_SHORT_NAME);
        arguments.add(GROUND_TRUTH_SEGMENTS_FILE);
        arguments.add("-" + CombineSegmentBreakpoints.LABELS_SHORT_NAME);
        arguments.add(TEST);
        arguments.add("-" + CombineSegmentBreakpoints.LABELS_SHORT_NAME);
        arguments.add(GT);
        arguments.add("-" + CombineSegmentBreakpoints.COLUMNS_OF_INTEREST_SHORT_NAME);
        arguments.add("Segment_Mean");
        arguments.add("-" + CombineSegmentBreakpoints.COLUMNS_OF_INTEREST_SHORT_NAME);
        arguments.add("Segment_Call");
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        runCommandLine(arguments);

        Assert.assertTrue(outputFile.exists());

        final String SEGMENT_CALL_1 = "Segment_Call_" + TEST;
        final String SEGMENT_MEAN_1 = "Segment_Mean_" + TEST;
        final String SEGMENT_MEAN_2 = "Segment_Mean_" + GT;
        final String SEGMENT_CALL_2 = "Segment_Call_" + GT;
        final List<SimpleAnnotatedGenomicRegion> regions = SimpleAnnotatedGenomicRegion.readAnnotatedRegions(outputFile, Sets.newHashSet(SEGMENT_MEAN_1, SEGMENT_CALL_1, SEGMENT_MEAN_2, SEGMENT_CALL_2));
        Assert.assertEquals(regions.size(), 13);
        Assert.assertTrue(regions.stream().allMatch(r -> r.getAnnotations().size() == 4));
        assertUnionedSegFiles(SEGMENT_CALL_1, SEGMENT_MEAN_1, SEGMENT_MEAN_2, SEGMENT_CALL_2, regions);
    }
}
