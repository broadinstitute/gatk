package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Unit tests for {@link SegmentedGenome}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class SegmentedGenomeUnitTest extends BaseTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/exome/";

    private static final File TARGET_FILE_SMALL_SEGMENT_MERGING
            = new File(TEST_SUB_DIR, "targets-for-small-segment-merging-base.tsv");
    private static final File SNP_FILE_SMALL_SEGMENT_MERGING
            = new File(TEST_SUB_DIR, "snps-for-small-segment-merging-base.tsv");
    private static final File SEGMENT_FILE_SMALL_SEGMENT_MERGING
            = new File(TEST_SUB_DIR, "segments-for-small-segment-merging-base.seg");
    private static final File SEGMENT_FILE_SMALL_SEGMENT_MERGING_NO_SMALL
            = new File(TEST_SUB_DIR, "segments-for-small-segment-merging-no-small.seg");
    private static final String SAMPLE_NAME = "test";

    //a segment is small if number of targets it contains is strictly less than SMALL_SEGMENT_TARGET_NUMBER_THRESHOLD
    private static final int SMALL_SEGMENT_TARGET_NUMBER_THRESHOLD = 3;

    //the base test case starts with identical sets of 10 segments on chromosomes 1 & 3 (with 7 small segments each)
    //and a single small segment (containing only a single target) on chromosome 2
    @Test
    public void testCountSmallSegments() {
        final Genome genome = new Genome(TARGET_FILE_SMALL_SEGMENT_MERGING, SNP_FILE_SMALL_SEGMENT_MERGING,
                SAMPLE_NAME);
        final SegmentedGenome model = new SegmentedGenome(SEGMENT_FILE_SMALL_SEGMENT_MERGING, genome);

        final int resultNumberOfSegmentsBeforeMerging = model.getSegments().size();
        final int resultNumberOfSmallSegmentsBeforeMerging =
                countSmallSegments(model, SMALL_SEGMENT_TARGET_NUMBER_THRESHOLD);
        final int expectedNumberOfSegmentsBeforeMerging = 21;
        final int expectedNumberOfSmallSegmentsBeforeMerging = 15;
        Assert.assertEquals(resultNumberOfSegmentsBeforeMerging, expectedNumberOfSegmentsBeforeMerging);
        Assert.assertEquals(resultNumberOfSmallSegmentsBeforeMerging, expectedNumberOfSmallSegmentsBeforeMerging);
    }

    @Test
    public void testSmallSegmentMergingBase() {
        final Genome genome = new Genome(TARGET_FILE_SMALL_SEGMENT_MERGING, SNP_FILE_SMALL_SEGMENT_MERGING,
                SAMPLE_NAME);
        final SegmentedGenome model = new SegmentedGenome(SEGMENT_FILE_SMALL_SEGMENT_MERGING, genome);

        final SegmentedGenome modelMerged = model.mergeSmallSegments(SMALL_SEGMENT_TARGET_NUMBER_THRESHOLD);

        final List<SimpleInterval> result = modelMerged.getSegments();

        final int resultNumberOfSegmentsAfterMerging = modelMerged.getSegments().size();
        final int resultNumberOfSmallSegmentsAfterMerging =
                countSmallSegments(modelMerged, SMALL_SEGMENT_TARGET_NUMBER_THRESHOLD);
        final int expectedNumberOfSegmentsAfterMerging = 8;
        final int expectedNumberOfSmallSegmentsAfterMerging = 0;
        Assert.assertEquals(resultNumberOfSegmentsAfterMerging, expectedNumberOfSegmentsAfterMerging);
        Assert.assertEquals(resultNumberOfSmallSegmentsAfterMerging, expectedNumberOfSmallSegmentsAfterMerging);

        final List<SimpleInterval> expected = new ArrayList<>(Arrays.asList(
                new SimpleInterval("1", 141, 278),
                new SimpleInterval("1", 311, 439),
                new SimpleInterval("1", 466, 589),
                new SimpleInterval("1", 612, 808),
                new SimpleInterval("3", 141, 278),
                new SimpleInterval("3", 311, 439),
                new SimpleInterval("3", 466, 589),
                new SimpleInterval("3", 612, 808)
        ));

        Assert.assertEquals(result, expected);
    }

    //test case with no small segments; model should not change after merging
    @Test
    public void testSmallSegmentMergingNoSmallSegments() {
        final Genome genome = new Genome(TARGET_FILE_SMALL_SEGMENT_MERGING, SNP_FILE_SMALL_SEGMENT_MERGING,
                SAMPLE_NAME);
        final SegmentedGenome model = new SegmentedGenome(SEGMENT_FILE_SMALL_SEGMENT_MERGING_NO_SMALL, genome);

        final List<SimpleInterval> expected = model.getSegments();

        final SegmentedGenome modelMerged = model.mergeSmallSegments(SMALL_SEGMENT_TARGET_NUMBER_THRESHOLD);

        final List<SimpleInterval> result = modelMerged.getSegments();

        Assert.assertEquals(expected, result);
    }

    //returns the number of segments that contain a number of targets below a given threshold.
    private static int countSmallSegments(final SegmentedGenome model,
                                          final int targetNumberThreshold) {
        return (int) model.getSegments().stream()
                .filter(s -> model.getGenome().getTargets().targetCount(s) < targetNumberThreshold)
                .count();
    }
}
