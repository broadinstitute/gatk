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
 * Unit tests for {@link SegmentedModel}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class SegmentedModelUnitTest extends BaseTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/tools/exome/";

    private static final File TARGET_FILE_SMALL_SEGMENT_MERGING
            = new File(TEST_SUB_DIR + "targets-small-segment-merging-base.tsv");
    private static final File SNP_FILE_SMALL_SEGMENT_MERGING
            = new File(TEST_SUB_DIR + "snps-small-segment-merging-base.tsv");
    private static final File SEGMENT_FILE_SMALL_SEGMENT_MERGING
            = new File(TEST_SUB_DIR + "segments-small-segment-merging-base.tsv");
    private static final File SEGMENT_FILE_SMALL_SEGMENT_MERGING_NO_SMALL
            = new File(TEST_SUB_DIR + "segments-small-segment-merging-no-small.tsv");
    private static final String SAMPLE_NAME = "TCGA-02-0001-01C-01D-0182-01";

    //a segment is small if number of targets it contains is strictly less than SMALL_SEGMENT_TARGET_NUMBER_THRESHOLD
    private static final int SMALL_SEGMENT_TARGET_NUMBER_THRESHOLD = 3;

    private static final SimpleInterval segment1 = new SimpleInterval("chr1", 1, 4);
    private static final SimpleInterval segment2 = new SimpleInterval("chr1", 5, 12);
    private static final SimpleInterval segment3 = new SimpleInterval("chr1", 11, 20);
    private static final SimpleInterval segment4 = new SimpleInterval("chr2", 1, 10);
    private static final SimpleInterval segment5 = new SimpleInterval("chr2", 2, 9);

    @Test
    public void testMergeSegments() {
        Assert.assertEquals(SegmentMergeUtils.mergeSegments(segment1, segment2), new SimpleInterval("chr1", 1, 12));
        Assert.assertEquals(SegmentMergeUtils.mergeSegments(segment2, segment1), new SimpleInterval("chr1", 1, 12));
        Assert.assertEquals(SegmentMergeUtils.mergeSegments(segment2, segment3), new SimpleInterval("chr1", 5, 20));
        Assert.assertEquals(SegmentMergeUtils.mergeSegments(segment5, segment4), new SimpleInterval("chr2", 1, 10));
    }

    @Test(expectedExceptions = RuntimeException.class)
    public void testJoinWithRuntimeException() {
        SegmentMergeUtils.mergeSegments(segment1, segment4);
    }

    //the base test case starts with identical sets of 10 segments on chromosomes 1 & 3 (with 7 small segments each)
    //and a single small segment (containing only a single target) on chromosome 2
    @Test
    public void testCountSmallSegments() {
        final Genome genome = new Genome(TARGET_FILE_SMALL_SEGMENT_MERGING, SNP_FILE_SMALL_SEGMENT_MERGING,
                SAMPLE_NAME);
        final SegmentedModel model = new SegmentedModel(SEGMENT_FILE_SMALL_SEGMENT_MERGING, genome);

        final int resultNumberOfSegmentsBeforeMerging = model.getSegments().size();
        final int resultNumberOfSmallSegmentsBeforeMerging =
                model.countSmallSegments(SMALL_SEGMENT_TARGET_NUMBER_THRESHOLD);
        final int expectedNumberOfSegmentsBeforeMerging = 21;
        final int expectedNumberOfSmallSegmentsBeforeMerging = 15;
        Assert.assertEquals(resultNumberOfSegmentsBeforeMerging, expectedNumberOfSegmentsBeforeMerging);
        Assert.assertEquals(resultNumberOfSmallSegmentsBeforeMerging, expectedNumberOfSmallSegmentsBeforeMerging);
    }

    @Test
    public void testSmallSegmentMergingBase() {
        final Genome genome = new Genome(TARGET_FILE_SMALL_SEGMENT_MERGING, SNP_FILE_SMALL_SEGMENT_MERGING,
                SAMPLE_NAME);
        final SegmentedModel model = new SegmentedModel(SEGMENT_FILE_SMALL_SEGMENT_MERGING, genome);

        final SegmentedModel modelMerged = model.mergeSmallSegments(SMALL_SEGMENT_TARGET_NUMBER_THRESHOLD);

        final List<SimpleInterval> result = modelMerged.getSegments();

        final int resultNumberOfSegmentsAfterMerging = modelMerged.getSegments().size();
        final int resultNumberOfSmallSegmentsAfterMerging =
                modelMerged.countSmallSegments(SMALL_SEGMENT_TARGET_NUMBER_THRESHOLD);
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
        final SegmentedModel model = new SegmentedModel(SEGMENT_FILE_SMALL_SEGMENT_MERGING_NO_SMALL, genome);

        final List<SimpleInterval> expected = model.getSegments();

        final SegmentedModel modelMerged = model.mergeSmallSegments(SMALL_SEGMENT_TARGET_NUMBER_THRESHOLD);

        final List<SimpleInterval> result = modelMerged.getSegments();

        Assert.assertEquals(expected, result);
    }
}
