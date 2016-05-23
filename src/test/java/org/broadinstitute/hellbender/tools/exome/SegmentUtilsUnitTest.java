package org.broadinstitute.hellbender.tools.exome;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


public class SegmentUtilsUnitTest extends BaseTest {
    //a common set of ModeledSegments for tests
    private static final ModeledSegment ci1 = new ModeledSegment(new SimpleInterval("chr1", 1, 4), "call", 100, 1);
    private static final ModeledSegment ci2 = new ModeledSegment(new SimpleInterval("chr1", 5, 12), "call", 200, 0);
    private static final ModeledSegment ci3 = new ModeledSegment(new SimpleInterval("chr2", 1, 10), "call", 300, 0.5);

    @Test
    public void testReadUncalledSegments() {
        final File TEST_DIR = new File("src/test/resources/org/broadinstitute/hellbender/tools/exome/caller");
        final File TEST_SEGMENTS = new File(TEST_DIR,"segments.tsv");

        List<SimpleInterval> segments = SegmentUtils.readIntervalsFromSegmentFile(TEST_SEGMENTS);
        Assert.assertEquals(segments.size(), 4);
        Assert.assertEquals(segments.get(0).getContig(), "chr");
        Assert.assertEquals(segments.get(1).getStart(), 300);
        Assert.assertEquals(segments.get(2).getEnd(), 550);
    }

    @Test
    public void testReadAndWriteCalledIntervals() {
        final File file = createTempFile("test",".txt");
        List<ModeledSegment> calledIntervals = Arrays.asList(ci1, ci2, ci3);

        SegmentUtils.writeModeledSegmentFile(file, calledIntervals, "sample");
        List<ModeledSegment> sameIntervals = SegmentUtils.readModeledSegmentsFromSegmentFile(file);
        Assert.assertEquals(calledIntervals, sameIntervals);
    }

    @Test
    public void testSegmentDifference() {
        final ReadCountRecord.SingleSampleRecord target1 =
                new ReadCountRecord.SingleSampleRecord(new Target("myTarget", new SimpleInterval("chr1", 1, 2)), 1.0);
        final ReadCountRecord.SingleSampleRecord target2 =
                new ReadCountRecord.SingleSampleRecord(new Target("myTarget", new SimpleInterval("chr1", 3, 4)), 2.0);
        final ReadCountRecord.SingleSampleRecord target3 =
                new ReadCountRecord.SingleSampleRecord(new Target("myTarget", new SimpleInterval("chr2", 3, 4)), 5.0);
        final ReadCountRecord.SingleSampleRecord target4 =
                new ReadCountRecord.SingleSampleRecord(new Target("myTarget", new SimpleInterval("chr2", 5, 6)), 6.0);
        final HashedListTargetCollection<ReadCountRecord.SingleSampleRecord> targets =
                new HashedListTargetCollection<>(Arrays.asList(target1, target2, target3, target4));

        final ModeledSegment seg = new ModeledSegment(new SimpleInterval("chr1", 1, 4), "", 200, 1.1);

        double [] gt_diffs = {-0.14354693, 1.85645307};
        List<Double> diffs = SegmentUtils.segmentMeanTargetDifference(seg, targets);

        for (int i = 0; i < diffs.size(); i++){
            Assert.assertEquals(gt_diffs[i], diffs.get(i), 0.0000001, "Predicted difference not the same as ground truth (within tolerance).");
        }
    }

    @Test
    public void testTrimInterval() {
        final SimpleInterval target1 = new SimpleInterval("chr1", 1, 10);
        final SimpleInterval target2 = new SimpleInterval("chr1", 20, 30);
        final SimpleInterval target3 = new SimpleInterval("chr1", 31, 40);
        final SimpleInterval target4 = new SimpleInterval("chr1", 90, 100);
        final SimpleInterval target5 = new SimpleInterval("chr1", 110, 120);
        final HashedListTargetCollection<SimpleInterval> targets =
                new HashedListTargetCollection<>(Arrays.asList(target1, target2, target3, target4, target5));

        final SimpleInterval snp1 = new SimpleInterval("chr1", 5, 5);
        final SimpleInterval snp2 = new SimpleInterval("chr1", 42, 42);
        final SimpleInterval snp3 = new SimpleInterval("chr1", 60, 60);
        final SimpleInterval snp4 = new SimpleInterval("chr1", 91, 91);
        final HashedListTargetCollection<SimpleInterval> snps =
                new HashedListTargetCollection<>(Arrays.asList(snp1, snp2, snp3, snp4));

        Assert.assertEquals(SegmentUtils.trimInterval(new SimpleInterval("chr1", 1, 10), targets, snps),
                new SimpleInterval("chr1", 1, 10));
        Assert.assertEquals(SegmentUtils.trimInterval(new SimpleInterval("chr1", 2, 4), targets, snps),
                new SimpleInterval("chr1", 2, 4));
        Assert.assertEquals(SegmentUtils.trimInterval(new SimpleInterval("chr1", 5, 8), targets, snps),
                new SimpleInterval("chr1", 5, 8));
        Assert.assertEquals(SegmentUtils.trimInterval(new SimpleInterval("chr1", 8, 25), targets, snps),
                new SimpleInterval("chr1", 8, 25));
        Assert.assertEquals(SegmentUtils.trimInterval(new SimpleInterval("chr1", 11, 31), targets, snps),
                new SimpleInterval("chr1", 20, 31));
        Assert.assertEquals(SegmentUtils.trimInterval(new SimpleInterval("chr1", 15, 35), targets, snps),
                new SimpleInterval("chr1", 20, 35));
        Assert.assertEquals(SegmentUtils.trimInterval(new SimpleInterval("chr1", 15, 45), targets, snps),
                new SimpleInterval("chr1", 20, 42));
        Assert.assertEquals(SegmentUtils.trimInterval(new SimpleInterval("chr1", 40, 42), targets, snps),
                new SimpleInterval("chr1", 40, 42));
        Assert.assertEquals(SegmentUtils.trimInterval(new SimpleInterval("chr1", 43, 89), targets, snps),
                new SimpleInterval("chr1", 60, 60));
        Assert.assertEquals(SegmentUtils.trimInterval(new SimpleInterval("chr1", 90, 114), targets, snps),
                new SimpleInterval("chr1", 90, 114));
        Assert.assertEquals(SegmentUtils.trimInterval(new SimpleInterval("chr1", 101, 200), targets, snps),
                new SimpleInterval("chr1", 110, 120));
        Assert.assertEquals(SegmentUtils.trimInterval(new SimpleInterval("chr2", 1, 10), targets, snps),
                new SimpleInterval("chr2", 1, 10));
    }

    /**
     * Test for {@link SegmentUtils#unionSegments}.  Expected behavior:
     * <p>
     * On chr1 {@link SegmentUtils#collectBreakpointsByContig} gives:
     * </p>
     *      <p>
     *      1, 5, 10, 20, 40, 40, 42, 90, 91, 115, 125, 140.
     *      </p>
     * <p>
     * Then {@link SegmentUtils#constructUntrimmedSegments} finds the segments:
     * </p>
     *      <p>
     *      [1, 4], [5, 10], [11, 19], [20, 40], [41, 41], [42, 89], [90, 91], [92, 114], [115, 125], [126, 140].
     *      </p>
     * <p>
     * and returns the non-empty segments:
     * </p>
     *      <p>
     *      [1, 4], [5, 10], [20, 40], [42, 89], [90, 91], [92, 114], [115, 125], [126, 140].
     *      </p>
     * <p>
     * Then {@link SegmentUtils#mergeSpuriousStartsAndEnds} merges the last segment left to form [115, 140],
     * and {@link SegmentMergeUtils#mergeSpuriousMiddles} randomly merges segment [92, 114] left or right.
     * </p>
     * <p>
     * Finally, {@link SegmentUtils#trimInterval} gives:
     * </p>
     *      <p>
     *      [1, 10], [20, 40], [42, 42], [90, 114], [115, 140] (if [92, 114] merged left) or
     *      </p>
     *      <p>
     *      [1, 10], [20, 40], [42, 42], [90, 91], [92, 140] (if [92, 114] merged right)
     *      </p>
     * <p>
     * The remaining empty segment on chr2 is retained.
     */
    @Test
    public void testUnionSegments() {
        final String sampleName = "placeholder_sample_name";
        final List<Target> targets = new ArrayList<Target>();
        targets.add(new Target("t1", new SimpleInterval("chr1", 1, 10)));
        targets.add(new Target("t2", new SimpleInterval("chr1", 20, 30)));
        targets.add(new Target("t3", new SimpleInterval("chr1", 31, 40)));
        targets.add(new Target("t4", new SimpleInterval("chr1", 90, 100)));
        targets.add(new Target("t5", new SimpleInterval("chr1", 110, 120)));
        targets.add(new Target("t6", new SimpleInterval("chr1", 130, 140)));

        final RealMatrix zeroCoverageMatrix = new Array2DRowRealMatrix(targets.size(), 1);
        final ReadCountCollection counts = new ReadCountCollection(targets, Arrays.asList(sampleName), zeroCoverageMatrix);

        final AllelicCount snp1 = new AllelicCount(new SimpleInterval("chr1", 5, 5), 0, 1);
        final AllelicCount snp2 = new AllelicCount(new SimpleInterval("chr1", 40, 40), 0, 1);
        final AllelicCount snp3 = new AllelicCount(new SimpleInterval("chr1", 42, 42), 0, 1);
        final AllelicCount snp4 = new AllelicCount(new SimpleInterval("chr1", 91, 91), 0, 1);
        final AllelicCount snp5 = new AllelicCount(new SimpleInterval("chr1", 115, 115), 0, 1);
        final AllelicCount snp6 = new AllelicCount(new SimpleInterval("chr1", 125, 125), 0, 1);
        final AllelicCount snp7 = new AllelicCount(new SimpleInterval("chr2", 10, 10), 0, 1);
        final List<AllelicCount> snps =
                Arrays.asList(snp1, snp2, snp3, snp4, snp5, snp6, snp7);

        final List<SimpleInterval> targetSegments = Arrays.asList(
                new SimpleInterval("chr1", 1, 10),
                new SimpleInterval("chr1", 20, 40),
                new SimpleInterval("chr1", 90, 140));

        final List<SimpleInterval> snpSegments = Arrays.asList(
                new SimpleInterval("chr1", 5, 40),
                new SimpleInterval("chr1", 42, 91),
                new SimpleInterval("chr1", 115, 125),
                new SimpleInterval("chr2", 10, 10));

        final List<SimpleInterval> unionedSegments =
                SegmentUtils.unionSegments(targetSegments, snpSegments, new Genome(counts, snps, "test"));
        final List<SimpleInterval> expectedLeft = Arrays.asList(
                new SimpleInterval("chr1", 1, 10),
                new SimpleInterval("chr1", 20, 40),
                new SimpleInterval("chr1", 42, 42),
                new SimpleInterval("chr1", 90, 114),
                new SimpleInterval("chr1", 115, 140),
                new SimpleInterval("chr2", 10, 10));
        final List<SimpleInterval> expectedRight = Arrays.asList(
                new SimpleInterval("chr1", 1, 10),
                new SimpleInterval("chr1", 20, 40),
                new SimpleInterval("chr1", 42, 42),
                new SimpleInterval("chr1", 90, 91),
                new SimpleInterval("chr1", 92, 140),
                new SimpleInterval("chr2", 10, 10));
        Assert.assertTrue(unionedSegments.equals(expectedLeft) || unionedSegments.equals(expectedRight));
    }
}
