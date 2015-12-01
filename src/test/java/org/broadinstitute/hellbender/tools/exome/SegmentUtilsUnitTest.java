package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
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


        List<SimpleInterval> segments = SegmentUtils.readIntervalsFromSegfile(TEST_SEGMENTS);
        Assert.assertEquals(segments.size(), 4);
        Assert.assertEquals(segments.get(0).getContig(), "chr");
        Assert.assertEquals(segments.get(1).getStart(), 300);
        Assert.assertEquals(segments.get(2).getEnd(), 550);
    }

    @Test
    public void testReadAndWriteCalledIntervals() {
        final File file = createTempFile("test",".txt");
        List<ModeledSegment> calledIntervals = Arrays.asList(ci1, ci2, ci3);

        SegmentUtils.writeModeledSegmentsToSegfile(file, calledIntervals, "sample");
        List<ModeledSegment> sameIntervals = SegmentUtils.readModeledSegmentsFromSegfile(file);
        Assert.assertEquals(calledIntervals, sameIntervals);
    }

    @Test
    public void testMeanTargetCoverage() {
        //a common set of targets for tests
        final TargetCoverage target1 = new TargetCoverage("myTarget", new SimpleInterval("chr1", 1, 2), 1.0);
        final TargetCoverage target2 = new TargetCoverage("myTarget", new SimpleInterval("chr1", 3, 4), 2.0);
        final TargetCoverage target3 = new TargetCoverage("myTarget", new SimpleInterval("chr1", 5, 6), 3.0);
        final TargetCoverage target4 = new TargetCoverage("myTarget", new SimpleInterval("chr2", 1, 2), 4.0);
        final TargetCoverage target5 = new TargetCoverage("myTarget", new SimpleInterval("chr2", 3, 4), 5.0);
        final TargetCoverage target6 = new TargetCoverage("myTarget", new SimpleInterval("chr2", 5, 6), 6.0);
        final HashedListTargetCollection<TargetCoverage> targets =
                new HashedListTargetCollection<>(Arrays.asList(target1, target2, target3, target4, target5, target6));

        Assert.assertEquals(SegmentUtils.meanTargetCoverage(ci1, targets), (1.0 + 2.0) / 2, 0.00001);
        Assert.assertEquals(SegmentUtils.meanTargetCoverage(ci2, targets), (3.0) / 1, 0.00001);
        Assert.assertEquals(SegmentUtils.meanTargetCoverage(ci3, targets), (4.0 + 5.0 + 6.0) / 3, 0.00001);
    }

    @Test
    public void testSegmentDifference() {
        //a common set of targets for tests
        final TargetCoverage target1 = new TargetCoverage("myTarget", new SimpleInterval("chr1", 1, 2), 1.0);
        final TargetCoverage target2 = new TargetCoverage("myTarget", new SimpleInterval("chr1", 3, 4), 2.0);
        final TargetCoverage target3 = new TargetCoverage("myTarget", new SimpleInterval("chr2", 3, 4), 5.0);
        final TargetCoverage target4 = new TargetCoverage("myTarget", new SimpleInterval("chr2", 5, 6), 6.0);
        final HashedListTargetCollection<TargetCoverage> targets =
                new HashedListTargetCollection<>(Arrays.asList(target1, target2, target3, target4));

        final ModeledSegment seg = new ModeledSegment(new SimpleInterval("chr1", 1, 4), "", 200, 1.1);

        double [] gt_diffs = {-0.14354693, 1.85645307};
        List<Double> diffs = SegmentUtils.segmentMeanTargetDifference(seg, targets);

        for (int i = 0; i < diffs.size(); i++){
            Assert.assertEquals(gt_diffs[i], diffs.get(i), 0.0000001, "Predicted difference not the same as ground truth (within tolerance).");
        }
    }

    @Test
    public void testParsimoniousInterval() {
        final SimpleInterval target1 = new SimpleInterval("chr1", 1, 10);
        final SimpleInterval target2 = new SimpleInterval("chr1", 20, 30);
        final SimpleInterval target3 = new SimpleInterval("chr1", 31, 40);
        final SimpleInterval target4 = new SimpleInterval("chr1", 90, 100);
        final HashedListTargetCollection<SimpleInterval> targets =
                new HashedListTargetCollection<>(Arrays.asList(target1, target2, target3, target4));

        final SimpleInterval snp1 = new SimpleInterval("chr1", 5, 5);
        final SimpleInterval snp2 = new SimpleInterval("chr1", 42, 42);
        final SimpleInterval snp3 = new SimpleInterval("chr1", 60, 60);
        final HashedListTargetCollection<SimpleInterval> snps =
                new HashedListTargetCollection<>(Arrays.asList(snp1, snp2, snp3));

        Assert.assertEquals(SegmentUtils.parsimoniousInterval(new SimpleInterval("chr1", 1, 10), targets, snps),
                new SimpleInterval("chr1",1,10));
        Assert.assertEquals(SegmentUtils.parsimoniousInterval(new SimpleInterval("chr1", 1, 4), targets, snps),
                new SimpleInterval("chr1",1,10));
        Assert.assertEquals(SegmentUtils.parsimoniousInterval(new SimpleInterval("chr1", 2, 4), targets, snps),
                new SimpleInterval("chr1",1,10));
        Assert.assertEquals(SegmentUtils.parsimoniousInterval(new SimpleInterval("chr1", 10, 25), targets, snps),
                new SimpleInterval("chr1",1,30));
        Assert.assertEquals(SegmentUtils.parsimoniousInterval(new SimpleInterval("chr1", 11, 31), targets, snps),
                new SimpleInterval("chr1",20,40));
        Assert.assertEquals(SegmentUtils.parsimoniousInterval(new SimpleInterval("chr1", 40, 42), targets, snps),
                new SimpleInterval("chr1",31,42));
        Assert.assertEquals(SegmentUtils.parsimoniousInterval(new SimpleInterval("chr1", 43, 89), targets, snps),
                new SimpleInterval("chr1",60,60));
        Assert.assertEquals(SegmentUtils.parsimoniousInterval(new SimpleInterval("chr1", 101, 200), targets, snps),
                new SimpleInterval("chr1",101,101));
        Assert.assertEquals(SegmentUtils.parsimoniousInterval(new SimpleInterval("chr2", 1, 10), targets, snps),
                new SimpleInterval("chr2",1,1));
    }

    @Test
    public void testSegmentUnion() {
        final SimpleInterval target1 = new SimpleInterval("chr1", 1, 10);
        final SimpleInterval target2 = new SimpleInterval("chr1", 20, 30);
        final SimpleInterval target3 = new SimpleInterval("chr1", 31, 40);
        final SimpleInterval target4 = new SimpleInterval("chr1", 90, 100);
        final HashedListTargetCollection<SimpleInterval> targets =
                new HashedListTargetCollection<>(Arrays.asList(target1, target2, target3, target4));

        final SimpleInterval snp1 = new SimpleInterval("chr1", 5, 5);
        final SimpleInterval snp2 = new SimpleInterval("chr1", 42, 42);
        final SimpleInterval snp3 = new SimpleInterval("chr1", 60, 60);
        final SimpleInterval snp4 = new SimpleInterval("chr2", 10, 10);
        final HashedListTargetCollection<SimpleInterval> snps =
                new HashedListTargetCollection<>(Arrays.asList(snp1, snp2, snp3, snp4));

        final List<SimpleInterval> targetSegments = Arrays.asList(
                new SimpleInterval("chr1", 1, 10),
                new SimpleInterval("chr1", 20, 40),
                new SimpleInterval("chr1", 42, 60),
                new SimpleInterval("chr1", 50, 60),
                new SimpleInterval("chr2", 5, 15));

        final List<SimpleInterval> snpSegments = Arrays.asList(
                new SimpleInterval("chr1", 4, 5),
                new SimpleInterval("chr1", 5, 6),
                new SimpleInterval("chr1", 80, 110),
                new SimpleInterval("chr1", 200, 200),
                new SimpleInterval("chr2", 5, 15));


        /**
         * Expected behavior: on chr1 we get breakpoints 1,4,5,6,10,20,40,42,50,60,80,110,200.  We now follow the logic of
         * SegmentUtils.segmentUnion.
         * [1,4] -> [1,10] under parsimoniousInterval, so the next position to start is 11.
         * The "intervals" [11,5]; [11,6]; and [11,10] end before they start and are trivially empty.
         * [11,20] -> [20,30] under parsimoniousInterval, so the next position to start is 31.
         * [31,40] -> [31,40] and the next start is 41.
         * [41,42] -> [42,42] and the next start is 43.
         * [43,50] is empty.
         * [43,60] -> [60,60] and the next start is 61.
         * [61,80] is empty.
         * [61,110] -> [90,100].
         * [101,200] is empty.
         * Thus we get [1,10]; [20,30]; [42,42]; [60,60]; [90,100] on chr1
         * More simply, we get [10,10] on chr2.
         */

        List<SimpleInterval> segUnion = SegmentUtils.segmentUnion(targetSegments, snpSegments, targets, snps);
        List<SimpleInterval> expected = Arrays.asList(
                new SimpleInterval("chr1", 1, 10),
                new SimpleInterval("chr1", 20, 30),
                new SimpleInterval("chr1", 31, 40),
                new SimpleInterval("chr1", 42, 42),
                new SimpleInterval("chr1", 60, 60),
                new SimpleInterval("chr1", 90, 100),
                new SimpleInterval("chr2", 10, 10));
        Assert.assertEquals(segUnion, expected);
    }



}
