package org.broadinstitute.hellbender.tools.exome;


import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;


/**
 * Tests for Segment and SegmentUtils classes
 *
 * @author David Benjamin
 */
public final class SegmentUnitTest extends BaseTest{

    private static boolean sameSegment(final Segment segment1, final Segment segment2) {
        if (segment1 == null || segment2 == null) {
            return false;
        } else {
            return segment1.getSample().equals(segment2.getSample()) && segment1.getInterval().equals(segment2.getInterval());
        }
    }

    //a common set of targets for Segment tests
    private static final TargetCoverage target1 = new TargetCoverage("myTarget", new SimpleInterval("chr1", 1, 2), 1.0);
    private static final TargetCoverage target2 = new TargetCoverage("myTarget", new SimpleInterval("chr1", 3, 4), 2.0);
    private static final TargetCoverage target3 = new TargetCoverage("myTarget", new SimpleInterval("chr1", 5, 6), 3.0);
    private static final TargetCoverage target4 = new TargetCoverage("myTarget", new SimpleInterval("chr2", 1, 2), 4.0);
    private static final TargetCoverage target5 = new TargetCoverage("myTarget", new SimpleInterval("chr2", 3, 4), 5.0);
    private static final TargetCoverage target6 = new TargetCoverage("myTarget", new SimpleInterval("chr2", 5, 6), 6.0);
    private static final HashedListExonCollection<TargetCoverage> targets =
            new HashedListExonCollection<TargetCoverage>(Arrays.asList(target1, target2, target3, target4, target5, target6));

    //a common set of segments for Segment tests
    private final Segment segment1 = new Segment("sample", new SimpleInterval("chr1", 1, 4));
    private final Segment segment2 = new Segment("sample", new SimpleInterval("chr1", 5, 12));
    private final Segment segment3 = new Segment("sample", new SimpleInterval("chr1", 11, 20));
    private final Segment segment4 = new Segment("sample", new SimpleInterval("chr2", 1, 10));

    @Test
    public void testSegmentConstructor() {
        SimpleInterval interval = new SimpleInterval("chr",1,2);
        String call = "call";

        Segment seg1 = new Segment("sample", interval);
        Assert.assertEquals(seg1.getInterval(), interval);

        Segment seg2 = new Segment("sample", interval, call);
        Assert.assertEquals(seg2.getCall(), call);
    }

    @Test
    public void testLocatableMethods() {
        Assert.assertEquals(segment1.getContig(), "chr1");
        Assert.assertEquals(segment1.getStart(), 1);
        Assert.assertEquals(segment1.getEnd(), 4);

        Assert.assertEquals(segment4.getContig(), "chr2");
        Assert.assertEquals(segment4.getStart(), 1);
        Assert.assertEquals(segment4.getEnd(), 10);
    }



    @Test
    public void testOverlapsSegment() {
        Assert.assertTrue(segment2.overlaps(segment3));
        Assert.assertFalse(segment1.overlaps(segment2));
        Assert.assertFalse(segment1.overlaps(segment3));
        Assert.assertFalse(segment1.overlaps(segment4));
    }


    @Test
    public void testSegmentMean() {
        Assert.assertEquals(SegmentUtils.meanTargetCoverage(segment1, targets), (1.0 + 2.0) / 2, 0.00001);
        Assert.assertEquals(SegmentUtils.meanTargetCoverage(segment2, targets), (3.0) / 1, 0.00001);
        Assert.assertEquals(SegmentUtils.meanTargetCoverage(segment4, targets), (4.0 + 5.0 + 6.0) / 3, 0.00001);
    }

    @Test
    public void testSegmentNumTargets() {
        Assert.assertEquals(SegmentUtils.numTargets(segment1, targets), 2);
        Assert.assertEquals(SegmentUtils.numTargets(segment2, targets), 1);
        Assert.assertEquals(SegmentUtils.numTargets(segment3, targets), 0);
        Assert.assertEquals(SegmentUtils.numTargets(segment4, targets), 3);
    }

    @Test
    public void testReadAndWriteSegments() {
        final File file = createTempFile("test",".txt");
        Segment segment1 = new Segment("sample", new SimpleInterval("chr1", 1, 10), "+");
        Segment segment2 = new Segment("sample", new SimpleInterval("chr2", 1, 10), "-");
        List<Segment> segments = Arrays.asList(segment1, segment2);

        try {
            SegmentUtils.writeCalledSegments(file, segments);
            List<Segment> sameSegments = SegmentUtils.readCalledSegments(file);
            Assert.assertEquals(segments.size(), sameSegments.size());
            Assert.assertTrue(sameSegment(segments.get(0), sameSegments.get(0)));
            Assert.assertTrue(sameSegment(segments.get(1), sameSegments.get(1)));
        } catch (IOException e) {
            Assert.fail(e.getMessage());
        }

    }

    @Test
    public void testReadUncalledSegments() {

        final File TEST_DIR = new File("src/test/resources/org/broadinstitute/tools/exome/caller");
        final File TEST_SEGMENTS = new File(TEST_DIR,"segments.tsv");

        try {
            List<Segment> segments = SegmentUtils.readUncalledSegments(TEST_SEGMENTS);
            Assert.assertEquals(segments.size(), 4);
            Assert.assertEquals(segments.get(0).getInterval().getContig(), "chr");
            Assert.assertEquals(segments.get(1).getInterval().getStart(), 300);
            Assert.assertEquals(segments.get(2).getInterval().getEnd(), 550);

        } catch (IOException e) {
            Assert.fail(e.getMessage());
        }

    }

}
