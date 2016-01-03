package org.broadinstitute.hellbender.tools.exome;


import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;


/**
 * Tests for Segment and SegmentUtils classes
 *
 */
public final class ModeledSegmentUnitTest extends BaseTest{

    //a common set of CalledIntervals for tests
    private final ModeledSegment ci1 = new ModeledSegment(new SimpleInterval("chr1", 1, 4), "call", 12, 0);
    private final ModeledSegment ci2 = new ModeledSegment(new SimpleInterval("chr1", 5, 12), 100, 0);
    private final ModeledSegment ci3 = new ModeledSegment(new SimpleInterval("chr2", 1, 10), 100, -0.0002);

    @Test
    public void testSegmentConstructor() {
        SimpleInterval interval = new SimpleInterval("chr",1,2);
        String call = "call";
        double segmentMean = -0.0001;
        final long originalProbeCount = 100;

        final ModeledSegment ci = new ModeledSegment(interval, call, originalProbeCount, segmentMean);
        Assert.assertEquals(ci.getSimpleInterval(), interval);
        Assert.assertEquals(ci.getCall(), call);
        Assert.assertEquals(ci.getTargetCount(), originalProbeCount);
        Assert.assertEquals(ci.getSegmentMean(), segmentMean);
        Assert.assertEquals(ci.getSegmentMeanInCRSpace(), Math.pow(2, segmentMean));

        ci.setTargetCount(5);
        Assert.assertEquals(ci.getTargetCount(), 5);

        double segmentMeanInCRSpace = 2;
        ci.setSegmentMeanInCRSpace(segmentMeanInCRSpace);
        Assert.assertEquals(ci.getSegmentMeanInCRSpace(), segmentMeanInCRSpace);
        Assert.assertEquals(ci.getSegmentMean(), 1.0);

        double segmentMeanAgain = 1;
        ci.setSegmentMean(segmentMeanAgain);
        Assert.assertEquals(ci.getSegmentMeanInCRSpace(), segmentMeanInCRSpace);
        Assert.assertEquals(ci.getSegmentMean(), 1.0);
    }

    @Test
    public void testSetInterval() {
        SimpleInterval interval = new SimpleInterval("chr",1,2);
        String call = "call";
        double segmentMean = -0.0001;
        final long originalProbeCount = 100;
        final ModeledSegment ci = new ModeledSegment(interval, call, originalProbeCount, segmentMean);

        SimpleInterval updated = new SimpleInterval("chr",1,3);
        ci.setSimpleInterval(updated);
        Assert.assertEquals(updated.getEnd(), ci.getEnd());
    }

    @Test
    public void testLocatableMethods() {
        Assert.assertEquals(ci1.getContig(), "chr1");
        Assert.assertEquals(ci1.getStart(), 1);
        Assert.assertEquals(ci1.getEnd(), 4);

        Assert.assertEquals(ci3.getContig(), "chr2");
        Assert.assertEquals(ci3.getStart(), 1);
        Assert.assertEquals(ci3.getEnd(), 10);
    }

    @Test
    public void testEquals() {
        Assert.assertFalse(ci1.equals(ci2));
        Assert.assertFalse(ci1.equals(ci3));
        Assert.assertEquals(ci1, ci1);
        Assert.assertEquals(ci1, new ModeledSegment(new SimpleInterval("chr1", 1, 4), "call", 12, 0));
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testOriginalProbeCountInvalid() {
        SimpleInterval interval = new SimpleInterval("chr",1,2);
        String call = "call";
        double segmentMean = -0.0001;
        final long originalProbeCount = -100;

        final ModeledSegment ci = new ModeledSegment(interval, call, originalProbeCount, segmentMean);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testOriginalProbeCountInvalidOnSet() {
        SimpleInterval interval = new SimpleInterval("chr",1,2);
        String call = "call";
        double segmentMean = -0.0001;
        final long originalProbeCount = 100;

        final ModeledSegment ci = new ModeledSegment(interval, call, originalProbeCount, segmentMean);
        ci.setTargetCount(-100);
    }
}
