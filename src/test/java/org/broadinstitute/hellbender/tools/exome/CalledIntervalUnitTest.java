package org.broadinstitute.hellbender.tools.exome;


import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;


/**
 * Tests for Segment and SegmentUtils classes
 *
 * @author David Benjamin
 */
public final class CalledIntervalUnitTest extends BaseTest{

    //a common set of CalledIntervals for tests
    private final CalledInterval ci1 = new CalledInterval(new SimpleInterval("chr1", 1, 4), "call");
    private final CalledInterval ci2 = new CalledInterval(new SimpleInterval("chr1", 5, 12), "call");
    private final CalledInterval ci3 = new CalledInterval(new SimpleInterval("chr2", 1, 10), "call");

    @Test
    public void testSegmentConstructor() {
        SimpleInterval interval = new SimpleInterval("chr",1,2);
        String call = "call";

        final CalledInterval ci = new CalledInterval(interval, call);
        Assert.assertEquals(ci.getInterval(), interval);
        Assert.assertEquals(ci.getCall(), call);
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
        Assert.assertEquals(ci1, new CalledInterval(new SimpleInterval("chr1", 1, 4), "call"));
    }

}
