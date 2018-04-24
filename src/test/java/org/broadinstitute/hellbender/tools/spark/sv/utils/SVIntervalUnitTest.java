package org.broadinstitute.hellbender.tools.spark.sv.utils;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

public class SVIntervalUnitTest extends GATKBaseTest {
    @Test(groups = "sv")
    public void testOverlapLen() {
        SVInterval container = new SVInterval(1, 1000, 2000);
        SVInterval containee = new SVInterval(1, 1100, 1200);
        Assert.assertTrue(container.overlaps(containee));
        Assert.assertTrue(container.overlapLen(containee) == containee.getLength());

    }
}