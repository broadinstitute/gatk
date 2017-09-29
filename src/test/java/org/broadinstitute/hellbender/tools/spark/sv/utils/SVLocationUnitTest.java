package org.broadinstitute.hellbender.tools.spark.sv.utils;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

public class SVLocationUnitTest extends BaseTest {
    @Test(groups = "sv")
    public void testGetters() {
        final SVLocation location = new SVLocation(1, 2);
        Assert.assertEquals(location.getContig(), 1);
        Assert.assertEquals(location.getPosition(), 2);
    }

    @Test(groups = "sv")
    public void testEqualsAndHashcode() {
        final SVLocation location1 = new SVLocation(1, 2);
        final SVLocation location2 = new SVLocation(1, 2);
        Assert.assertEquals(location1, location2);
        Assert.assertEquals(location1.hashCode(), location2.hashCode());
        final SVLocation location3 = new SVLocation(2, 2);
        Assert.assertNotEquals(location1, location3);
        Assert.assertNotEquals(location1.hashCode(), location3.hashCode());
        final SVLocation location4 = new SVLocation(1, 3);
        Assert.assertNotEquals(location1, location4);
        Assert.assertNotEquals(location1, null);
        Assert.assertNotEquals(location1, this);
    }

    @Test(groups = "sv")
    public void testComparison() {
        final SVLocation location1 = new SVLocation(1, 2);
        final SVLocation location2 = new SVLocation(1, 2);
        Assert.assertEquals(location1.compareTo(location1), 0);
        Assert.assertEquals(location1.compareTo(location2), 0);
        Assert.assertEquals(location2.compareTo(location1), 0);
        final SVLocation location3 = new SVLocation(2, 2);
        Assert.assertEquals(location1.compareTo(location3), -1);
        Assert.assertEquals(location3.compareTo(location1), 1);
        final SVLocation location4 = new SVLocation(1, 3);
        Assert.assertEquals(location1.compareTo(location4), -1);
        Assert.assertEquals(location4.compareTo(location1), 1);
    }
}
