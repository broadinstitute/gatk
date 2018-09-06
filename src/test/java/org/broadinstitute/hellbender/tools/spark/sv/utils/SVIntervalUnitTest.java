package org.broadinstitute.hellbender.tools.spark.sv.utils;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

public class SVIntervalUnitTest extends GATKBaseTest {
    @Test(groups = "sv")
    public void testOverlapLen() {
        SVInterval container = new SVInterval(1, 1000, 2000);
        SVInterval containee = new SVInterval(1, 1100, 1200);
        Assert.assertTrue(container.overlaps(containee));
        Assert.assertTrue(container.overlapLen(containee) == containee.getLength());

    }

    @DataProvider
    private Object[][] forTestContainment() {
        final List<Object[]> data = new ArrayList<>(20);
        SVInterval one = new SVInterval(1, 1000, 2000);
        SVInterval two = new SVInterval(1, 1100, 1200);
        data.add(new Object[]{one, two, true});
        data.add(new Object[]{two, one, false});

        two = new SVInterval(1, 1000, 1200);
        data.add(new Object[]{one, two, true});
        data.add(new Object[]{two, one, false});

        two = new SVInterval(1, 1100, 2000);
        data.add(new Object[]{one, two, true});
        data.add(new Object[]{two, one, false});

        two = new SVInterval(1, 1000, 2000);
        data.add(new Object[]{one, two, true});
        data.add(new Object[]{two, one, true});

        return data.toArray(new Object[data.size()][]);
    }

    @Test(groups = "sv", dataProvider = "forTestContainment")
    public void testContainment(final SVInterval one, final SVInterval two, final boolean expected) {
        Assert.assertEquals(one.contains(two), expected);
    }

    @DataProvider
    private Object[][] forTestCtorArgsValidators() {
        final List<Object[]> data = new ArrayList<>(20);

        data.add(new Object[]{-1,  0,  0, SVInterval.SVIntervalConstructorArgsValidator.RefuseNegativeContigs});
        data.add(new Object[]{ 0, -1,  0, SVInterval.SVIntervalConstructorArgsValidator.RefuseNegativeCoordinates});
        data.add(new Object[]{ 0,  0, -1, SVInterval.SVIntervalConstructorArgsValidator.RefuseNegativeCoordinates});

        data.add(new Object[]{-1,  0,  0, SVInterval.SVIntervalConstructorArgsValidator.RefuseNegativeContigAndCoordinates});
        data.add(new Object[]{ 0, -1,  0, SVInterval.SVIntervalConstructorArgsValidator.RefuseNegativeContigAndCoordinates});
        data.add(new Object[]{ 0,  0, -1, SVInterval.SVIntervalConstructorArgsValidator.RefuseNegativeContigAndCoordinates});
        data.add(new Object[]{-1, -1, -1, SVInterval.SVIntervalConstructorArgsValidator.RefuseNegativeContigAndCoordinates});

        return data.toArray(new Object[data.size()][]);
    }
    @Test(groups = "sv", dataProvider = "forTestCtorArgsValidators", expectedExceptions = IllegalArgumentException.class)
    public void test(final int contig, final int start, final int end, final SVInterval.SVIntervalConstructorArgsValidator argsValidator) {
        new SVInterval(contig, start, end, argsValidator);
    }
}