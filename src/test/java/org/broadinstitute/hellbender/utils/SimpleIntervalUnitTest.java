package org.broadinstitute.hellbender.utils;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public final class SimpleIntervalUnitTest extends BaseTest {

    @DataProvider(name = "badIntervals")
    public Object[][] badIntervals(){
        return new Object[][]{
                {null,1,12, "null contig"},
                {"1", 0, 10, "start==0"},
                {"1", -10, 10, "negative start"},
                {"1", 10, 9, "end < start"}
        };
    }

    @Test(dataProvider = "badIntervals", expectedExceptions = IllegalArgumentException.class)
    public void badIntervals(String contig, int start, int end, String name){
        SimpleInterval interval = new SimpleInterval(contig, start, end);
    }

    @DataProvider(name = "goodIntervals")
    public Object[][] goodIntervals(){
        return new Object[][]{
                {"1", "1", 1, Integer.MAX_VALUE},
                {"1:2", "1", 2, 2},
                {"1:3+", "1", 3, Integer.MAX_VALUE},
                {"1:4-5", "1", 4, 5},
                {"1:2,000", "1", 2000, 2000},
                {"1:3,000+", "1", 3000, Integer.MAX_VALUE},
                {"1:4,000-5,000", "1", 4000, 5000},
                {"1:4,0,0,0-5,0,0,0", "1", 4000, 5000}, //this is OK too, we just remove commas wherever they are
        };
    }

    @Test(dataProvider = "goodIntervals")
    public void testGoodIntervals(String str, String contig, int start, int end){
        SimpleInterval interval = new SimpleInterval(str);
        Assert.assertEquals(interval.getContig(), contig, "contig");
        Assert.assertEquals(interval.getStart(), start, "start");
        Assert.assertEquals(interval.getEnd(), end, "end");
    }

    @Test
    public void testEquality(){
        final SimpleInterval i1 = new SimpleInterval("1",1,100);
        final SimpleInterval i2 = new SimpleInterval("1",1,100);
        final SimpleInterval i3 = new SimpleInterval("chr1", 1, 100);
        final SimpleInterval i4 = new SimpleInterval("chr1", 2, 100);
        final SimpleInterval i5 = new SimpleInterval("chr1", 1, 200);
        final SimpleInterval i6 = new SimpleInterval("1:1-100");

        Assert.assertTrue(i1.equals(i1));
        Assert.assertTrue(i1.equals(i2));
        Assert.assertTrue(i2.equals(i1));
        Assert.assertTrue(i1.equals(i6));
        Assert.assertTrue(i6.equals(i1));

        Assert.assertFalse(i1.equals(i3));
        Assert.assertFalse(i1.equals(i4));
        Assert.assertFalse(i1.equals(i5));
    }

    @Test
    public void testToString(){
        final SimpleInterval i1 = new SimpleInterval("1",1,100);
        Assert.assertEquals(i1.toString(), "1:1-100");
    }

    @DataProvider(name = "IntervalSizeData")
    public Object[][] getIntervalSizeData() {
        // Intervals + expected sizes
        return new Object[][]{
                { new SimpleInterval("1", 1, 1), 1 },
                { new SimpleInterval("1", 1, 2), 2 },
                { new SimpleInterval("1", 1, 10), 10 },
                { new SimpleInterval("1", 2, 10), 9 }
        };
    }

    @Test(dataProvider = "IntervalSizeData")
    public void testGetSize( final SimpleInterval interval, final int expectedSize ) {
        Assert.assertEquals(interval.size(), expectedSize, "size() incorrect for interval " + interval);
    }

    @DataProvider(name = "IntervalOverlapData")
    public Object[][] getIntervalOverlapData() {
        final SimpleInterval standardInterval = new SimpleInterval("1", 10, 20);
        final SimpleInterval oneBaseInterval = new SimpleInterval("1", 10, 10);

        return new Object[][] {
                { standardInterval, new SimpleInterval("2", 10, 20), false },
                { standardInterval, new SimpleInterval("1", 1, 5), false },
                { standardInterval, new SimpleInterval("1", 1, 9), false },
                { standardInterval, new SimpleInterval("1", 1, 10), true },
                { standardInterval, new SimpleInterval("1", 1, 15), true },
                { standardInterval, new SimpleInterval("1", 10, 10), true },
                { standardInterval, new SimpleInterval("1", 10, 15), true },
                { standardInterval, new SimpleInterval("1", 10, 20), true },
                { standardInterval, new SimpleInterval("1", 15, 20), true },
                { standardInterval, new SimpleInterval("1", 15, 25), true },
                { standardInterval, new SimpleInterval("1", 20, 20), true },
                { standardInterval, new SimpleInterval("1", 20, 25), true },
                { standardInterval, new SimpleInterval("1", 21, 25), false },
                { standardInterval, new SimpleInterval("1", 25, 30), false },
                { oneBaseInterval, new SimpleInterval("2", 10, 10), false },
                { oneBaseInterval, new SimpleInterval("1", 1, 5), false },
                { oneBaseInterval, new SimpleInterval("1", 1, 9), false },
                { oneBaseInterval, new SimpleInterval("1", 1, 10), true },
                { oneBaseInterval, new SimpleInterval("1", 10, 10), true },
                { oneBaseInterval, new SimpleInterval("1", 10, 15), true },
                { oneBaseInterval, new SimpleInterval("1", 11, 15), false },
                { oneBaseInterval, new SimpleInterval("1", 15, 20), false },
                { standardInterval, null, false },
                { standardInterval, standardInterval, true }
        };
    }

    @Test(dataProvider = "IntervalOverlapData")
    public void testOverlap( final SimpleInterval firstInterval, final SimpleInterval secondInterval, final boolean expectedOverlapResult ) {
        Assert.assertEquals(firstInterval.overlaps(secondInterval), expectedOverlapResult,
                            "overlap() returned incorrect result for intervals " + firstInterval + " and " + secondInterval);
    }

    @DataProvider(name = "IntervalContainsData")
    public Object[][] getIntervalContainsData() {
        final SimpleInterval containingInterval = new SimpleInterval("1", 10, 20);

        return new Object[][] {
                { containingInterval, new SimpleInterval("2", 10, 20), false },
                { containingInterval, new SimpleInterval("1", 1, 5), false },
                { containingInterval, new SimpleInterval("1", 1, 10), false },
                { containingInterval, new SimpleInterval("1", 5, 15), false },
                { containingInterval, new SimpleInterval("1", 9, 10), false },
                { containingInterval, new SimpleInterval("1", 9, 20), false },
                { containingInterval, new SimpleInterval("1", 10, 10), true },
                { containingInterval, new SimpleInterval("1", 10, 15), true },
                { containingInterval, new SimpleInterval("1", 10, 20), true },
                { containingInterval, new SimpleInterval("1", 10, 21), false },
                { containingInterval, new SimpleInterval("1", 15, 25), false },
                { containingInterval, new SimpleInterval("1", 20, 20), true },
                { containingInterval, new SimpleInterval("1", 20, 21), false },
                { containingInterval, new SimpleInterval("1", 20, 25), false },
                { containingInterval, new SimpleInterval("1", 21, 25), false },
                { containingInterval, new SimpleInterval("1", 25, 30), false },
                { containingInterval, null, false },
                { containingInterval, containingInterval, true }
        };
    }

    @Test(dataProvider = "IntervalContainsData")
    public void testContains( final SimpleInterval firstInterval, final SimpleInterval secondInterval, final boolean expectedContainsResult ) {
        Assert.assertEquals(firstInterval.contains(secondInterval), expectedContainsResult,
                            "contains() returned incorrect result for intervals " + firstInterval + " and " + secondInterval);
    }
}
