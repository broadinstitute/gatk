package org.broadinstitute.hellbender.utils;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;


public class ClosedSVIntervalUnitTest extends GATKBaseTest {

    private final String[] contigIDToName = new String[]{"chr1", "chr2"};

    @DataProvider(name="singleIntervals")
    public Object[][] getSingleIntervalData() {
        return new Object[][] {
                { new ClosedSVInterval(0, 1, 100), 100, new SimpleInterval("chr1", 1, 100)},
                { new ClosedSVInterval(1, 201, 201), 1, new SimpleInterval("chr2", 201, 201)},
                { new ClosedSVInterval(2, 2, 3), 2, null}
        };
    }


    @Test(dataProvider = "singleIntervals")
    public void testToSimpleInterval(
            final ClosedSVInterval interval,
            final int expectedLength,
            final SimpleInterval expectedSimpleInterval)
    {
        try {
            Assert.assertEquals(interval.toSimpleInterval(contigIDToName), expectedSimpleInterval);
        } catch (ArrayIndexOutOfBoundsException e) {
            Assert.assertNull(expectedSimpleInterval);
        }
    }

    @Test(dataProvider = "singleIntervals")
    public void testGetLen(
            final ClosedSVInterval interval,
            final int expectedLength,
            final SimpleInterval expectedSimpleInterval)
    {
        Assert.assertEquals(interval.getLength(), expectedLength);
    }

    @DataProvider(name="doubleIntervals")
    public Object[][] getDoubleIntervalData() {
        return new Object[][] {
                // disjoint on same contig, A upstream of B by several bases
                { new ClosedSVInterval(0, 1, 100), new ClosedSVInterval(0, 200, 300), false, 0, true, 99 },
                // disjoint on same contig, A upstream of B by 1 base (0-length gap)
                { new ClosedSVInterval(0, 1, 100), new ClosedSVInterval(0, 101, 300), false, 0, true, 0 },
                // disjoint on same contig, A downstream of B
                { new ClosedSVInterval(0, 1000, 1100), new ClosedSVInterval(0, 200, 300), false, 0, false, -1 },
                // disjoint on different contigs, A upstream of B
                { new ClosedSVInterval(0, 1, 100), new ClosedSVInterval(1, 200, 300), false, 0, true, Integer.MAX_VALUE },
                // disjoint on different contigs, A downstream of B
                { new ClosedSVInterval(1, 1, 100), new ClosedSVInterval(0, 200, 300), false, 0, false, Integer.MAX_VALUE },
                // overlapping by several bases, A starts before B
                { new ClosedSVInterval(0, 1, 100), new ClosedSVInterval(0, 51, 200), true, 50, false, -1 },
                // overlapping by 1 base, A starts before B
                { new ClosedSVInterval(0, 1, 100), new ClosedSVInterval(0, 100, 200), true, 1, false, -1 },
                // overlapping by several bases, B starts before A
                { new ClosedSVInterval(0, 100, 200), new ClosedSVInterval(0, 50, 200), true, 101, false, -1 }
        };
    }


    @Test(dataProvider = "doubleIntervals")
    public void testCompareIntervals(
            final ClosedSVInterval A,
            final ClosedSVInterval B,
            final boolean expectedOverlapBool,
            final int expectedOverlapLen,
            final boolean expectedIsUpstreamBool,
            final int expectedGapLen)
    {
        Assert.assertEquals(A.overlaps(B), expectedOverlapBool);
        Assert.assertEquals(A.isDisjointFrom(B), !expectedOverlapBool);
        Assert.assertEquals(A.overlapLen(B), expectedOverlapLen);
        Assert.assertEquals(A.isUpstreamOf(B), expectedIsUpstreamBool);
        if (expectedIsUpstreamBool) {
            Assert.assertEquals(A.gapLen(B), expectedGapLen);
        }

    }
}
