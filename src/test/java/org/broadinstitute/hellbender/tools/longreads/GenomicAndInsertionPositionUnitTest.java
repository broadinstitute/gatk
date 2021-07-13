package org.broadinstitute.hellbender.tools.longreads;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.longreads.graph.GenomicAndInsertionPosition;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

/**
 * Unit tests for the {@link org.broadinstitute.hellbender.tools.longreads.graph.GenomicAndInsertionPosition} class.
 */
public class GenomicAndInsertionPositionUnitTest extends GATKBaseTest {

    private GenomicAndInsertionPosition makePos(final String contig, final int start, final int insertionOffset) {
        return new GenomicAndInsertionPosition(contig, start, insertionOffset);
    }

    @DataProvider
    private Object[][] provideForIsAdjacentTo() {
        return new Object[][] {
            { makePos("x", 1, 0), makePos("x", 1, 0), false },
            { makePos("x", 1, 1), makePos("x", 1, 0), true },
            { makePos("x", 1, 0), makePos("x", 2, 0), true },
            { makePos("x", 1, 0), makePos("y", 2, 0), true },
            { makePos("x", 1, 1), makePos("x", 2, 0), true },
            { makePos("x", 1, 1), makePos("x", 2, 2), false },
            { makePos("x", 1, 1), makePos("x", 2, 1), false },
        };
    }

    @Test(dataProvider = "provideForIsAdjacentTo")
    public void testIsAdjacentTo(final GenomicAndInsertionPosition pos1, final GenomicAndInsertionPosition pos2, final boolean expected) {
        Assert.assertEquals(pos1.isAdjacentTo(pos2), expected);
        Assert.assertEquals(pos2.isAdjacentTo(pos1), expected);
    }

}
