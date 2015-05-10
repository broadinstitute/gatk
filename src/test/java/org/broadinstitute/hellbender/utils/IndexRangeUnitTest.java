package org.broadinstitute.hellbender.utils;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

/**
 * Unit tests for {@link IndexRange}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class IndexRangeUnitTest extends BaseTest {

    @Test(dataProvider = "correctFromToData")
    public void testCorrectConstruction(final int from, final int to) {
        final IndexRange range = new IndexRange(from,to);
        Assert.assertEquals(range.from, from);
        Assert.assertEquals(range.to, to);
    }

    @Test(dataProvider = "correctFromToData", dependsOnMethods = "testCorrectConstruction")
    public void testSize(final int from, final int to) {
        final IndexRange range = new IndexRange(from,to);
        Assert.assertEquals(range.size(),to - from);
    }

    @Test(dataProvider = "correctFromToData", dependsOnMethods = "testCorrectConstruction")
    public void testValidLength(final int from, final int to) {
        final IndexRange range = new IndexRange(from,to);
        Assert.assertTrue(range.isValidLength(to));
        Assert.assertFalse(range.isValidLength(to - 1));
        Assert.assertTrue(range.isValidLength(to + 1));
        Assert.assertFalse(range.isValidLength(-1));
        Assert.assertFalse(range.isValidLength(-10));
    }

    @Test(dataProvider = "correctFromToData", dependsOnMethods = "testCorrectConstruction")
    public void testForEach(final int from, final int to) {
        final IndexRange range = new IndexRange(from,to);
        final List<Integer> indexes = new ArrayList<>(to -from);
        range.forEach(i -> indexes.add(i));
        Assert.assertEquals(indexes.size(), to - from);
        for (int i = 0; i < indexes.size(); i++) {
            Assert.assertEquals(indexes.get(i).intValue(),from + i);
        }
    }

    @Test(dataProvider = "correctFromToData", dependsOnMethods = "testCorrectConstruction",
            expectedExceptions = IllegalArgumentException.class)
    public void testNullForEachConsumer(final int from, final int to) {
        final IndexRange range = new IndexRange(from,to);
        range.forEach(null);
    }

    @Test(dataProvider = "wrongFromToData",
            expectedExceptions = IllegalArgumentException.class)
    public void testWrongConstruction(final int from, final int to) {
        new IndexRange(from,to);
    }

    @Test(dependsOnMethods = "testCorrectConstruction")
    public void testEquals() {
        final IndexRange range01 = new IndexRange(0,1);
        final IndexRange range11 = new IndexRange(1,1);
        final IndexRange range00 = new IndexRange(0,0);
        final IndexRange range12 = new IndexRange(1,2);

        Assert.assertNotEquals(range01,range00);
        Assert.assertNotEquals(range01,range11);
        Assert.assertNotEquals(range01,range12);
        Assert.assertFalse(range01.equals(null));
        Assert.assertEquals(range01,new IndexRange(range01.from,range01.to));
        Assert.assertFalse(range01.equals(new Object()));
    }

    @DataProvider(name= "correctFromToData")
    public Object[][] correctFromToData() {
        return new Object[][] {
                { 0, 0},
                { 0, 11},
                { 0, 1},
                { 1, 14},
                {14, 14}
        };
    }

    @DataProvider(name= "wrongFromToData")
    public Object[][] wrongFromToData() {
        return new Object[][] {
                { -1, -1},
                { -1, -2},
                { 1, 0},
                { 15, 14},
                { -1, 0},
        };
    }
}
