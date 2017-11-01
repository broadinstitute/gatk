package org.broadinstitute.hellbender.utils.collections;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.util.Arrays;

/**
 * Basic unit test for ExpandingArrayList
 */
public class ExpandingArrayListUnitTest extends GATKBaseTest {
    ExpandingArrayList<Integer> empty, initCap10, hasOne, hasTen;

    @BeforeMethod
    public void before() {
        empty = new ExpandingArrayList<>();

        initCap10 = new ExpandingArrayList<>(10);

        hasOne = new ExpandingArrayList<>();
        hasOne.add(1);

        hasTen = new ExpandingArrayList<>();
        hasTen.addAll(Arrays.asList(1, 2, 3, 4, 5, 6, 7, 8, 9, 10));
    }

    @Test
    public void testBasicSizes() {
        logger.warn("Executing testBasicSizes");

        Assert.assertEquals(0, empty.size());
        Assert.assertEquals(0, initCap10.size());
        Assert.assertEquals(1, hasOne.size());
        Assert.assertEquals(10, hasTen.size());
    }

    @Test
    public void testTenElements() {
        logger.warn("Executing testTenElements");

        for ( int i = 0; i < 10; i++ ) {
            Assert.assertEquals(i+1, (int)hasTen.get(i));
        }
    }

    @Test
    public void testSettingTenElements() {
        logger.warn("Executing testSettingTenElements");

        for ( int i = 0; i < 10; i++ ) {
            Assert.assertEquals(i+1, (int)hasTen.set(i, 2*i));
        }

        Assert.assertEquals(10, hasTen.size());
        for ( int i = 0; i < 10; i++ ) {
            Assert.assertEquals(2*i, (int)hasTen.get(i));
        }
    }

    @Test
    public void testAdd() {
        logger.warn("Executing testAdd");

        Assert.assertEquals(0, empty.size());
        empty.add(1);
        Assert.assertEquals(1, empty.size());
        Assert.assertEquals(1, (int)empty.get(0));
        empty.add(2);
        Assert.assertEquals(2, empty.size());
        Assert.assertEquals(2, (int)empty.get(1));
    }

    @Test
    public void testSet1() {
        logger.warn("Executing testSet1");

        Assert.assertEquals(0, empty.size());
        empty.set(0, 1);
        Assert.assertEquals(1, empty.size());
        Assert.assertEquals(1, (int)empty.get(0));

        empty.set(1, 2);
        Assert.assertEquals(2, empty.size());
        Assert.assertEquals(2, (int)empty.get(1));

        // doesn't expand
        empty.set(0, 3);
        Assert.assertEquals(2, empty.size());
        Assert.assertEquals(3, (int)empty.get(0));
    }

    @Test
    public void testSetExpanding() {
        logger.warn("Executing testSetExpanding");

        Assert.assertEquals(0, empty.size());
        empty.set(3, 1);
        Assert.assertEquals(4, empty.size());
        Assert.assertEquals(empty.get(0), null);
        Assert.assertEquals(empty.get(1), null);
        Assert.assertEquals(empty.get(2), null);
        Assert.assertEquals(1, (int)empty.get(3));
    }

    @Test
    public void testSetExpandingReset() {
        logger.warn("Executing testSetExpandingReset");

        Assert.assertEquals(0, empty.size());
        empty.set(3, 3);
        empty.set(2, 2);
        empty.set(1, 1);
        empty.set(0, 0);
        Assert.assertEquals(4, empty.size());
        for ( int i = 0; i < 4; i++ )
            Assert.assertEquals(i, (int)empty.get(i));
    }

    @Test
    public void testSetExpandingBig() {
        logger.warn("Executing testSetExpandingBig");

        Assert.assertEquals(0, empty.size());
        empty.set(1000, 1000);
        Assert.assertEquals(1001, empty.size());
        for ( int i = 0; i < 1000; i++ )
            Assert.assertEquals(empty.get(i), null);
        Assert.assertEquals(1000, (int)empty.get(1000));
    }

    @Test (expectedExceptions=IndexOutOfBoundsException.class )
    public void testSetBadGetNegative() {
        logger.warn("Executing testSetBadGetNegative");
        empty.get(-1);
    }

    @Test
    public void testSetBadGetPost() {
        logger.warn("Executing testSetBadGetPost");
        empty.set(1, 1);
        Assert.assertEquals(empty.get(0), null);
        Assert.assertEquals(1, (int)empty.get(1));
        Assert.assertEquals(empty.get(2), null);
    }
}
