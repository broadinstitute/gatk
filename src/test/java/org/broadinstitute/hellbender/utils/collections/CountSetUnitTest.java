package org.broadinstitute.hellbender.utils.collections;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.NoSuchElementException;
import java.util.Random;

/**
 * Unit tests for {@link CountSet}
 */
public final class CountSetUnitTest extends BaseTest {

    @Test(expectedExceptions = NoSuchElementException.class)
    public void testMinNoElement() throws Exception {
        new CountSet(1).min();
    }

    @Test(expectedExceptions = NoSuchElementException.class)
    public void testMaxNoElement() throws Exception {
        new CountSet(1).max();
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInvalidCtorArg() throws Exception {
        new CountSet(-1);
    }

    @Test
    public void testSize() {
        final int capacityLow = 0;
        final int capacityHigh = 10;
        for (int capacity = capacityLow; capacity < capacityHigh; capacity++) {
            final CountSet empty = new CountSet(capacity);
            Assert.assertEquals(empty.size(), 0);
            Assert.assertTrue(empty.isEmpty());
            CountSet nonEmpty = new CountSet(capacity);
            for (int i = 0; i < capacity * 3; i++) {
                nonEmpty.add(i);
                Assert.assertEquals(nonEmpty.size(), i + 1);
                Assert.assertFalse(nonEmpty.isEmpty());
            }
        }
    }

    @Test
    public void testToString() {
        final int capacityLow = 0;
        final int capacityHigh = 10;
        for (int capacity = capacityLow; capacity < capacityHigh; capacity++) {
            CountSet set = new CountSet(capacity);
            Assert.assertNotNull(set.toString());
            for (int i = 0; i < capacity * 3; i++) {
                set.add(i);
                Assert.assertEquals(set.size(), i + 1);
            }
            Assert.assertNotNull(set.toString());
        }
    }


    @Test
    public void testMinMax() {
        final int CAPACITY = 10;
        final CountSet subject = new CountSet(CAPACITY);
        final int REPEATS = 1000;
        final Random rnd = new Random(13);
        final int[] values = new int[REPEATS];
        for (int i = 0; i < REPEATS; i++) {
            int newInt = rnd.nextInt(Integer.MAX_VALUE) * (rnd.nextBoolean() ? -1 : 1);
            values[i] = newInt;
        }
        for (int i = 0; i < values.length; i++) {
            subject.add(values[i]);
        }
        Arrays.sort(values);
        Assert.assertEquals(subject.min(), values[0]);
        Assert.assertEquals(subject.max(), values[REPEATS - 1]);
    }
}
