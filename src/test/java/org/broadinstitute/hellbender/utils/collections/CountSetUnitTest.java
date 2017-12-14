package org.broadinstitute.hellbender.utils.collections;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Unit tests for {@link CountSet}
 */
public final class CountSetUnitTest extends GATKBaseTest {

    @Test(dataProvider="capacities")
    public void testSize(final int capacity) {
        final CountSet empty = new CountSet(capacity);
        Assert.assertEquals(empty.size(), 0);
        CountSet nonEmpty = new CountSet(capacity);
        for (int i = 0; i < capacity * 3; i++) {
            nonEmpty.add(i);
            Assert.assertEquals(nonEmpty.size(), i + 1);
        }
    }

    @Test
    public void testSingleValueAdd() {
        final int CAPACITY = 10;
        final CountSet subject = new CountSet(CAPACITY);
        final Set<Integer> reasuranceSet = new LinkedHashSet<>(CAPACITY);
        final int REPEATS = 1000;
        final Random rnd = new Random(13);
        for (int i = 0; i < REPEATS; i++) {
            int newInt = rnd.nextInt(500);
            boolean expectedResult = reasuranceSet.add(newInt);
            boolean result = subject.add(newInt);
            Assert.assertEquals(result, expectedResult);
            Assert.assertEquals(subject.size(), reasuranceSet.size());
        }
        for (final int j : reasuranceSet)
            Assert.assertTrue(subject.contains(j));
        for (int j = 0; j < 501; j++)
            Assert.assertEquals(subject.contains(j), reasuranceSet.contains(j));
    }

    @Test
    public void testToIntArray() {
        final CountSet subject = new CountSet(10);
        subject.addAll(1,4,7);
        final int[] intArray = subject.toIntArray();
        Assert.assertEquals(intArray.length, 3);
        Assert.assertEquals(intArray[0], 1);
        Assert.assertEquals(intArray[1], 4);
        Assert.assertEquals(intArray[2], 7);
    }

    @Test
    public void testCopyTo() {
        final CountSet subject = new CountSet(10);
        subject.addAll(1,4,7);
        final int[] intArray = new int[3];
        subject.copyTo(intArray);
        Assert.assertEquals(intArray[0], 1);
        Assert.assertEquals(intArray[1], 4);
        Assert.assertEquals(intArray[2], 7);
    }

    @Test
    public void testSetToSingleValue() {
        final CountSet subject = new CountSet(10);
        subject.setTo(-31);
        Assert.assertEquals(subject.size(), 1);
        Assert.assertEquals(subject.min(), -31);
        Assert.assertEquals(subject.max(), -31);
        Assert.assertTrue(subject.contains(-31));
        Assert.assertFalse(subject.contains(-21));
    }

    @Test
    void testSetToArrayOfValues() {
        final int CAPACITY = 10;
        final CountSet subject = new CountSet(CAPACITY);
        final int REPEATS = 1000;
        final Random rnd = new Random(13);
        final int[] values = new int[REPEATS];
        for (int i = 0; i < REPEATS; i++) {
            int newInt = rnd.nextInt(Integer.MAX_VALUE) * (rnd.nextBoolean() ? -1 : 1);
            values[i] = newInt;
        }
        subject.setTo(values);
        Arrays.sort(values);
        Assert.assertEquals(subject.size(), REPEATS);
        Assert.assertEquals(subject.min(), values[0]);
        Assert.assertEquals(subject.max(), values[REPEATS - 1]);
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
        subject.addAll(values);
        Arrays.sort(values);
        Assert.assertEquals(subject.min(), values[0]);
        Assert.assertEquals(subject.max(), values[REPEATS - 1]);
    }

    @Test
    public void testIncrease() {
        final int CAPACITY = 10;
        final CountSet subject = new CountSet(CAPACITY);
        final Set<Integer> reasuranceSet = new LinkedHashSet<>(CAPACITY);
        final int REPEATS = 1000;
        final Random rnd = new Random(13);
        final int[] values = new int[REPEATS];
        final Integer[] valueWrappers = new Integer[REPEATS];
        for (int i = 0; i < REPEATS; i++) {
            int newInt = rnd.nextInt(500);
            values[i] = newInt;
            valueWrappers[i] = newInt;
        }

        subject.incAll(3);

        for (final int j : reasuranceSet)
            Assert.assertTrue(subject.contains(j + 3));
        for (int j = 0; j < 501; j++)
            Assert.assertEquals(subject.contains(j + 3), reasuranceSet.contains(j));

    }

    @Test
    public void testArrayValueAdd() {
        final int CAPACITY = 10;
        final CountSet subject = new CountSet(CAPACITY);
        final Set<Integer> reasuranceSet = new LinkedHashSet<>(CAPACITY);
        final int REPEATS = 1000;
        final Random rnd = new Random(13);
        final int[] values = new int[REPEATS];
        final Integer[] valueWrappers = new Integer[REPEATS];
        for (int i = 0; i < REPEATS; i++) {
            int newInt = rnd.nextInt(500);
            values[i] = newInt;
            valueWrappers[i] = newInt;
        }

        boolean expectedResult = reasuranceSet.addAll(Arrays.asList(valueWrappers));
        boolean result = subject.addAll(values);
        Assert.assertEquals(result, expectedResult);
        Assert.assertEquals(subject.size(), reasuranceSet.size());

        for (final int j : reasuranceSet)
            Assert.assertTrue(subject.contains(j));
        for (int j = 0; j < 501; j++)
            Assert.assertEquals(subject.contains(j), reasuranceSet.contains(j));

    }

    @Test
    public void testAddRange() {
        final CountSet subject = new CountSet(10);
        subject.addRange(10,21);
        Assert.assertEquals(subject.size(), 12);
        for (int i = 10; i < 22; i++)
            Assert.assertTrue(subject.contains(i));
        for (int i = -1; i < 10; i++)
            Assert.assertFalse(subject.contains(i));
        for (int i = 22; i < 31; i++)
            Assert.assertFalse(subject.contains(i));
    }

    @Test
    public void tesIterator() {
        final CountSet subject = new CountSet(10);
        final int[] sortedArray = new int[]{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        subject.addAll(sortedArray);
        final Iterator<Integer> iterator = subject.iterator();
        for (int i = 0; i < sortedArray.length; i++) {
            Assert.assertTrue(iterator.hasNext());
            Assert.assertEquals(iterator.next().intValue(), sortedArray[i]);
        }
        Assert.assertFalse(iterator.hasNext());
    }

    @DataProvider(name="capacities")
    public Iterator<Object[]> capacities() {
        final int MIN = 0;
        final int MAX = 255;
        return new Iterator<Object[]>() {
                private int current = MIN;


            @Override
            public boolean hasNext() {
                return current < MAX;
            }

            @Override
            public Object[] next() {
                return new Object[] { Integer.valueOf(current++) };
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }
}
