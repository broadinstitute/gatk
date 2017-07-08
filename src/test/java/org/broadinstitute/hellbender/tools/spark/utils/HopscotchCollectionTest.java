package org.broadinstitute.hellbender.tools.spark.utils;

import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public final class HopscotchCollectionTest extends BaseTest {
    private static final List<Integer> testVals = Arrays.asList(1, 4, 3, 3, 2, 5, 1, 2, 2, 4, -1, -2894765);
    private static final Integer notInTestVals = 6;

    @Test
    void legalCapacitiesTest() {
        final int[] caps = SetSizeUtils.legalSizes;
        final int nCaps = caps.length;
        // test that they're spaced properly -- each is supposed to be about sqrt(2) bigger than the previous one
        for ( int idx = 1; idx < nCaps; ++idx ) {
            final double err = Math.abs(1. - Math.sqrt(2.)*caps[idx-1]/caps[idx]);
            Assert.assertTrue(err < .015, "testing capacity "+caps[idx]+" at index "+idx);
        }
        // test that they're all primes
        for ( int idx = 0; idx < nCaps; ++idx ) {
            Assert.assertTrue(isPrime(caps[idx]), "testing capacity "+caps[idx]+" at index "+idx);
        }
    }

    private static boolean isPrime( final long iii ) {
        if ( iii % 2 == 0 || iii % 3 == 0 ) return false;
        long iFact = 5;
        while ( iFact*iFact <= iii ) {
            if ( iii % iFact == 0 || iii % (iFact+2) == 0 ) return false;
            iFact += 6;
        }
        return true;
    }

    @Test
    void createFromCollectionTest() {
        final HopscotchCollection<Integer> hopscotchCollection = new HopscotchCollection<>(testVals);
        Assert.assertEquals(hopscotchCollection.size(), testVals.size());
        Assert.assertTrue(hopscotchCollection.containsAll(testVals));
    }

    @Test
    void addTest() {
        final HopscotchCollection<Integer> hopscotchCollection = new HopscotchCollection<>(testVals.size());
        testVals.stream().forEach(hopscotchCollection::add);
        Assert.assertEquals(hopscotchCollection.size(), testVals.size());
        Assert.assertTrue(hopscotchCollection.containsAll(testVals));
    }

    @Test
    void clearTest() {
        final HopscotchCollection<Integer> hopscotchCollection = new HopscotchCollection<>(testVals);
        hopscotchCollection.clear();
        Assert.assertEquals(hopscotchCollection.size(), 0);
        Assert.assertEquals(SVUtils.iteratorSize(hopscotchCollection.iterator()), 0);
    }

    @Test
    void capacityTest() {
        for ( final int size : new int[]{1000, 2000, 5000, 7000, 12000}) {
            final int capacity = new HopscotchCollection<>(size).capacity();
            Assert.assertTrue(capacity >= size);
            Assert.assertTrue(capacity < 2*size);
            final List<Integer> legalSizes = IntStream.of(SetSizeUtils.legalSizes).boxed().collect(Collectors.toList());
            Assert.assertTrue(legalSizes.contains(capacity));
        }
    }

    @Test
    void containsTest() {
        final HopscotchCollection<Integer> hopscotchCollection = new HopscotchCollection<>(testVals);
        Assert.assertTrue(hopscotchCollection.containsAll(testVals));
        Assert.assertFalse(hopscotchCollection.contains(notInTestVals));
    }

    @Test
    void findTest() {
        final HopscotchCollection<Integer> hopscotchCollection = new HopscotchCollection<>(testVals);
        for ( final Integer value : testVals ) {
            Assert.assertEquals(hopscotchCollection.find(value), value);
        }
        Assert.assertNull(hopscotchCollection.find(notInTestVals));
    }

    @Test
    void findEachTest() {
        final HopscotchCollection<Integer> hopscotchCollection = new HopscotchCollection<>(testVals);
        for ( final Integer value : testVals ) {
            Assert.assertEquals(SVUtils.iteratorSize(hopscotchCollection.findEach(value)),
                    testVals.stream().filter(i -> i.equals(value)).count());
        }
        Assert.assertEquals(SVUtils.iteratorSize(hopscotchCollection.findEach(notInTestVals)), 0);
    }

    @Test
    void isEmptyTest() {
        final HopscotchCollection<Integer> hopscotchCollection = new HopscotchCollection<>();
        Assert.assertTrue(hopscotchCollection.isEmpty());
        hopscotchCollection.add(1);
        Assert.assertFalse(hopscotchCollection.isEmpty());
    }

    @Test
    void iteratorTest() {
        final HopscotchCollection<Integer> hopscotchCollection = new HopscotchCollection<>(testVals);
        Assert.assertEquals(SVUtils.iteratorSize(hopscotchCollection.iterator()), testVals.size());
        final Integer KEY_TO_REMOVE = 1;
        final Iterator<Integer> itr1 = hopscotchCollection.iterator();
        while ( itr1.hasNext() ) {
            if ( itr1.next().equals(KEY_TO_REMOVE) )
                itr1.remove();
        }
        final int onesCount = (int)testVals.stream().filter(i -> i.equals(KEY_TO_REMOVE)).count();
        Assert.assertEquals(hopscotchCollection.size(), testVals.size()-onesCount);
        Assert.assertEquals(SVUtils.iteratorSize(hopscotchCollection.iterator()), hopscotchCollection.size());
        final Iterator<Integer> itr2 = hopscotchCollection.iterator();
        while ( itr2.hasNext() ) {
            itr2.next();
            itr2.remove();
        }
        Assert.assertTrue(hopscotchCollection.isEmpty());
        Assert.assertEquals(SVUtils.iteratorSize(hopscotchCollection.iterator()), hopscotchCollection.size());
    }

    @Test
    void removeTest() {
        final HopscotchCollection<Integer> hopscotchCollection = new HopscotchCollection<>(testVals);
        Assert.assertFalse(hopscotchCollection.remove(notInTestVals));
        for ( final Integer value : testVals ) {
            Assert.assertTrue(hopscotchCollection.remove(value));
        }
        Assert.assertTrue(hopscotchCollection.isEmpty());
    }

    @Test
    void removeEachTest() {
        final HopscotchCollection<Integer> hopscotchCollection = new HopscotchCollection<>(testVals);
        Assert.assertFalse(hopscotchCollection.removeEach(notInTestVals));
        for ( final Integer value : new HashSet<Integer>(testVals) ) {
            Assert.assertTrue(hopscotchCollection.removeEach(value));
        }
        Assert.assertTrue(hopscotchCollection.isEmpty());
    }

    @Test
    void removeAllTest() {
        final HopscotchCollection<Integer> hopscotchCollection = new HopscotchCollection<>(testVals);
        Assert.assertFalse(hopscotchCollection.removeAll(Arrays.asList(notInTestVals)));
        Assert.assertTrue(hopscotchCollection.removeAll(Arrays.asList(testVals.get(0))));
        Assert.assertTrue(hopscotchCollection.removeAll(new HashSet<Integer>(testVals)));
        Assert.assertTrue(hopscotchCollection.isEmpty());
    }

    @Test
    void sizeTest() {
        final HopscotchCollection<Integer> hopscotchCollection = new HopscotchCollection<>();
        Assert.assertEquals(hopscotchCollection.size(), 0);
        hopscotchCollection.addAll(testVals);
        Assert.assertEquals(hopscotchCollection.size(), testVals.size());
        hopscotchCollection.add(notInTestVals);
        Assert.assertEquals(hopscotchCollection.size(), testVals.size()+1);
        hopscotchCollection.removeEach(notInTestVals);
        Assert.assertEquals(hopscotchCollection.size(), testVals.size());
        hopscotchCollection.clear();
        Assert.assertEquals(hopscotchCollection.size(), 0);
    }

    @Test
    void resizeTest() {
        final HopscotchCollection<Integer> hopscotchCollection = new HopscotchCollection<>(1);
        final int N_ENTRIES = 1000000;
        for ( int idx = 0; idx != N_ENTRIES; ++idx ) {
            hopscotchCollection.add(idx);
        }
        Assert.assertEquals(hopscotchCollection.size(), N_ENTRIES);
        Assert.assertEquals(hopscotchCollection.stream().mapToInt(i->i).min().orElse(-1), 0);
        Assert.assertEquals(hopscotchCollection.stream().mapToInt(i->i).max().orElse(-1), N_ENTRIES-1);
    }
}
