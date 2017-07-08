package org.broadinstitute.hellbender.tools.spark.utils;

import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchMapTest.IntPair;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public final class HopscotchMultiMapTest extends BaseTest {

    // a HopscotchMultiMap is just like a HopscotchCollection, which is separately tested,
    //    except that retrieval is by key rather than by entry.
    // so we'll just test that new behavior here.

    private static final List<IntPair> testVals = Arrays.asList(
            new IntPair(1,2), new IntPair(1,2), // repeated entries are legal in a multimap
            new IntPair(1,3), // as are repeated keys
            new IntPair(3,4), new IntPair(5,6));
    private static final IntPair notInTestVals = new IntPair(7,8);

    @Test
    void addTest() {
        final HopscotchMultiMap<Integer, Integer, IntPair> hopscotchMultiMap = new HopscotchMultiMap<>();
        for ( final IntPair entry : testVals ) {
            Assert.assertTrue(hopscotchMultiMap.add(entry));
        }
        Assert.assertEquals(hopscotchMultiMap.size(), testVals.size());
    }

    @Test
    void containsTest() {
        final HopscotchMultiMap<Integer, Integer, IntPair> hopscotchMultiMap = new HopscotchMultiMap<>(testVals);
        Assert.assertEquals(hopscotchMultiMap.size(), testVals.size());
        for ( final IntPair entry : testVals ) {
            Assert.assertTrue(hopscotchMultiMap.contains(entry.getKey()));
        }
        Assert.assertFalse(hopscotchMultiMap.contains(notInTestVals));
    }

    @Test
    void findTest() {
        final HopscotchMultiMap<Integer, Integer, IntPair> hopscotchMultiMap = new HopscotchMultiMap<>(testVals);
        for ( final IntPair entry : testVals ) {
            Assert.assertEquals(hopscotchMultiMap.find(entry.getKey()).getKey(), entry.getKey());
        }
        Assert.assertNull(hopscotchMultiMap.find(notInTestVals));
    }

    @Test
    void findEachTest() {
        final HopscotchMultiMap<Integer, Integer, IntPair> hopscotchMultiMap = new HopscotchMultiMap<>(testVals);
        for ( final Integer key : testVals.stream().map(IntPair::getKey).collect(Collectors.toSet()) ) {
            Assert.assertEquals(SVUtils.iteratorSize(hopscotchMultiMap.findEach(key)),
                    testVals.stream().filter(intPair -> intPair.getKey().equals(key)).count());
        }
        Assert.assertEquals(SVUtils.iteratorSize(hopscotchMultiMap.findEach(notInTestVals)), 0);
    }

    @Test
    void removeTest() {
        final HopscotchMultiMap<Integer, Integer, IntPair> hopscotchMultiMap = new HopscotchMultiMap<>(testVals);
        for ( final IntPair entry : testVals ) {
            Assert.assertTrue(hopscotchMultiMap.remove(entry.getKey()));
        }
        Assert.assertEquals(hopscotchMultiMap.size(), 0);
    }

    @Test
    void removeEachTest() {
        final HopscotchMultiMap<Integer, Integer, IntPair> hopscotchMultiMap = new HopscotchMultiMap<>(testVals);
        for ( final Integer key : testVals.stream().map(IntPair::getKey).collect(Collectors.toSet()) ) {
            Assert.assertTrue(hopscotchMultiMap.removeEach(key));
            Assert.assertFalse(hopscotchMultiMap.removeEach(key));
        }
        Assert.assertEquals(hopscotchMultiMap.size(), 0);
    }

    // found a bug in removeEach when multiple keys hash to the same bucket
    // (non-equivalent keys were getting deleted)
    // this test demonstrates the fix
    @Test
    void removeEachBugTest() {
        final HopscotchMultiMap<Integer, Integer, IntPair> hopscotchMultiMap = new HopscotchMultiMap<>();
        final int capacity = hopscotchMultiMap.capacity();
        hopscotchMultiMap.add(new IntPair(1, 1));
        hopscotchMultiMap.add(new IntPair(capacity+1, 1));
        hopscotchMultiMap.add(new IntPair(1, 2));
        hopscotchMultiMap.add(new IntPair(capacity+1, 2));
        hopscotchMultiMap.add(new IntPair(1, 3));
        Assert.assertTrue(hopscotchMultiMap.removeEach(1));
        Assert.assertEquals(hopscotchMultiMap.size(),2);
    }
}
