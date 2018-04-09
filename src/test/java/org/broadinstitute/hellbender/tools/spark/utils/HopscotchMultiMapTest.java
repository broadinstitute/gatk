package org.broadinstitute.hellbender.tools.spark.utils;

import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchMapTest.IntPair;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;

public final class HopscotchMultiMapTest extends GATKBaseTest {

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

    static class CollidingHashMapKey {

        Integer value;

        public CollidingHashMapKey(final Integer value) {
            this.value = value;
        }

        @Override
        public int hashCode() {
            return 10;
        }

        public Integer getValue() {
            return value;
        }
    }

    static final class CollidingHashMapPair implements Map.Entry<CollidingHashMapKey, Integer> {
        private final CollidingHashMapKey key;
        private int value;

        CollidingHashMapPair(final int key, final int value) {
            this.key = new CollidingHashMapKey(key);
            this.value = value;
        }

        @Override
        public CollidingHashMapKey getKey() {
            return key;
        }

        @Override
        public Integer getValue() {
            return value;
        }

        @Override
        public Integer setValue(final Integer value) {
            final int oldValue = this.value;
            this.value = value;
            return oldValue;
        }
    }

    @Test
    public void testGroupingIterator() {
        HopscotchUniqueMultiMap<CollidingHashMapKey, Integer, CollidingHashMapPair> qNamesMultiMap = new HopscotchUniqueMultiMap<>(8);
        qNamesMultiMap.add(new CollidingHashMapPair(1, 1));
        qNamesMultiMap.add(new CollidingHashMapPair(2, 1));
        qNamesMultiMap.add(new CollidingHashMapPair(3, 1));
        qNamesMultiMap.add(new CollidingHashMapPair(4, 2));
        qNamesMultiMap.add(new CollidingHashMapPair(5, 2));
        qNamesMultiMap.add(new CollidingHashMapPair(2, 3));
        qNamesMultiMap.add(new CollidingHashMapPair(3, 3));
        qNamesMultiMap.add(new CollidingHashMapPair(6, 3));

        final Iterator<CollidingHashMapPair> iterator = qNamesMultiMap.getGroupingIterator();
        final Set<Integer> keysSeen = new HashSet<>(5);
        Integer prevQname = null;
        while (iterator.hasNext()) {
            CollidingHashMapPair next = iterator.next();
            if (prevQname != null) {
                Assert.assertTrue(! keysSeen.contains(next.getKey().value) || next.getKey().getValue().equals(prevQname));

            }
            prevQname = next.getKey().getValue();
            keysSeen.add(next.getKey().getValue());
        }
    }

}
