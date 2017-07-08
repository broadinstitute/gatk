package org.broadinstitute.hellbender.tools.spark.utils;

import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Objects;

public final class HopscotchMapTest extends BaseTest {

    // a HopscotchMap is just like a HopscotchSet, which is separately tested, except for two things:
    //   the uniqueness criterion is on key rather than on the entire entry,
    //   retrieval is by key rather than by entry.
    // so we'll just test that new behavior here.

    private static final List<IntPair> testVals = Arrays.asList(new IntPair(1,2), new IntPair(3,4), new IntPair(5,6));
    private static final IntPair notInTestVals = new IntPair(7,8);

    @Test
    void addTest() {
        final HopscotchMap<Integer, Integer, IntPair> hopscotchMap = new HopscotchMap<>();
        Assert.assertTrue(hopscotchMap.add(new IntPair(1,2)));
        Assert.assertFalse(hopscotchMap.add(new IntPair(1,2)));
        Assert.assertFalse(hopscotchMap.add(new IntPair(1,3)));
        Assert.assertTrue(hopscotchMap.add(new IntPair(2,2)));
        Assert.assertEquals(hopscotchMap.size(), 2);
    }

    @Test
    void containsTest() {
        final HopscotchMap<Integer, Integer, IntPair> hopscotchMap = new HopscotchMap<>(testVals);
        Assert.assertEquals(hopscotchMap.size(), testVals.size());
        for ( final IntPair entry : testVals ) {
            Assert.assertTrue(hopscotchMap.contains(entry.getKey()));
        }
        Assert.assertFalse(hopscotchMap.contains(notInTestVals));
    }

    @Test
    void findTest() {
        final HopscotchMap<Integer, Integer, IntPair> hopscotchMap = new HopscotchMap<>(testVals);
        for ( final IntPair entry : testVals ) {
            Assert.assertEquals(hopscotchMap.find(entry.getKey()), entry);
        }
        Assert.assertNull(hopscotchMap.find(notInTestVals));
    }

    @Test
    void findEachTest() {
        final HopscotchMap<Integer, Integer, IntPair> hopscotchMap = new HopscotchMap<>(testVals);
        for ( final IntPair entry : testVals ) {
            Assert.assertEquals(SVUtils.iteratorSize(hopscotchMap.findEach(entry.getKey())), 1);
        }
        Assert.assertEquals(SVUtils.iteratorSize(hopscotchMap.findEach(notInTestVals)), 0);
    }

    @Test
    void removeTest() {
        final HopscotchMap<Integer, Integer, IntPair> hopscotchMap = new HopscotchMap<>(testVals);
        for ( final IntPair entry : testVals ) {
            Assert.assertTrue(hopscotchMap.remove(entry.getKey()));
            Assert.assertFalse(hopscotchMap.remove(entry.getKey()));
        }
        Assert.assertEquals(hopscotchMap.size(), 0);
    }

    @Test
    void removeEachTest() {
        final HopscotchMap<Integer, Integer, IntPair> hopscotchMap = new HopscotchMap<>(testVals);
        for ( final IntPair entry : testVals ) {
            Assert.assertTrue(hopscotchMap.removeEach(entry.getKey()));
            Assert.assertFalse(hopscotchMap.removeEach(entry.getKey()));
        }
        Assert.assertEquals(hopscotchMap.size(), 0);
    }

    static final class IntPair implements Map.Entry<Integer, Integer> {
        private final int key;
        private int value;

        IntPair( final int key, final int value ) {
            this.key = key;
            this.value = value;
        }

        @Override
        public Integer getKey() { return key; }

        @Override
        public Integer getValue() { return value; }

        @Override
        public Integer setValue( final Integer value ) {
            final int oldValue = this.value;
            this.value = value;
            return oldValue;
        }

        @Override
        public boolean equals( final Object obj ) {
            if ( !(obj instanceof IntPair) ) return false;
            final IntPair that = (IntPair)obj;
            return Objects.equals(this.key, that.key) && Objects.equals(this.value, that.value);
        }

        @Override
        public int hashCode() { return 47 * (47 * key + value); }
    }
}
