package org.broadinstitute.hellbender.tools.spark.utils;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public final class LongHopscotchSetTest extends BaseTest {
    private static final long[] testVals = {0, 1, 2, 8, 16, 42, 97, 100, 2894765};
    private static final long[] notAllTestVals = {0, 1, 2, 3, 7, 22, 61};
    private static final long notInTestVals = 6;
    private static final int RAND_SEED = 0xdeadf00;
    private static final int HHASH_NVALS = 100000;
    private static final int N_TRIALS = 10;

    private static long randomLong(Random rng) {
        return (((long) rng.nextInt()) | (((long) rng.nextInt()) << 31)) & ~Long.MIN_VALUE;
    }

    @Test
    void createFromCollectionTest() {
        final LongHopscotchSet hopscotchSet = new LongHopscotchSet(testVals);
        Assert.assertEquals(hopscotchSet.size(), testVals.length);
        Assert.assertTrue(hopscotchSet.containsAll(testVals));
    }

    @Test
    void getLegalSizeBelowTest() {
        Assert.assertEquals(SetSizeUtils.getLegalSizeBelow((long)(1000 * LongHopscotchSet.LOAD_FACTOR)), 719);
        Assert.assertEquals(SetSizeUtils.getLegalSizeBelow((long)(3e9 * LongHopscotchSet.LOAD_FACTOR)), 2147483629);
    }

    @Test
    void bytesPerBucketTest() {
        Assert.assertTrue(LongHopscotchSet.bytesPerEntry > 0);
    }

    @Test
    void addTest() {
        final LongHopscotchSet hopscotchSet = new LongHopscotchSet(testVals.length);
        for (final long val : testVals) {
            hopscotchSet.add(val);
        }
        Assert.assertEquals(hopscotchSet.size(), testVals.length);
        Assert.assertTrue(hopscotchSet.containsAll(testVals));
    }

    @Test
    void clearTest() {
        final LongHopscotchSet hopscotchSet = new LongHopscotchSet(testVals);
        hopscotchSet.clear();
        Assert.assertEquals(hopscotchSet.size(), 0);
        Assert.assertEquals(SVUtils.iteratorSize(hopscotchSet.iterator()), 0);
    }

    @Test
    void capacityTest() {
        for (final int size : new int[]{1000, 2000, 5000, 7000, 12000}) {
            final long capacity = new LongHopscotchSet(size).capacity();
            Assert.assertTrue(capacity >= size);
            Assert.assertTrue(capacity < 2 * size);
            final List<Integer> legalSizes = IntStream.of(SetSizeUtils.legalSizes).boxed().collect(Collectors.toList());
            Assert.assertTrue(legalSizes.contains(new Integer((int) capacity)));
        }
    }

    @Test
    void containsTest() {
        final LongHopscotchSet hopscotchSet = new LongHopscotchSet(testVals);
        Assert.assertTrue(hopscotchSet.containsAll(testVals));
        Assert.assertFalse(hopscotchSet.contains(notInTestVals));
        Assert.assertFalse(hopscotchSet.containsAll(notAllTestVals));
    }

    @Test
    void isEmptyTest() {
        final LongHopscotchSet hopscotchSet = new LongHopscotchSet();
        Assert.assertTrue(hopscotchSet.isEmpty());
        hopscotchSet.add(1);
        Assert.assertFalse(hopscotchSet.isEmpty());
    }

    @Test
    void iteratorTest() {
        final LongHopscotchSet hopscotchSet = new LongHopscotchSet(testVals);
        Assert.assertEquals(SVUtils.iteratorSize(hopscotchSet.iterator()), testVals.length);
        final Integer KEY_TO_REMOVE = 1;
        final LongIterator itr1 = hopscotchSet.iterator();
        while (itr1.hasNext()) {
            if (itr1.next() == KEY_TO_REMOVE)
                itr1.remove();
        }

        final List<Long> filteredVals = new LinkedList<>();
        for (long testVal : testVals) {
            if (KEY_TO_REMOVE == testVal) {
                filteredVals.add(new Long(testVal));
            }
        }
        final int onesCount = filteredVals.size();
        Assert.assertEquals(hopscotchSet.size(), testVals.length - onesCount);
        Assert.assertEquals(SVUtils.iteratorSize(hopscotchSet.iterator()), hopscotchSet.size());
        final LongIterator itr2 = hopscotchSet.iterator();
        while (itr2.hasNext()) {
            itr2.next();
            itr2.remove();
        }
        Assert.assertTrue(hopscotchSet.isEmpty());
        Assert.assertEquals(SVUtils.iteratorSize(hopscotchSet.iterator()), hopscotchSet.size());
    }

    @Test
    void removeTest() {
        final LongHopscotchSet hopscotchSet = new LongHopscotchSet(testVals);
        Assert.assertFalse(hopscotchSet.remove(notInTestVals));
        for (final long value : testVals) {
            Assert.assertTrue(hopscotchSet.remove(value));
        }
        Assert.assertTrue(hopscotchSet.isEmpty());
    }

    @Test
    void removeAllTest() {
        final long[] notInTestValsArray = {notInTestVals};
        final LongHopscotchSet hopscotchSet = new LongHopscotchSet(testVals);
        Assert.assertFalse(hopscotchSet.removeAll(notInTestValsArray));
        Assert.assertTrue(hopscotchSet.removeAll(testVals));
        hopscotchSet.addAll(testVals);
        Assert.assertTrue(hopscotchSet.removeAll(new LongHopscotchSet(testVals)));
        Assert.assertTrue(hopscotchSet.isEmpty());
    }

    @Test
    void sizeTest() {
        final LongHopscotchSet hopscotchSet = new LongHopscotchSet();
        Assert.assertEquals(hopscotchSet.size(), 0);
        hopscotchSet.addAll(testVals);
        Assert.assertEquals(hopscotchSet.size(), testVals.length);
        hopscotchSet.add(notInTestVals);
        Assert.assertEquals(hopscotchSet.size(), testVals.length + 1);
        hopscotchSet.remove(notInTestVals);
        Assert.assertEquals(hopscotchSet.size(), testVals.length);
        hopscotchSet.clear();
        Assert.assertEquals(hopscotchSet.size(), 0);
    }

    @Test
    void resizeTest() {
        final LongHopscotchSet hopscotchSet = new LongHopscotchSet(1);
        final int N_ENTRIES = 1000000;
        final LongHopscotchSet hopscotchSet_2 = new LongHopscotchSet(N_ENTRIES);
        for (int idx = 0; idx != N_ENTRIES; ++idx) {
            hopscotchSet.add(idx);
            hopscotchSet_2.add(idx);
        }
        Assert.assertEquals(hopscotchSet, hopscotchSet_2);
    }

    @Test
    void noDupsTest() {
        final LongHopscotchSet hopscotchSet = new LongHopscotchSet();
        Assert.assertTrue(hopscotchSet.add(1));
        Assert.assertFalse(hopscotchSet.add(1));
        Assert.assertEquals(hopscotchSet.size(), 1);
    }

    @Test
    void loadRandomLongsTest() {
        final Random rng = new Random(RAND_SEED);
        for (int trialNo = 0; trialNo != N_TRIALS; ++trialNo) {
            final HashSet<Long> hashSet = new HashSet<>();
            final LongHopscotchSet hopscotchSet = new LongHopscotchSet(HHASH_NVALS);
            final String trialMsg = "trialNo=" + trialNo;
            for (int valNo = 0; valNo != HHASH_NVALS; ++valNo) {
                final long randLong = randomLong(rng);
                hopscotchSet.add(randLong);
                hashSet.add(new Long(randLong));
            }
            Assert.assertEquals(hashSet.size(), hopscotchSet.size(), trialMsg);
            for (final Long val : hashSet) {
                Assert.assertTrue(hopscotchSet.contains(val), trialMsg + ", testVal=" + val);
            }
            final LongIterator itr = hopscotchSet.iterator();
            while (itr.hasNext()) {
                final long val = itr.next();
                Assert.assertTrue(hashSet.contains(val), trialMsg + ", testVal=" + val);
            }

            for (final Long val : hashSet) {
                Assert.assertTrue(hopscotchSet.remove(val));
            }
            hashSet.clear();

            Assert.assertEquals(hashSet.size(), hopscotchSet.size(), trialMsg);
            for (final Long val : hashSet) {
                Assert.assertTrue(hopscotchSet.contains(val), trialMsg + ", testVal=" + val);
            }
            final LongIterator itr2 = hopscotchSet.iterator();
            while (itr2.hasNext()) {
                final long val = itr2.next();
                Assert.assertTrue(hashSet.contains(val), trialMsg + ", testVal=" + val);
            }
        }
    }

    @Test
    void removeAllWithIteratorTest() {
        final Random rng = new Random(RAND_SEED);
        final LongHopscotchSet hopscotchSet = new LongHopscotchSet(HHASH_NVALS);
        for (int valNo = 0; valNo != HHASH_NVALS; ++valNo) {
            final long longVal = ((long) rng.nextInt()) & ~Long.MIN_VALUE;
            hopscotchSet.add(longVal);
        }
        final LongIterator itr = hopscotchSet.iterator();
        while (itr.hasNext()) {
            final long next = itr.next();
            itr.remove();
            Assert.assertFalse(hopscotchSet.contains(next));
        }
        Assert.assertEquals(hopscotchSet.size(), 0);
    }

    @Test
    void equalsAndHashcodeTest() {
        final LongHopscotchSet hopscotchSet1 = new LongHopscotchSet(testVals.length);
        final LongHopscotchSet hopscotchSet2 = new LongHopscotchSet(testVals.length);
        final LongHopscotchSet hopscotchSet3 = new LongHopscotchSet(notAllTestVals.length);
        hopscotchSet1.addAll(testVals);
        hopscotchSet2.addAll(testVals);
        hopscotchSet3.addAll(notAllTestVals);
        Assert.assertEquals(hopscotchSet1.hashCode(), hopscotchSet2.hashCode());
        Assert.assertEquals(hopscotchSet1, hopscotchSet2);
        Assert.assertNotEquals(hopscotchSet1.hashCode(), hopscotchSet3.hashCode());
        Assert.assertNotEquals(hopscotchSet1, hopscotchSet3);
    }

    @Test
    void serializationTest() {
        final Random rng = new Random(RAND_SEED);
        final LongHopscotchSet hopscotchSet = new LongHopscotchSet(HHASH_NVALS);
        for (int valNo = 0; valNo != HHASH_NVALS; ++valNo) {
            final long randLong = randomLong(rng);
            hopscotchSet.add(randLong);
        }

        final ByteArrayOutputStream bos = new ByteArrayOutputStream();
        final Output out = new Output(bos);
        final Kryo kryo = new Kryo();
        kryo.writeObject(out, hopscotchSet);
        out.flush();

        final ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());
        final Input in = new Input(bis);
        final LongHopscotchSet hopscotchSet2 = kryo.readObject(in, LongHopscotchSet.class);

        Assert.assertEquals(hopscotchSet.size(), hopscotchSet2.size());
        final LongIterator itr1 = hopscotchSet.iterator();
        while (itr1.hasNext()) {
            final long val = itr1.next();
            Assert.assertTrue(hopscotchSet2.contains(val));
        }
        final LongIterator itr2 = hopscotchSet2.iterator();
        while (itr2.hasNext()) {
            final long val = itr2.next();
            Assert.assertTrue(hopscotchSet.contains(val));
        }
    }
}
