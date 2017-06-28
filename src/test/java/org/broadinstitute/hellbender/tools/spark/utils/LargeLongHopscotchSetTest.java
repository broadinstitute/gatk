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
import java.util.Collection;
import java.util.HashSet;
import java.util.Random;

public final class LargeLongHopscotchSetTest extends BaseTest {

    private static final long[] testVals = {0, 1, 2, 8, 16, 42, 97, 100, 2894765};
    private static final long[] notAllTestVals = {0, 1, 2, 3, 7, 22, 61};
    private static final long notInTestVals = 6;
    private static final int RAND_SEED = 0xdeadf00;
    private static final int HHASH_NVALS = 100000;

    private static long randomLong(Random rng) {
        return (((long) rng.nextInt()) | (((long) rng.nextInt()) << 31)) & Long.MAX_VALUE;
    }

    @Test
    void getSetsTest() {
        final LargeLongHopscotchSet hopscotchSet = new LargeLongHopscotchSet(1);
        hopscotchSet.add(1);
        Assert.assertEquals(hopscotchSet.getSets().size(), 1);
        final LargeLongHopscotchSet hopscotchSet2 = new LargeLongHopscotchSet(64000);
        Assert.assertTrue(hopscotchSet2.getSets().size() > 1);
    }

    @Test
    void createFromCollectionTest() {
        final LargeLongHopscotchSet hopscotchSet = new LargeLongHopscotchSet(testVals.length);
        hopscotchSet.addAll(testVals);
        Assert.assertEquals(hopscotchSet.size(), testVals.length);
        Assert.assertTrue(hopscotchSet.containsAll(testVals));
    }

    @Test
    void addTest() {
        final LargeLongHopscotchSet hopscotchSet = new LargeLongHopscotchSet(testVals.length);
        for (final long val : testVals) {
            hopscotchSet.add(val);
        }
        Assert.assertEquals(hopscotchSet.size(), testVals.length);
        Assert.assertTrue(hopscotchSet.containsAll(testVals));
    }

    @Test
    void capacityTest() {
        for (final int size : new int[]{1000, 2000, 5000, 7000, 12000}) {
            final long capacity = new LargeLongHopscotchSet(size).capacity();
            Assert.assertTrue(capacity >= size);
        }
    }

    @Test
    void containsTest() {
        final LargeLongHopscotchSet hopscotchSet = new LargeLongHopscotchSet(testVals.length);
        hopscotchSet.addAll(testVals);
        Assert.assertTrue(hopscotchSet.containsAll(testVals));
        Assert.assertFalse(hopscotchSet.contains(notInTestVals));
        Assert.assertFalse(hopscotchSet.containsAll(notAllTestVals));
    }

    @Test
    void isEmptyTest() {
        final LargeLongHopscotchSet hopscotchSet = new LargeLongHopscotchSet(1);
        Assert.assertTrue(hopscotchSet.isEmpty());
        hopscotchSet.add(1);
        Assert.assertFalse(hopscotchSet.isEmpty());
    }

    @Test
    void iteratorTest() {
        final LargeLongHopscotchSet hopscotchSet = new LargeLongHopscotchSet(testVals.length);
        hopscotchSet.addAll(testVals);
        Assert.assertEquals(SVUtils.iteratorSize(hopscotchSet.iterator()), testVals.length);
    }

    @Test
    void removeTest() {
        final LargeLongHopscotchSet hopscotchSet = new LargeLongHopscotchSet(testVals.length);
        hopscotchSet.addAll(testVals);
        Assert.assertFalse(hopscotchSet.remove(notInTestVals));
        for (final long value : testVals) {
            Assert.assertTrue(hopscotchSet.remove(value));
        }
        Assert.assertTrue(hopscotchSet.isEmpty());
    }

    @Test
    void removeAllTest() {
        final long[] notInTestValsArray = {notInTestVals};
        final LargeLongHopscotchSet hopscotchSet = new LargeLongHopscotchSet(testVals.length);
        hopscotchSet.addAll(testVals);
        Assert.assertFalse(hopscotchSet.removeAll(notInTestValsArray));
        Assert.assertTrue(hopscotchSet.containsAll(testVals));
        Assert.assertTrue(hopscotchSet.removeAll(testVals));
        Assert.assertTrue(hopscotchSet.isEmpty());
    }

    @Test
    void sizeTest() {
        final LargeLongHopscotchSet hopscotchSet = new LargeLongHopscotchSet(testVals.length);
        Assert.assertEquals(hopscotchSet.size(), 0);
        hopscotchSet.addAll(testVals);
        Assert.assertEquals(hopscotchSet.size(), testVals.length);
        hopscotchSet.add(notInTestVals);
        Assert.assertEquals(hopscotchSet.size(), testVals.length + 1);
        hopscotchSet.remove(notInTestVals);
        Assert.assertEquals(hopscotchSet.size(), testVals.length);
    }

    @Test
    void resizeTest() {
        final Random rng = new Random(RAND_SEED);
        final LargeLongHopscotchSet hopscotchSet = new LargeLongHopscotchSet(HHASH_NVALS / 10);
        final Collection<Long> hashSet = new HashSet<>();
        for (int i = 0; i < HHASH_NVALS; i++) {
            final long val = randomLong(rng);
            hopscotchSet.add(val);
            hashSet.add(val);
        }
        Assert.assertTrue(hopscotchSet.size() == hashSet.size());
        LongIterator itrL = hopscotchSet.iterator();
        while (itrL.hasNext()) {
            Assert.assertTrue(hashSet.contains(new Long(itrL.next())));
        }
        for (Long val : hashSet) {
            Assert.assertTrue(hopscotchSet.contains(val.longValue()));
        }
    }

    @Test
    void noDupsTest() {
        final LargeLongHopscotchSet hopscotchSet = new LargeLongHopscotchSet(1);
        Assert.assertTrue(hopscotchSet.add(1));
        Assert.assertFalse(hopscotchSet.add(1));
        Assert.assertEquals(hopscotchSet.size(), 1);
    }

    @Test
    void loadRandomLongsTest() {
        final Random rng = new Random(RAND_SEED);
        final HashSet<Long> hashSet = new HashSet<>();
        final LargeLongHopscotchSet hopscotchSet = new LargeLongHopscotchSet(HHASH_NVALS);
        for (int valNo = 0; valNo != HHASH_NVALS; ++valNo) {
            final long randLong = randomLong(rng);
            hopscotchSet.add(randLong);
            hashSet.add(new Long(randLong));
        }
        Assert.assertEquals(hashSet.size(), hopscotchSet.size());
        for (final Long val : hashSet) {
            Assert.assertTrue(hopscotchSet.contains(val), "testVal=" + val);
        }
        final LongIterator itr = hopscotchSet.iterator();
        while (itr.hasNext()) {
            final long val = itr.next();
            Assert.assertTrue(hashSet.contains(val), "testVal=" + val);
        }

        for (final Long val : hashSet) {
            Assert.assertTrue(hopscotchSet.remove(val));
        }
        hashSet.clear();

        Assert.assertEquals(hashSet.size(), hopscotchSet.size());
        for (final Long val : hashSet) {
            Assert.assertTrue(hopscotchSet.contains(val), "testVal=" + val);
        }
        final LongIterator itr2 = hopscotchSet.iterator();
        while (itr2.hasNext()) {
            final long val = itr2.next();
            Assert.assertTrue(hashSet.contains(val), "testVal=" + val);
        }
    }

    @Test
    void equalsAndHashcodeTest() {
        final LargeLongHopscotchSet hopscotchSet1 = new LargeLongHopscotchSet(testVals.length);
        final LargeLongHopscotchSet hopscotchSet2 = new LargeLongHopscotchSet(testVals.length);
        final LargeLongHopscotchSet hopscotchSet3 = new LargeLongHopscotchSet(notAllTestVals.length);
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
        final LargeLongHopscotchSet hopscotchSet = new LargeLongHopscotchSet(HHASH_NVALS);
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
        final LargeLongHopscotchSet hopscotchSet2 = kryo.readObject(in, LargeLongHopscotchSet.class);

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
