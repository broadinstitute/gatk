package org.broadinstitute.hellbender.tools.spark.utils;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Random;

public final class HopscotchSetTest extends GATKBaseTest {
    private static final int RAND_SEED = 0xdeadf00;
    private static final int HHASH_NVALS = 10000;
    private static final int N_TRIALS = 100;

    @Test
    void noDupsTest() {
        final HopscotchSet<Integer> hopscotchSet = new HopscotchSet<>();
        Assert.assertTrue(hopscotchSet.add(1));
        Assert.assertFalse(hopscotchSet.add(1));
        Assert.assertEquals(hopscotchSet.size(), 1);
    }

    @Test
    void findTest() {
        final HopscotchSet<Integer> hopscotchSet = new HopscotchSet<>();
        Assert.assertTrue(hopscotchSet.add(1));
        Assert.assertEquals(hopscotchSet.find(1), Integer.valueOf(1));
        Assert.assertNull(hopscotchSet.find(2));
    }

    @Test
    void loadRandomIntsTest() {
        final Random rng = new Random(RAND_SEED);
        for ( int trialNo = 0; trialNo != N_TRIALS; ++trialNo ) {
            final HashSet<Integer> hashSet = new HashSet<>();
            final HopscotchSet<Integer> hopscotchSet = new HopscotchSet<>(HHASH_NVALS);
            final String trialMsg = "trialNo="+trialNo;
            for ( int valNo = 0; valNo != HHASH_NVALS; ++valNo ) {
                final Integer randVal = rng.nextInt();
                hopscotchSet.add(randVal);
                hashSet.add(randVal);
            }
            Assert.assertEquals(hashSet.size(), hopscotchSet.size(), trialMsg);
            for ( final Integer val : hashSet ) {
                Assert.assertTrue(hopscotchSet.contains(val), trialMsg+", testVal="+val);
            }
            for ( final Integer val : hopscotchSet ) {
                Assert.assertTrue(hashSet.contains(val), trialMsg+", testVal="+val);
            }

            for ( final Integer val : hashSet ) {
                Assert.assertTrue(hopscotchSet.remove(val));
            }
            hashSet.clear();

            Assert.assertEquals(hashSet.size(), hopscotchSet.size(), trialMsg);
            for ( final Integer val : hashSet ) {
                Assert.assertTrue(hopscotchSet.contains(val), trialMsg+", testVal="+val);
            }
            for ( final Integer val : hopscotchSet ) {
                Assert.assertTrue(hashSet.contains(val), trialMsg+", testVal="+val);
            }
        }
    }

    @Test
    void removeAllWithIteratorTest() {
        final Random rng = new Random(RAND_SEED);
        final HopscotchSet<Integer> hopscotchSet = new HopscotchSet<>(HHASH_NVALS);
        for ( int valNo = 0; valNo != HHASH_NVALS; ++valNo ) {
            hopscotchSet.add(rng.nextInt());
        }
        final Iterator<Integer> itr = hopscotchSet.iterator();
        while ( itr.hasNext() ) {
            final Integer next = itr.next();
            itr.remove();
            Assert.assertFalse(hopscotchSet.contains(next));
        }
        Assert.assertEquals(hopscotchSet.size(), 0);
    }

    @Test
    void serializationTest() {
        final Random rng = new Random(RAND_SEED);
        final HopscotchSet<Integer> hopscotchSet = new HopscotchSet<>(HHASH_NVALS);
        for ( int valNo = 0; valNo != HHASH_NVALS; ++valNo ) {
            hopscotchSet.add(rng.nextInt());
        }

        final ByteArrayOutputStream bos = new ByteArrayOutputStream();
        final Output out = new Output(bos);
        final Kryo kryo = new Kryo();
        kryo.writeClassAndObject(out, hopscotchSet);
        out.flush();

        final ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());
        final Input in = new Input(bis);
        @SuppressWarnings("unchecked")
        final HopscotchSet<Integer> hopscotchSet2 = (HopscotchSet<Integer>)kryo.readClassAndObject(in);

        Assert.assertEquals(hopscotchSet.size(), hopscotchSet2.size());
        for ( final Integer val : hopscotchSet ) {
            Assert.assertTrue(hopscotchSet2.contains(val));
        }
        for ( final Integer val : hopscotchSet2 ) {
            Assert.assertTrue(hopscotchSet.contains(val));
        }
    }
}
