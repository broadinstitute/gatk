package org.broadinstitute.hellbender.tools.spark.utils;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.util.Log;
import org.broadinstitute.hellbender.utils.LoggingUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.util.HashSet;
import java.util.Random;

public final class LongBloomFilterTest extends BaseTest {

    private static final long[] testVals = {0, 1, 2, 8, 16, 42, 97, 100, 2894765};
    private static final long[] notAllTestVals = {0, 1, 2, 3, 7, 22, 61};
    private static final long notInTestVals = 6;
    private static final int RAND_SEED = 0xdeadf00;
    private static final int HHASH_NVALS = 10000;
    private static final int FPR_NVALS = 10000;
    private static final int N_TRIALS = 100;
    private static final float FPP = 0.01F;

    private static long randomLong(Random rng) {
        return (((long) rng.nextInt()) | (((long) rng.nextInt()) << 31)) & ~Long.MIN_VALUE;
    }

    @Test
    void createFromCollectionTest() {
        final LongBloomFilter bloomFilter = new LongBloomFilter(testVals.length, FPP);
        bloomFilter.addAll(testVals);
        Assert.assertTrue(bloomFilter.containsAll(testVals));
    }

    @Test
    void addTest() {
        LoggingUtils.setLoggingLevel(Log.LogLevel.DEBUG);
        final LongBloomFilter bloomFilter = new LongBloomFilter(testVals.length, FPP);
        for (final long val : testVals) {
            bloomFilter.add(val);
        }
        Assert.assertTrue(bloomFilter.containsAll(testVals));
    }

    @Test
    void clearTest() {
        final LongBloomFilter bloomFilter = new LongBloomFilter(testVals.length, FPP);
        bloomFilter.addAll(testVals);
        bloomFilter.add(1L);
        bloomFilter.clear();
        Assert.assertFalse(bloomFilter.contains(1L));
        Assert.assertTrue(bloomFilter.isEmpty());
    }

    @Test
    void containsTest() {
        final LongBloomFilter bloomFilter = new LongBloomFilter(testVals.length, FPP);
        bloomFilter.addAll(testVals);
        Assert.assertTrue(bloomFilter.containsAll(testVals));
        Assert.assertFalse(bloomFilter.contains(notInTestVals));
        Assert.assertFalse(bloomFilter.containsAll(notAllTestVals));
    }

    @Test
    void equalsAndHashcodeTest() {
        final LongBloomFilter bloomFilter1 = new LongBloomFilter(testVals.length, FPP);
        final LongBloomFilter bloomFilter2 = new LongBloomFilter(testVals.length, FPP);
        final LongBloomFilter bloomFilter3 = new LongBloomFilter(notAllTestVals.length, FPP);
        bloomFilter1.addAll(testVals);
        bloomFilter2.addAll(testVals);
        bloomFilter3.addAll(notAllTestVals);
        Assert.assertEquals(bloomFilter1.hashCode(), bloomFilter2.hashCode());
        Assert.assertEquals(bloomFilter1, bloomFilter2);
        Assert.assertNotEquals(bloomFilter1.hashCode(), bloomFilter3.hashCode());
        Assert.assertNotEquals(bloomFilter1, bloomFilter3);
    }

    @Test
    void getLegalSizeBelowTest() {
        Assert.assertEquals(LongBloomFilter.getLegalSizeBelow(1000),719);
        Assert.assertEquals(LongBloomFilter.getLegalSizeBelow(2147483630),2147483629);
    }

    @Test
    void loadRandomIntsTest() {
        final Random rng = new Random(RAND_SEED);
        for (int trialNo = 0; trialNo != N_TRIALS; ++trialNo) {
            final HashSet<Long> hashSet = new HashSet<>();
            final LongBloomFilter bloomFilter = new LongBloomFilter(HHASH_NVALS, FPP);
            final String trialMsg = "trialNo=" + trialNo;
            for (int valNo = 0; valNo != HHASH_NVALS; ++valNo) {
                final long randLong = randomLong(rng);
                bloomFilter.add(randLong);
                hashSet.add(new Long(randLong));
            }
            for (final Long val : hashSet) {
                Assert.assertTrue(bloomFilter.contains(val), trialMsg + ", testVal=" + val);
            }
            int num_false_pos = 0;
            int num_total = 0;
            for (int valNo = 0; valNo != FPR_NVALS; ++valNo) {
                final long randLong = randomLong(rng);
                if (!hashSet.contains(new Long(randLong))) {
                    num_total++;
                    if (bloomFilter.contains(randLong)) {
                        num_false_pos++;
                    }
                }
            }
            Assert.assertTrue(num_false_pos < num_total * FPP * 10);
        }
    }

    @Test
    void serializationTest() {
        final Random rng = new Random(RAND_SEED);
        final LongBloomFilter bloomFilter = new LongBloomFilter(HHASH_NVALS, FPP);
        final HashSet<Long> hashSet = new HashSet<>(HHASH_NVALS);
        for (int valNo = 0; valNo != HHASH_NVALS; ++valNo) {
            final long randLong = randomLong(rng);
            bloomFilter.add(randLong);
            hashSet.add(randLong);
        }

        final ByteArrayOutputStream bos = new ByteArrayOutputStream();
        final Output out = new Output(bos);
        final Kryo kryo = new Kryo();
        kryo.writeClassAndObject(out, bloomFilter);
        out.flush();

        final ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());
        final Input in = new Input(bis);
        @SuppressWarnings("unchecked")
        final LongBloomFilter bloomFilter2 = (LongBloomFilter) kryo.readClassAndObject(in);

        Assert.assertEquals(bloomFilter, bloomFilter2);
        for (Long val : hashSet) {
            Assert.assertTrue(bloomFilter2.contains(val));
        }
    }
}
