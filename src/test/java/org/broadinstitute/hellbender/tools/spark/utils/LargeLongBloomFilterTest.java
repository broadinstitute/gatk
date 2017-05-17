package org.broadinstitute.hellbender.tools.spark.utils;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.util.HashSet;
import java.util.Random;

public final class LargeLongBloomFilterTest extends BaseTest {

    private static final long[] testVals = {0, 1, 2, 8, 16, 42, 97, 100, 2894765};
    private static final long[] notAllTestVals = {0, 1, 2, 3, 7, 22, 61};
    private static final long notInTestVals = 6;
    private static final int RAND_SEED = 0xdeadf00;
    private static final long HHASH_NVALS = 100000;
    private static final int N_TRIALS = 10;
    private static final double FPP = 0.01;
    private static final int BYTES_PER_PARTITION = 10000;

    private static long randomLong(Random rng) {
        return (((long) rng.nextInt()) | (((long) rng.nextInt()) << 31)) & ~Long.MIN_VALUE;
    }

    @Test
    void getSetsTest() {
        final LargeLongBloomFilter bloomFilter = new LargeLongBloomFilter(BYTES_PER_PARTITION, 1, FPP);
        Assert.assertTrue(bloomFilter.getSets().size() == 1);
    }

    @Test
    void createFromCollectionTest() {
        final LargeLongBloomFilter bloomFilter = new LargeLongBloomFilter(BYTES_PER_PARTITION, testVals.length, FPP);
        bloomFilter.addAll(testVals);
        Assert.assertTrue(bloomFilter.containsAll(testVals));
    }

    @Test
    void addTest() {
        final LargeLongBloomFilter bloomFilter = new LargeLongBloomFilter(BYTES_PER_PARTITION, testVals.length, FPP);
        for (final long val : testVals) {
            bloomFilter.add(val);
        }
        Assert.assertTrue(bloomFilter.containsAll(testVals));
    }

    @Test
    void clearTest() {
        final LargeLongBloomFilter bloomFilter = new LargeLongBloomFilter(BYTES_PER_PARTITION, testVals.length, FPP);
        bloomFilter.addAll(testVals);
        bloomFilter.clear();
        Assert.assertTrue(bloomFilter.isEmpty());
        for (final long val : testVals) {
            Assert.assertFalse(bloomFilter.contains(val));
        }
    }

    @Test
    void containsTest() {
        final LargeLongBloomFilter bloomFilter = new LargeLongBloomFilter(BYTES_PER_PARTITION, testVals.length, FPP);
        bloomFilter.addAll(testVals);
        Assert.assertTrue(bloomFilter.containsAll(testVals));
        Assert.assertFalse(bloomFilter.contains(notInTestVals));

        //The probability of all 4 being FPs should be nearly 0
        Assert.assertFalse(bloomFilter.containsAll(notAllTestVals));
    }

    @Test
    public void testBloomFilterFPP() {
        final long numElements = 1024 * 1024L;
        final double fpp = 0.5;
        final LargeLongBloomFilter bf = new LargeLongBloomFilter(1024L, numElements, fpp);
        final Random rand = new Random(8387349L);
        for (int i = 0; i < numElements; i++) {
            bf.add(rand.nextLong() >>> 1);
        }
        final long numTrials = 10000;
        double fppTest = LargeLongBloomFilter.estimateFalsePositiveProbability(bf, numTrials, 8394983L);
        Assert.assertTrue(Math.abs(fpp - fppTest) < 0.1 * fpp);
    }

    @Test
    void loadRandomIntsTest() {
        final Random rng = new Random(RAND_SEED);
        for (int trialNo = 0; trialNo != N_TRIALS; ++trialNo) {
            final HashSet<Long> hashSet = new HashSet<>((int)HHASH_NVALS);
            final LargeLongBloomFilter bloomFilter = new LargeLongBloomFilter(BYTES_PER_PARTITION, HHASH_NVALS, FPP);
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
            for (int valNo = 0; valNo != HHASH_NVALS; ++valNo) {
                final long randLong = randomLong(rng);
                if (!hashSet.contains(new Long(randLong))) {
                    num_total++;
                    if (bloomFilter.contains(randLong)) {
                        num_false_pos++;
                    }
                }
            }
            final double fpr = ((double)num_false_pos)/num_total;
            Assert.assertTrue(fpr < FPP * 3);
        }
    }

    @Test
    void equalsAndHashcodeTest() {
        final LargeLongBloomFilter bloomFilter1 = new LargeLongBloomFilter(BYTES_PER_PARTITION, testVals.length, FPP);
        final LargeLongBloomFilter bloomFilter2 = new LargeLongBloomFilter(BYTES_PER_PARTITION, testVals.length, FPP);
        final LargeLongBloomFilter bloomFilter3 = new LargeLongBloomFilter(BYTES_PER_PARTITION, notAllTestVals.length, FPP);
        bloomFilter1.addAll(testVals);
        bloomFilter2.addAll(testVals);
        bloomFilter3.addAll(notAllTestVals);
        Assert.assertEquals(bloomFilter1.hashCode(), bloomFilter2.hashCode());
        Assert.assertEquals(bloomFilter1, bloomFilter2);
        Assert.assertNotEquals(bloomFilter1.hashCode(), bloomFilter3.hashCode());
        Assert.assertNotEquals(bloomFilter1, bloomFilter3);
    }

    @Test
    void serializationTest() {
        final Random rng = new Random(RAND_SEED);
        final LargeLongBloomFilter bloomFilter = new LargeLongBloomFilter(BYTES_PER_PARTITION, HHASH_NVALS, FPP);
        final HashSet<Long> hashSet = new HashSet<>((int)HHASH_NVALS);
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
        final LargeLongBloomFilter bloomFilter2 = (LargeLongBloomFilter) kryo.readClassAndObject(in);

        Assert.assertEquals(bloomFilter, bloomFilter2);
        for (Long val : hashSet) {
            Assert.assertTrue(bloomFilter2.contains(val));
        }
    }
}
