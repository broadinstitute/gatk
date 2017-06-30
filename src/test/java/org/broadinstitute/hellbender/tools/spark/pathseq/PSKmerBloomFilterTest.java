package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.broadinstitute.hellbender.tools.spark.sv.SVKmerShort;
import org.broadinstitute.hellbender.tools.spark.utils.LongBloomFilter;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.TestException;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;

import java.io.*;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Random;
import java.util.Set;

public class PSKmerBloomFilterTest extends BaseTest {

    private final double falsePositiveProb = 0.1;
    private final int setSize = 100000;
    private final long seed = 48393943L;
    private final int kmerSize = 31;
    private final SVKmerShort mask = SVKmerShort.getMask(new byte[]{15}, kmerSize);

    private PSKmerBloomFilter kmerSet;

    private PSKmerBloomFilter createRandomSet(final long randomSeed) {
        final LongBloomFilter longSet = new LongBloomFilter(setSize, falsePositiveProb);
        final Random rand = new Random(randomSeed);
        for (int i = 0; i < setSize; i++) {
            longSet.add(rand.nextLong() >>> 2);
        }
        return new PSKmerBloomFilter(longSet, kmerSize, mask, 10000);
    }

    @BeforeTest
    private void setup() {
        kmerSet = createRandomSet(seed);
    }

    @Test
    public void test() {

        final int kSize = 31;
        final LongBloomFilter bloomFilter = new LongBloomFilter(setSize, falsePositiveProb);
        final Set<Long> longSet = new HashSet<>(setSize);
        final SVKmerShort testMask = SVKmerShort.getMask(new byte[]{2, 18}, kSize);

        final Random rand = new Random(93848383L);
        for (int i = 0; i < setSize; i++) {
            final long val = rand.nextLong() >>> 2;
            longSet.add(PSKmerCollection.canonicalizeAndMask(new SVKmerShort(val), kSize, testMask));
        }
        for (final Long val : longSet) {
            bloomFilter.add(val);
        }

        final PSKmerBloomFilter testKmerSet = new PSKmerBloomFilter(bloomFilter, kSize, testMask, setSize);

        Assert.assertEquals(testKmerSet.getMask(), testMask);
        Assert.assertEquals(testKmerSet.kmerSize(), kSize);

        final Iterator<Long> iter = longSet.iterator();
        while (iter.hasNext()) {
            Assert.assertTrue(testKmerSet.contains(new SVKmerShort(iter.next())));
        }
        Assert.assertTrue(testKmerSet.getFalsePositiveProbability() <= falsePositiveProb);
        Assert.assertTrue(testKmerSet.getFalsePositiveProbability() > 0.5 * falsePositiveProb);
    }

    @Test
    public void testSerializeDeserialize() {

        final Kryo kryo = new Kryo();

        final File tempFile = createTempFile("serializeBloom", ".bfi");
        try (final OutputStream outputStream = new FileOutputStream(tempFile)) {
            final Output output = new Output(outputStream);
            kryo.writeObject(output, kmerSet);
            output.close();
        } catch (final IOException e) {
            throw new TestException(e.getMessage());
        }

        try (final InputStream inputStream = new FileInputStream(tempFile)) {
            final Input input = new Input(inputStream);
            final PSKmerBloomFilter testSet = kryo.readObject(input, PSKmerBloomFilter.class);
            Assert.assertEquals(kmerSet, testSet);
        } catch (final IOException e) {
            throw new TestException(e.getMessage());
        }
    }

    @Test
    public void testHashCodeAndEquals() {

        final PSKmerBloomFilter kmerSetEqual = createRandomSet(seed);
        final PSKmerBloomFilter kmerSet2 = createRandomSet(10593934L);

        Assert.assertEquals(kmerSetEqual, kmerSet);
        Assert.assertNotEquals(kmerSet2, kmerSet);

        Assert.assertEquals(kmerSetEqual.hashCode(), kmerSet.hashCode());
        Assert.assertNotEquals(kmerSet2.hashCode(), kmerSet.hashCode());
    }

}