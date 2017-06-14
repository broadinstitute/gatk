package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.broadinstitute.hellbender.tools.spark.sv.SVKmerShort;
import org.broadinstitute.hellbender.tools.spark.utils.LargeLongHopscotchSet;
import org.broadinstitute.hellbender.tools.spark.utils.LongIterator;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.TestException;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;

import java.io.*;
import java.util.Random;

public class PSKmerSetTest extends BaseTest {

    private final int setSize = 1000;
    private final long seed = 48393943L;
    private final int kmerSize = 31;
    private final SVKmerShort mask = SVKmerShort.getMask(new byte[]{15}, kmerSize);

    private PSKmerSet kmerSet;

    private PSKmerSet createRandomSet(final long randomSeed) {
        final LargeLongHopscotchSet longSet = new LargeLongHopscotchSet(setSize);
        final Random rand = new Random(randomSeed);
        for (int i = 0; i < setSize; i++) {
            longSet.add(rand.nextLong() >>> 2);
        }
        return new PSKmerSet(longSet, kmerSize, mask);
    }

    @BeforeTest
    private void setup() {
        kmerSet = createRandomSet(seed);
    }

    @Test
    public void test() {
        final int kSize = 31;
        final LargeLongHopscotchSet testSet = new LargeLongHopscotchSet(4);
        final long[] longValues = new long[]{483L, 943L, 2L, 493L};
        testSet.addAll(longValues);
        final SVKmerShort testMask = SVKmerShort.getMask(new byte[]{}, kSize);
        final PSKmerSet testKmerSet = new PSKmerSet(testSet, kSize, testMask);

        Assert.assertEquals(testKmerSet.getMask(), testMask);
        Assert.assertEquals(testKmerSet.kmerSize(), kSize);
        Assert.assertEquals(testKmerSet.setSize(), testSet.size());

        LongIterator iter = testSet.iterator();
        while (iter.hasNext()) {
            Assert.assertTrue(testKmerSet.contains(new SVKmerShort(iter.next())));
        }

        iter = testKmerSet.iterator();
        while (iter.hasNext()) {
            Assert.assertTrue(testSet.contains(iter.next()));
        }
    }

    @Test
    public void testSerializeDeserialize() {

        final Kryo kryo = new Kryo();

        final File tempFile = createTempFile("serializeSet", ".hss");
        try (final OutputStream outputStream = new FileOutputStream(tempFile)) {
            final Output output = new Output(outputStream);
            kryo.writeObject(output, kmerSet);
            output.close();
        } catch (final IOException e) {
            throw new TestException(e.getMessage());
        }

        try (final InputStream inputStream = new FileInputStream(tempFile)) {
            final Input input = new Input(inputStream);
            final PSKmerSet testSet = kryo.readObject(input, PSKmerSet.class);
            Assert.assertEquals(kmerSet, testSet);
        } catch(final IOException e) {
            throw new TestException(e.getMessage());
        }
    }

    @Test
    public void testHashCodeAndEquals() {

        final PSKmerSet kmerSetEqual = createRandomSet(seed);
        final PSKmerSet kmerSet2 = createRandomSet(10593934L);

        Assert.assertEquals(kmerSetEqual, kmerSet);
        Assert.assertNotEquals(kmerSet2, kmerSet);

        Assert.assertEquals(kmerSetEqual.hashCode(), kmerSet.hashCode());
        Assert.assertNotEquals(kmerSet2.hashCode(), kmerSet.hashCode());
    }

}