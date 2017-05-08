package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.commons.io.FileUtils;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.tools.spark.sv.SVKmerShort;
import org.broadinstitute.hellbender.tools.spark.utils.LargeLongBloomFilter;
import org.broadinstitute.hellbender.tools.spark.utils.LargeLongHopscotchSet;
import org.broadinstitute.hellbender.tools.spark.utils.LargeQueryableLongSet;
import org.broadinstitute.hellbender.tools.spark.utils.LongIterator;
import org.testng.Assert;
import org.testng.TestException;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.LongStream;

import static org.broadinstitute.hellbender.tools.spark.pathseq.PSKmerLibUtils.getKmersFromReference;

public class PSKmerLibUtilsTest extends CommandLineProgramTest {

    final static long MAX_PARTITION_BYTES = 1024 * 1024L;
    final static long SEED = 6383943L;
    final static int K_SIZE = 31;

    protected final static LargeLongHopscotchSet randomTestSet(final long size) {
        final LargeLongHopscotchSet kmerSet = new LargeLongHopscotchSet(MAX_PARTITION_BYTES, size);
        final Random rand = new Random(SEED);
        for (int i = 0; i < size; i++) {
            kmerSet.add(rand.nextLong() >>> 1);
        }
        return kmerSet;
    }

    @Test(groups="spark")
    public void testGetKmersFromReference() {
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        final ReferenceMultiSource ref = new ReferenceMultiSource((PipelineOptions)null, hg19MiniReference, ReferenceWindowFunctions.IDENTITY_FUNCTION);
        final SAMSequenceDictionary refDict = ref.getReferenceSequenceDictionary(null);
        final LargeLongHopscotchSet set = getKmersFromReference(ctx, 1024 * 1024, 31, 1, new byte[0], ref, null, refDict);

        Assert.assertNotNull(set);
        Assert.assertFalse(set.isEmpty());

        final String expectedFile = "src/test/resources/" + PathSeqKmerSpark.class.getPackage().getName().replace(".", "/") + "/PathSeqKmerSpark/kmer.hss";
        try {
            final Input input_expected = new Input(FileUtils.openInputStream(new File(expectedFile)));

            final Kryo kryo = new Kryo();
            kryo.setReferences(false);

            final LargeLongHopscotchSet expectedKmerLib = (LargeLongHopscotchSet) kryo.readClassAndObject(input_expected);

            Assert.assertEquals(expectedKmerLib.size(),set.size());
            LongIterator iter = set.iterator();
            while (iter.hasNext()) {
                Assert.assertTrue(expectedKmerLib.contains(iter.next()));
            }

        } catch (IOException e) {
            throw new TestException("Could not open kmers file " + expectedFile, e);
        }
    }

    @Test
    public void testLongArrayCollectionToKmerSet() {
        final int numKmerArrays = 1000;
        final int arraySize = 1000;
        final Collection<long[]> kmerLists = new ArrayList<>(numKmerArrays);
        final Set<Long> truthSet = new HashSet<>(numKmerArrays * arraySize);
        final Random rand = new Random(SEED);
        for (int i = 0; i < arraySize; i++) {
            final long[] kmers = rand.longs(arraySize).map(val -> val >>> 1).toArray();
            kmerLists.add(kmers);
            truthSet.addAll(LongStream.of(kmers).boxed().collect(Collectors.toList()));
        }
        final LargeLongHopscotchSet result = PSKmerLibUtils.longArrayCollectionToKmerSet(kmerLists, MAX_PARTITION_BYTES);

        Assert.assertEquals(result.size(), truthSet.size());
        Assert.assertTrue(result.containsAll(truthSet.stream().mapToLong(val -> val.longValue()).toArray()));
    }

    @Test
    public void testMaskKmers() {
        final int numKmers = 1000;
        final byte[] mask = {0, 15};
        final LargeLongHopscotchSet kmerSet = randomTestSet(numKmers);
        final LargeLongHopscotchSet result = PSKmerLibUtils.maskKmers(kmerSet, K_SIZE, mask, MAX_PARTITION_BYTES);
        LongIterator itr = kmerSet.iterator();
        while (itr.hasNext()) {
            final long maskedKmer = new SVKmerShort(itr.next()).mask(mask, K_SIZE).canonical(K_SIZE - mask.length).getLong();
            Assert.assertTrue(result.contains(maskedKmer));
        }

        final LargeLongHopscotchSet result_identical = PSKmerLibUtils.maskKmers(kmerSet, K_SIZE, new byte[0], MAX_PARTITION_BYTES);
        itr = kmerSet.iterator();
        while (itr.hasNext()) {
            Assert.assertTrue(result_identical.contains(itr.next()));
        }
    }

    @Test
    public void testDownsampleKmerSet() {
        final int numKmers = 10000;
        final LargeLongHopscotchSet kmerSet = randomTestSet(numKmers);
        final double downsampleProbability = 0.5;
        final LargeLongHopscotchSet result = PSKmerLibUtils.downsampleKmerSet(kmerSet, MAX_PARTITION_BYTES,
                downsampleProbability, 6384783L);
        Assert.assertTrue(result.size() > 0.8 * downsampleProbability * kmerSet.size(),
                "Downsampled set much smaller than expected");
        Assert.assertTrue(result.size() < 1.2 * downsampleProbability * kmerSet.size(),
                "Downsampled set much larger than expected");

        final LargeLongHopscotchSet result_identical = PSKmerLibUtils.downsampleKmerSet(kmerSet, MAX_PARTITION_BYTES,
                1.0, 6384783L);
        Assert.assertEquals(result_identical.size(), kmerSet.size(), "Downsampling should not have done anything");
    }

    @Test
    public void testCreateBloomFilterFromSet() {
        final int numKmers = 10000;
        final LargeLongHopscotchSet kmerSet = randomTestSet(numKmers);
        final double fpp = 0.5;
        final LargeLongBloomFilter result = PSKmerLibUtils.createBloomFilterFromSet(kmerSet, MAX_PARTITION_BYTES, fpp);
        final LongIterator itr = kmerSet.iterator();
        while (itr.hasNext()) {
            Assert.assertTrue(result.contains(itr.next()));
        }

        final Random rand = new Random(38483L + SEED * 7);
        final int numTrials = 1000;
        int numPositive = 0;
        for (int i = 0; i < numTrials; i++) {
            if (result.contains(rand.nextLong() >>> 1)) {
                numPositive++;
            }
        }
        final double testFPR = numPositive / (double) numTrials;
        Assert.assertTrue(testFPR < fpp * 1.1, "New Bloom filter gave too many false positives");
        Assert.assertTrue(testFPR > fpp * 0.5, "New Bloom filter gave too few false positives");
    }

    @Test
    public void testReadWriteSets() {
        final long maxParitionBytes = 1024 * 1024L;
        final long numElements = 1024 * 1024L;
        final LargeLongHopscotchSet hssOut = new LargeLongHopscotchSet(maxParitionBytes, numElements);
        final Random rand = new Random(738489373L);
        for (long i = 0; i < numElements; i++) {
            hssOut.add(rand.nextLong() >>> 1);
        }
        final File hssFile = createTempFile("set", ".hss");
        PSKmerLibUtils.writeLargeLongHopscotchSet(hssOut, hssFile.getPath(), null);

        final LargeLongHopscotchSet hssIn = PSKmerLibUtils.readLargeLongHopscotchSet(hssFile.getPath(), null);
        Assert.assertEquals(hssIn, hssOut, "Hopscotch set changed after writing/reading");

        final double bloomFPP = 0.5;
        final LargeLongBloomFilter bfOut = new LargeLongBloomFilter(maxParitionBytes, numElements, bloomFPP);
        LongIterator hssIter = hssOut.iterator();
        while (hssIter.hasNext()) {
            bfOut.add(hssIter.next());
        }
        final File bfFile = createTempFile("set", ".bfi");
        PSKmerLibUtils.writeLargeLongBloomFilter(bfOut, bfFile.getPath(), null);

        final LargeQueryableLongSet queryableIn = PSKmerLibUtils.readLargeQueryableSet(bfFile.getPath(), null);
        hssIter = hssOut.iterator();
        while (hssIter.hasNext()) {
            Assert.assertTrue(queryableIn.contains(hssIter.next()), "Bloom filter changed after writing/reading");
        }
        for (int i = 0; i < 10000; i++) {
            final long val = rand.nextLong() >>> 1;
            Assert.assertEquals(queryableIn.contains(val), bfOut.contains(val), "Bloom filter changed after writing/reading");
        }
    }

}