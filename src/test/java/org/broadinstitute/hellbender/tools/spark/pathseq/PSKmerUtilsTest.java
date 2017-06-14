package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.datasources.ReferenceFileSource;
import org.broadinstitute.hellbender.tools.spark.sv.SVKmerShort;
import org.broadinstitute.hellbender.tools.spark.utils.LargeLongHopscotchSet;
import org.broadinstitute.hellbender.tools.spark.utils.LongBloomFilter;
import org.broadinstitute.hellbender.tools.spark.utils.LongIterator;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public class PSKmerUtilsTest extends CommandLineProgramTest {

    final static long SEED = 6383943L;

    @Test
    public void testGetKmersFromReference() throws IOException {
        final ReferenceFileSource ref = new ReferenceFileSource( hg19MiniReference);
        final Collection<long[]> longCollection = PSKmerUtils.getMaskedKmersFromLocalReference(ref, 31, 1, SVKmerShort.getMask(new byte[0], 31));

        Assert.assertNotNull(longCollection);
        Assert.assertFalse(longCollection.isEmpty());
        int totalKmers = 0;
        for (final long[] list : longCollection) {
            Assert.assertNotNull(list);
            totalKmers += list.length;
        }

        final Set<Long> longSet = new HashSet<>(totalKmers);
        for (final long[] list : longCollection) {
            for (final long val : list) {
                longSet.add(val);
            }
        }

        final String expectedFile = "src/test/resources/" + PathSeqBuildKmers.class.getPackage().getName().replace(".", "/") + "/hg19mini.hss";

        final Input expectedInput = new Input(FileUtils.openInputStream(new File(expectedFile)));
        final Kryo kryo = new Kryo();
        final PSKmerSet expectedKmerLib = kryo.readObject(expectedInput, PSKmerSet.class);

        Assert.assertEquals(longSet.size(),expectedKmerLib.setSize());
        final Iterator<Long> iter = longSet.iterator();
        while (iter.hasNext()) {
            Assert.assertTrue(expectedKmerLib.contains(new SVKmerShort(iter.next())));
        }
    }

    @Test
    public void testLongArrayCollectionSize() {
        final Collection<long[]> longLists = new ArrayList<>(10);
        Assert.assertEquals(PSKmerUtils.longArrayCollectionSize(longLists), 0L);
        longLists.add(new long[]{1L, 2L, 3L, 4L});
        Assert.assertEquals(PSKmerUtils.longArrayCollectionSize(longLists), 4L);
        longLists.add(new long[]{1L, 2L,});
        Assert.assertEquals(PSKmerUtils.longArrayCollectionSize(longLists), 6L);
        longLists.add(new long[]{1L, 2L, 3L});
        Assert.assertEquals(PSKmerUtils.longArrayCollectionSize(longLists), 9L);
        longLists.add(new long[]{1L});
        Assert.assertEquals(PSKmerUtils.longArrayCollectionSize(longLists), 10L);
        longLists.add(new long[]{});
        Assert.assertEquals(PSKmerUtils.longArrayCollectionSize(longLists), 10L);
    }

    @Test
    public void testLongArrayCollectionToSet() {
        final int numLists = 1000;
        final int listSize = 1000;
        final Collection<long[]> longArrays = new ArrayList<>(numLists);
        final Set<Long> truthSet = new HashSet<>(numLists * listSize);
        final Random rand = new Random(SEED);
        for (int i = 0; i < listSize; i++) {
            final List<Long> longList = rand.longs(listSize).map(val -> val >>> 2).boxed().collect(Collectors.toList());
            longArrays.add(longList.stream().mapToLong(Long::longValue).toArray());
            truthSet.addAll(longList.stream().collect(Collectors.toList()));
        }
        final LargeLongHopscotchSet result = PSKmerUtils.longArrayCollectionToSet(longArrays, PSKmerUtils.longArrayCollectionSize(longArrays));

        Assert.assertEquals(result.size(), truthSet.size());
        Assert.assertTrue(result.containsAll(truthSet.stream().mapToLong(val -> val.longValue()).toArray()));
    }

    @Test
    public void testLongArrayCollectionToBloomFilter() {
        final int numLists = 1000;
        final int listSize = 1000;
        final Collection<long[]> longArrays = new ArrayList<>(numLists);
        final Set<Long> truthSet = new HashSet<>(numLists * listSize);
        final Random rand = new Random(SEED);
        for (int i = 0; i < listSize; i++) {
            final List<Long> longList = rand.longs(listSize).map(val -> val >>> 2).boxed().collect(Collectors.toList());
            longArrays.add(longList.stream().mapToLong(Long::longValue).toArray());
            truthSet.addAll(longList.stream().collect(Collectors.toList()));
        }

        final double fpp = 0.5;
        final LongBloomFilter result = PSKmerUtils.longArrayCollectionToBloomFilter(longArrays, PSKmerUtils.longArrayCollectionSize(longArrays), fpp);

        final Iterator<Long> itr = truthSet.iterator();
        while (itr.hasNext()) {
            Assert.assertTrue(result.contains(itr.next()));
        }

        final int numTrials = 1000;
        int falsePositives = 0;
        for (int i = 0; i < numTrials; i++) {
            final long randomVal = rand.nextLong();
            if (result.contains(randomVal >>> 2) && !truthSet.contains(randomVal)) {
                falsePositives++;
            }
        }
        final double testFPR = falsePositives / (double) numTrials;
        Assert.assertTrue(testFPR < fpp * 1.1, "New Bloom filter gave too many false positives");
        Assert.assertTrue(testFPR > fpp * 0.5, "New Bloom filter gave too few false positives");
    }

    @SuppressWarnings("unchecked")
    @Test
    public void testReadWriteSets() {
        final long numElements = 1024 * 1024L;
        final int kSize = 31;
        final SVKmerShort mask = SVKmerShort.getMask(new byte[]{3, 20, 25}, kSize);

        final LargeLongHopscotchSet hssUnmasked = new LargeLongHopscotchSet(numElements);
        final Random rand = new Random(738489373L);
        for (long i = 0; i < numElements; i++) {
            hssUnmasked.add(rand.nextLong() >>> 2);
        }

        final LargeLongHopscotchSet hssMasked = new LargeLongHopscotchSet(numElements);
        final LongIterator hssUnmaskedIter = hssUnmasked.iterator();
        while (hssUnmaskedIter.hasNext()) {
            hssMasked.add(PSKmerCollection.canonicalizeAndMask(new SVKmerShort(hssUnmaskedIter.next()), kSize, mask));
        }

        final File hssFile = createTempFile("set", ".bin");
        final PSKmerSet truthSet = new PSKmerSet(hssMasked, kSize, mask);
        PSKmerUtils.writeKmerSet(hssFile.getPath(), truthSet);

        final PSKmerSet hssIn = (PSKmerSet) PSKmerUtils.readKmerFilter(hssFile.getPath() + PSKmerUtils.HOPSCOTCH_SET_EXTENSION);
        Assert.assertEquals(hssIn, truthSet, "Hopscotch set changed after writing/reading");

        final double bloomFPP = 0.5;
        final LongBloomFilter bfOut = new LongBloomFilter(numElements, bloomFPP);
        LongIterator hssIter = hssMasked.iterator();
        while (hssIter.hasNext()) {
            bfOut.add(hssIter.next());
        }

        final File bfFile = createTempFile("set", ".bin");
        PSKmerUtils.writeKmerBloomFilter(bfFile.getPath(), new PSKmerBloomFilter(bfOut, kSize, mask, 1000));

        final PSKmerCollection bloomIn = PSKmerUtils.readKmerFilter(bfFile.getPath() + PSKmerUtils.BLOOM_FILTER_EXTENSION);
        hssIter = hssMasked.iterator();
        while (hssIter.hasNext()) {
            Assert.assertTrue(bloomIn.contains(new SVKmerShort(hssIter.next())), "Bloom filter changed after writing/reading");
        }
        for (int i = 0; i < 10000; i++) {
            final long val = rand.nextLong() >>> 2;
            Assert.assertEquals(bloomIn.contains(new SVKmerShort(val)), bfOut.contains(PSKmerCollection.canonicalizeAndMask(new SVKmerShort(val), kSize, mask)), "Bloom filter changed after writing/reading");
        }
    }

}