package org.broadinstitute.hellbender.tools.dragstr;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.utils.RandomDNA;
import org.broadinstitute.hellbender.utils.BinaryTableReader;
import org.broadinstitute.hellbender.utils.BinaryTableWriter;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Unit test for {@code DragstrLocus}.
 */
public final class DragstrLocusUnitTest {

    private static final int SEED_1 = 13;
    private static final int SEED_2 = 31;
    private static final int TEST_DICTIONARY_CHR_COUNT = 25;
    private static final int TEST_INTERVAL_COUNT = 10_000;

    private static final SAMSequenceDictionary testDictionary = createTestDictionary();

    @Test
    private void testWriteRead() throws IOException {
        final File binaryFile = File.createTempFile("test-dl", ".bin");
        final File textFile = File.createTempFile("test-dl", ".tab");
        binaryFile.deleteOnExit();
        textFile.deleteOnExit();
        final BinaryTableWriter<DragstrLocus> binaryWriter = DragstrLocusUtils.binaryWriter(binaryFile);
        final TableWriter<DragstrLocus> textWriter = DragstrLocusUtils.textWriter(new FileOutputStream(textFile), testDictionary);
        final List<DragstrLocus> loci = randomLoci();
        binaryWriter.writeAll(loci);
        binaryWriter.close();
        textWriter.writeAllRecords(loci);
        textWriter.close();
        final BinaryTableReader<DragstrLocus> binaryReader = DragstrLocusUtils.binaryReader(new FileInputStream(binaryFile));
        final List<DragstrLocus> loci2 = binaryReader.readAll();

        Assert.assertEquals(loci, loci2);
    }

    private List<DragstrLocus> randomLoci() {
        final Random rdn = new Random(SEED_2);
        final RandomDNA randomDNA = new RandomDNA(rdn);
        final List<DragstrLocus> result = new ArrayList<>(TEST_INTERVAL_COUNT);
        for (int i = 0; i < TEST_INTERVAL_COUNT; i++) {
            final int chrIdx = rdn.nextInt(TEST_DICTIONARY_CHR_COUNT);
            final int unitLength = rdn.nextInt(10) + 1;
            final int repeatCount = rdn.nextInt(40) + 1;
            final long start = rdn.nextInt(testDictionary.getSequence(chrIdx).getSequenceLength()) + 1;
            final DragstrLocus locus = DragstrLocus.make(chrIdx, start, (byte) unitLength, (short)(repeatCount * unitLength), i);
            result.add(locus);
        }
        Collections.sort(result);
        return result;
    }

    private static SAMSequenceDictionary createTestDictionary() {

        final Random rdn = new Random(SEED_1);
        final int[] lengths = new int[TEST_DICTIONARY_CHR_COUNT];
        lengths[0] = (int) Math.round(rdn.nextGaussian() * 10_000_000 + 250_000_000);
        for (int i = 1; i < TEST_DICTIONARY_CHR_COUNT; i++) {
            lengths[i] = (int) Math.round(rdn.nextGaussian() * lengths[i - 1] * 0.05 + lengths[i - 1] * .8);
        }
        return new SAMSequenceDictionary(IntStream.range(0, TEST_DICTIONARY_CHR_COUNT)
                .mapToObj(idx -> new SAMSequenceRecord("chr" + (idx + 1), lengths[idx]))
                .collect(Collectors.toList()));
    }
}
