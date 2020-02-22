package org.broadinstitute.hellbender.utils.pairhmm;

import org.broadinstitute.hellbender.utils.RandomDNA;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;

public class DragstrLocusUnitTest {

    @Test()
    private void testWriteRead() throws IOException {
        final File binaryFile = File.createTempFile("test-dl", ".bin");
        final File textFile = File.createTempFile("test-dl", ".tab");
        //binaryFile.deleteOnExit();
        //textFile.deleteOnExit();
        final BinaryTableWriter<DragstrLocus> binaryWriter = DragstrLocus.binaryWriter(binaryFile);
        final TableWriter<DragstrLocus> textWriter = DragstrLocus.textWriter(new FileOutputStream(textFile));
        final List<DragstrLocus> loci = Arrays.stream(randomLoci()).map(oo -> (DragstrLocus) oo[0]).collect(Collectors.toList());
        binaryWriter.writeAll(loci);
        binaryWriter.close();
        textWriter.writeAllRecords(loci);
        textWriter.close();
        System.err.println(binaryFile);
        System.err.println(textFile);
        System.err.println(binaryFile.length());
        System.err.println(textFile.length());
        final BinaryTableReader<DragstrLocus> binaryReader = DragstrLocus.binaryReader(new FileInputStream(binaryFile));
        final List<DragstrLocus> loci2 = binaryReader.readAll();

        Assert.assertEquals(loci, loci2);
    }


    @DataProvider(name="randomLoci")
    public Object[][] randomLoci() {
        final int seed = 13;
        final Random rdn = new Random(seed);
        final RandomDNA randomDNA = new RandomDNA(rdn);
        final int chrCount = 25;
        final int[] lengths = new int[chrCount];
        lengths[0] = (int) Math.round(rdn.nextGaussian() * 10_000_000 + 250_000_000);
        for (int i = 1; i <chrCount; i++) {
            lengths[i] = (int) Math.round(rdn.nextGaussian() * lengths[i - 1] * 0.05 + lengths[i - 1] * .8);
        }
        final List<DragstrLocus> result = new ArrayList<>(10_000);
        for (int i = 0; i < 10_000; i++) {
            final int chrIdx = rdn.nextInt(chrCount);
            final int unitLength = rdn.nextInt(10) + 1;
            final byte[] unit = randomDNA.nextBases(unitLength);
            final int repeatCount = rdn.nextInt(40) + 1;
            final long start = rdn.nextInt(lengths[chrIdx]) + 1;
            final DragstrLocus locus = DragstrLocus.make(chrIdx, start, unit, repeatCount);
            result.add(locus);

        }
        return result.stream().map(dl -> new Object[] { dl }).toArray(Object[][]::new);
    }
}
