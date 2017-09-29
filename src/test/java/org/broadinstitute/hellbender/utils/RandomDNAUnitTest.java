package org.broadinstitute.hellbender.utils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.stat.inference.ChiSquareTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.IntStream;

public final class RandomDNAUnitTest {

    private static final int TEST_BASES_PER_LINE = 73;

    @Test(dataProvider="dictionaries")
    public void testRandomFasta(final SAMSequenceDictionary dict) throws IOException {
        final RandomDNA randomDNA = new RandomDNA(111);
        File fastaFile = null;
        try {
            fastaFile = randomDNA.nextFasta(dict, TEST_BASES_PER_LINE);
            assertFastaFileAndDictMatch(fastaFile, TEST_BASES_PER_LINE, dict);
        } finally {
            try { if (fastaFile != null) fastaFile.delete(); } catch (final RuntimeException ex) {};
        }
    }

    private void assertFastaFileAndDictMatch(final File fastaFile, final int basesPerLine, final SAMSequenceDictionary dict) {
        try (final BufferedReader reader = new BufferedReader(new FileReader(fastaFile))) {
            final List<Pair<String, Nucleotide.Counter>> nameLengthAndFreqs = new ArrayList<>();
            String line = reader.readLine();
            if (dict.getSequences().isEmpty()) {
                Assert.assertNull(line);
            } else {
                Assert.assertNotNull(line);
                do {
                    Assert.assertTrue(line.matches("^>\\S.*$"));
                    final String name = line.substring(1).split("\\s+")[0];
                    final Nucleotide.Counter frequencies = new Nucleotide.Counter();
                    line = reader.readLine();
                    while (line != null && !line.matches("^>.*$")) {
                        final String lineBases = line.trim();
                        final String nextLine = reader.readLine();
                        for (final byte base : lineBases.getBytes()) {
                            final Nucleotide nuc = Nucleotide.valueOf(base);
                            Assert.assertTrue(nuc.isConcrete());
                            frequencies.add(nuc);
                        }
                        if (nextLine != null && !nextLine.matches("^>.*$")){
                            Assert.assertEquals(lineBases.length(), basesPerLine);
                        } else {
                            Assert.assertTrue(lineBases.length() <= basesPerLine);
                        }
                        line = nextLine;
                    }
                    nameLengthAndFreqs.add(new ImmutablePair<>(name, frequencies));
                } while (line != null);
                Assert.assertEquals(nameLengthAndFreqs.size(), dict.getSequences().size());
                for (int i = 0; i < nameLengthAndFreqs.size(); i++) {
                    Assert.assertEquals(nameLengthAndFreqs.get(i).getLeft(), dict.getSequence(i).getSequenceName());
                    Assert.assertEquals(nameLengthAndFreqs.get(i).getRight().sum(), dict.getSequence(i).getSequenceLength());
                }
            }
        } catch (final IOException ex) {
            Assert.fail("exception thrown when openning fastaFile", ex);
        }
    }

    private int[] counts(final byte[] bytes){
        final int[] b= new int[4];
        for(int i=0; i < bytes.length; i++){
            switch (bytes[i]){
                case 'A': b[0]++; break;
                case 'C': b[1]++; break;
                case 'G': b[2]++; break;
                case 'T': b[3]++; break;
                default: throw new IllegalStateException("illegal base:" + bytes[i]);
            }
        }
        return b;
    }
    @Test
    public void testBases1(){
        int[] results = new int[4];

        final int n = 1000;
        final int m = 13;
        for (int i= 0; i < n; i++) {
            final byte[] b = new RandomDNA().nextBases(m);
            final int[] b0 = counts(b);
            results = pairwiseAdd(results, b0);
        }

        checkResults(results, n, m);
    }

    @Test
    public void testBases(){
        int[] results = new int[4];

        final int n = 1000;
        final int m = 13;
        for (int i= 0; i < n; i++) {
            final byte[] b = new byte[m];
            new RandomDNA().nextBases(b);
            final int[] b0 = counts(b);
            results = pairwiseAdd(results, b0);
        }

        checkResults(results, n, m);
    }

    public void checkResults(final int[] results, final int n, final int m) {
        final double[] dresults = MathUtils.promote(results);
        final double mean = MathUtils.mean(dresults, 0, dresults.length);
        final double std = new StandardDeviation().evaluate(dresults);
        final double expectedMean = (n*m)/4.0;
        final double s = std; // not really because it's the population not the sample dtd but it'll do
        Assert.assertTrue(mean < expectedMean + 2 * s / Math.sqrt(n * m), "unexpected mean:" + mean);
        Assert.assertTrue(mean > expectedMean-2*s/Math.sqrt(n*m), "unexpected mean:" +mean);
    }

    public void assertUniformity(final Nucleotide.Counter freqs) {
        final long[] observed = new long[] { freqs.get(Nucleotide.A), freqs.get(Nucleotide.C), freqs.get(Nucleotide.G), freqs.get(Nucleotide.T)};
        final ChiSquareTest test = new ChiSquareTest();
        final double pValue = test.chiSquareTest(new double[] { 0.25, 0.25, 0.25, 0.25 }, observed );
        Assert.assertTrue( pValue > 0.001);
    }

    private int[] pairwiseAdd(int[] a, int[] b) {
        Utils.validateArg(a.length == b.length, "lengths must be equal");
        return IntStream.range(0, a.length).map(n -> a[n] + b[n]).toArray();
    }

    @DataProvider
    public Object[][] dictionaries() {
        return Arrays.asList(
                new SAMSequenceDictionary(),
                new SAMSequenceDictionary(Arrays.asList(
                        new SAMSequenceRecord("seq1", 1000))),
                new SAMSequenceDictionary(Arrays.asList(
                        new SAMSequenceRecord("chr20", 10_000),
                        new SAMSequenceRecord("chrX", 1_000),
                        new SAMSequenceRecord("MT_unknown", 10))),
                new SAMSequenceDictionary(Arrays.asList(
                        new SAMSequenceRecord("1", 0),
                        new SAMSequenceRecord("2", TEST_BASES_PER_LINE),
                        new SAMSequenceRecord("3", TEST_BASES_PER_LINE - 1),
                        new SAMSequenceRecord("mmm", TEST_BASES_PER_LINE + 1),
                        new SAMSequenceRecord("xxx", TEST_BASES_PER_LINE * 13)))).stream()
        .map(dict -> new Object[] { dict})
        .toArray(Object[][]::new);
    }
}
