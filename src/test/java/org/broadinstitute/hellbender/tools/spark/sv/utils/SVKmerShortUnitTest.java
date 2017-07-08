package org.broadinstitute.hellbender.tools.spark.sv.utils;

import org.broadinstitute.hellbender.utils.RandomDNA;
import org.testng.Assert;
import org.testng.TestException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

/**
 * Unit tests for SVKmerShort and SVKmerizer<SVKmerShort>.
 */
public class SVKmerShortUnitTest {
    @Test
    public void testDefaultConstruction() {
        Assert.assertEquals(new SVKmerShort(10).toString(10), "AAAAAAAAAA");
        Assert.assertEquals(new SVKmerShort(11).toString(11), "AAAAAAAAAAA");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testDefaultConstructionWithTooLargeK() {
        final SVKmerShort tooBigK = new SVKmerShort(32);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testDefaultConstructionWithTooSmallK() {
        final SVKmerShort tooSmallK = new SVKmerShort(0);
    }

    @DataProvider(name = "sequenceStrings")
    public Object[][] getSequenceStrings() {
        return new Object[][]{
                {"ACGTACGTACGT"}, {"ACGTACGTACGTC"}, {"ACGTACGTACGTCC"}, {"ACGTACGTACGTCCC"}, {"ACGTACGTACGTCCCC"}
        };
    }

    @Test(dataProvider = "sequenceStrings")
    public void testConstructionAndToString(final String str) {
        Assert.assertEquals(str, SVKmerizer.toKmer(str, new SVKmerShort(str.length())).toString(str.length()));
        Assert.assertEquals(str, SVKmerizer.toKmer(str.getBytes(), new SVKmerShort(str.length())).toString(str.length()));
    }

    @Test(dataProvider = "sequenceStrings")
    public void testSuccessor(final String str) {
        final int K = str.length();
        final StringBuilder sb = new StringBuilder(str);
        sb.delete(0, 1);
        sb.append('A');
        final SVKmer kkk = SVKmerizer.toKmer(str, new SVKmerShort(str.length()));
        Assert.assertEquals(sb.toString(), kkk.successor(SVKmerLong.Base.A, K).toString(K));
        sb.setCharAt(K - 1, 'C');
        Assert.assertEquals(sb.toString(), kkk.successor(SVKmerLong.Base.C, K).toString(K));
        sb.setCharAt(K - 1, 'G');
        Assert.assertEquals(sb.toString(), kkk.successor(SVKmerLong.Base.G, K).toString(K));
        sb.setCharAt(K - 1, 'T');
        Assert.assertEquals(sb.toString(), kkk.successor(SVKmerLong.Base.T, K).toString(K));
    }

    @DataProvider(name = "maskData")
    public Object[][] getMaskData() {
        return new Object[][]{
                {"AGCTAGCTGCAGTCTGCTAG", "AGCTAGCTGCAGTCTGCTAG", new byte[0], Boolean.TRUE},
                {"AGCTAGCTGCAGTCTGCTAG", "AGCTAGCTGCTGTCTGCTAG", new byte[0], Boolean.FALSE},
                {"AGCTAGCTGCAGTCTGCTAG", "AGCTAGTTGCAGTCTGCTAG", new byte[]{6}, Boolean.TRUE},
                {"AGCTAGCTGCAGTCTGCTAG", "AGCTAGTTGCAGTCTGCTAG", new byte[]{5}, Boolean.FALSE},
                {"AGCTAGCTGCAGTCTGCTAG", "AGCTAGTTGCAGTCTGCTAG", new byte[]{7}, Boolean.FALSE}
        };
    }

    @Test(dataProvider = "maskData")
    public void testMaskedKmers(final String seq1, final String seq2, final byte[] maskBytes, final boolean equal) {

        if (seq1.length() != seq2.length()) {
            throw new TestException("Kmer sequences must be of same length");
        }
        final int K = seq1.length();

        final SVKmerShort kmer1 = (SVKmerShort)SVKmerizer.toKmer(seq1, new SVKmerShort(K));
        final SVKmerShort kmer2 = (SVKmerShort)SVKmerizer.toKmer(seq2, new SVKmerShort(K));
        final SVKmerShort mask = SVKmerShort.getMask(maskBytes, K);

        if (equal) {
            Assert.assertEquals(kmer1.mask(mask), kmer2.mask(mask));
        } else {
            Assert.assertNotEquals(kmer1.mask(mask), kmer2.mask(mask));
        }
    }

    @Test
    public void testRandomMaskedKmers() {

        final int K = 31;
        final RandomDNA randDNA = new RandomDNA(379483L);
        final int seqLen = 1000;
        final byte[] baseBytes = randDNA.nextBases(seqLen);
        final String bases = new String(baseBytes);
        final Random rand = new Random(2738493L);
        for (int i = 0; i < seqLen - K; i++) {

            //Generate random indices to mask
            List<Byte> maskIndices = new ArrayList<>(K);
            for (int j = 0; j < K; j++) {
                maskIndices.add(new Byte((byte) j));
            }
            Collections.shuffle(maskIndices);
            int maskSize = rand.nextInt(K - 3) + 2;
            maskIndices = maskIndices.subList(0, maskSize);
            final byte[] maskIndicesBytes = new byte[maskIndices.size()];
            for (int j = 0; j < maskIndices.size(); j++) {
                maskIndicesBytes[j] = maskIndices.get(j).byteValue();
            }
            final SVKmerShort mask = SVKmerShort.getMask(maskIndicesBytes, K);

            //Create distinct kmers
            String substr = bases.substring(i, i + K);
            final byte[] maskedBasesA = substr.getBytes();
            final byte[] maskedBasesT = substr.getBytes();
            for (int j = 0; j < maskIndicesBytes.length; j++) {
                maskedBasesA[maskIndicesBytes[j]] = 'A';
                maskedBasesT[maskIndicesBytes[j]] = 'T';
            }
            final SVKmerShort kmerA = (SVKmerShort)SVKmerizer.toKmer(maskedBasesA, new SVKmerShort(K));
            final SVKmerShort kmerT = (SVKmerShort)SVKmerizer.toKmer(maskedBasesT, new SVKmerShort(K));

            Assert.assertEquals(kmerA.mask(mask), kmerT.mask(mask));
        }
    }

    @Test(dataProvider = "sequenceStrings")
    public void testPredecessor(final String str) {
        final int K = str.length();
        final StringBuilder sb = new StringBuilder(str);
        sb.insert(0, 'A');
        sb.setLength(K);
        final SVKmer kkk = SVKmerizer.toKmer(str, new SVKmerShort(str.length()));
        Assert.assertEquals(sb.toString(), kkk.predecessor(SVKmerLong.Base.A, K).toString(K));
        sb.setCharAt(0, 'C');
        Assert.assertEquals(sb.toString(), kkk.predecessor(SVKmerLong.Base.C, K).toString(K));
        sb.setCharAt(0, 'G');
        Assert.assertEquals(sb.toString(), kkk.predecessor(SVKmerLong.Base.G, K).toString(K));
        sb.setCharAt(0, 'T');
        Assert.assertEquals(sb.toString(), kkk.predecessor(SVKmerLong.Base.T, K).toString(K));
    }

    @Test
    public void testReverseComplementation() {
        Assert.assertEquals(SVKmerizer.toKmer("ACGTACGA", new SVKmerShort(8)).reverseComplement(8), SVKmerizer.toKmer("TCGTACGT", new SVKmerShort(8)));
        Assert.assertEquals(SVKmerizer.toKmer("ACGTACGTC", new SVKmerShort(9)).reverseComplement(9), SVKmerizer.toKmer("GACGTACGT", new SVKmerShort(9)));
        Assert.assertEquals(SVKmerizer.toKmer("ACGTCCGTC", new SVKmerShort(9)).reverseComplement(9), SVKmerizer.toKmer("GACGGACGT", new SVKmerShort(9)));
        Assert.assertEquals(SVKmerizer.toKmer("ACGTGCGTC", new SVKmerShort(9)).reverseComplement(9), SVKmerizer.toKmer("GACGCACGT", new SVKmerShort(9)));
        Assert.assertEquals(SVKmerizer.toKmer("ACGTTCGTC", new SVKmerShort(9)).reverseComplement(9), SVKmerizer.toKmer("GACGAACGT", new SVKmerShort(9)));
    }

    @Test
    public void testCanonicalization() {
        Assert.assertEquals(SVKmerizer.toKmer("ACGTACGTC", new SVKmerShort(9)).canonical(9), SVKmerizer.toKmer("ACGTACGTC", new SVKmerShort(9)));
        Assert.assertEquals(SVKmerizer.toKmer("GACGTACGT", new SVKmerShort(9)).canonical(9), SVKmerizer.toKmer("ACGTACGTC", new SVKmerShort(9)));
    }

    @Test
    public void testComparison() {
        final SVKmerShort kkk1 = (SVKmerShort) SVKmerizer.toKmer("ACGTA", new SVKmerShort(5));
        final SVKmerShort kkk2 = (SVKmerShort) SVKmerizer.toKmer("ACGTC", new SVKmerShort(5));
        Assert.assertTrue(kkk1.compareTo(kkk1) == 0);
        Assert.assertTrue(kkk2.compareTo(kkk2) == 0);
        Assert.assertTrue(kkk1.compareTo(kkk2) < 0);
        Assert.assertTrue(kkk2.compareTo(kkk1) > 0);
    }

    @Test
    public void testHashCode() {
        Assert.assertNotEquals(SVKmerizer.toKmer("TAGCGTA", new SVKmerShort(7)).hashCode(), SVKmerizer.toKmer("TAGCGTC", new SVKmerShort(7)).hashCode());
        Assert.assertEquals(SVKmerizer.toKmer("TAGGGTC", new SVKmerShort(7)).hashCode(), SVKmerizer.toKmer("TAGGGTC", new SVKmerShort(7)).hashCode());
    }

    @Test
    public void testKmerization() {
        final SVKmerizer kmerizer = new SVKmerizer("AAAAATT", 5, 1, new SVKmerShort(7));
        Assert.assertTrue(kmerizer.hasNext());
        Assert.assertEquals(kmerizer.next(), SVKmerizer.toKmer("AAAAA", new SVKmerShort(5)));
        Assert.assertTrue(kmerizer.hasNext());
        Assert.assertEquals(kmerizer.next(), SVKmerizer.toKmer("AAAAT", new SVKmerShort(5)));
        Assert.assertTrue(kmerizer.hasNext());
        Assert.assertEquals(kmerizer.next(), SVKmerizer.toKmer("AAATT", new SVKmerShort(5)));
        Assert.assertTrue(!kmerizer.hasNext());
    }

    @Test
    public void testKmerizationAcrossN() {
        final SVKmerizer kmerizer = new SVKmerizer("AAAAANTTTTT", 5, 1, new SVKmerShort(11));
        Assert.assertTrue(kmerizer.hasNext());
        Assert.assertEquals(kmerizer.next(), SVKmerizer.toKmer("AAAAA", new SVKmerShort(5)));
        Assert.assertTrue(kmerizer.hasNext());
        Assert.assertEquals(kmerizer.next(), SVKmerizer.toKmer("TTTTT", new SVKmerShort(5)));
        Assert.assertTrue(!kmerizer.hasNext());
    }
}