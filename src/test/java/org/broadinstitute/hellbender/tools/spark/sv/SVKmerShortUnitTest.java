package org.broadinstitute.hellbender.tools.spark.sv;

import org.broadinstitute.hellbender.utils.RandomDNA;
import org.testng.Assert;
import org.testng.TestException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

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

    @Test
    public void testMaskedKmers() {

        final int K = 31;
        final RandomDNA randDNA = new RandomDNA(379483L);
        final int seqLen = 1000;
        final byte[] baseBytes = randDNA.nextBases(seqLen);
        final String bases = new String(baseBytes);
        final Random rand = new Random(2738493L);
        for (int i = 0; i < seqLen - K; i++) {

            //Random k-mer
            String substr = bases.substring(i, i + K);

            //For generating the correct answer
            final StringBuilder sb = new StringBuilder(substr);

            //Generate random indices to mask, an even number of them so that the k-mer will still be odd
            List<Byte> maskIndices = new ArrayList<>(K);
            for (int j = 0; j < K; j++) {
                maskIndices.add(new Byte((byte) j));
            }
            Collections.shuffle(maskIndices);
            int maskSize = rand.nextInt(K - 3) + 2;
            if ((maskSize & 1) == 1) {
                maskSize--;
            }
            maskIndices = maskIndices.subList(0, maskSize);
            Collections.sort(maskIndices);
            final byte[] maskIndicesBytes = new byte[maskIndices.size()];
            for (int j = 0; j < maskIndices.size(); j++) {
                maskIndicesBytes[j] = maskIndices.get(j).byteValue();
                sb.setCharAt(maskIndicesBytes[j], 'X');
            }

            //Expected masked kmer
            final byte[] maskedBases = sb.toString().replaceAll("X", "").getBytes();
            final SVKmer expectedKmer = SVKmerizer.toKmer(maskedBases, new SVKmerShort(maskedBases.length));

            final SVKmer randKmer = SVKmerizer.toKmer(substr.getBytes(), new SVKmerShort(substr.length())).mask(maskIndicesBytes, K);
            Assert.assertEquals(randKmer.toString(K - maskIndicesBytes.length), expectedKmer.toString(maskedBases.length));

        }
    }

    @Test
    public void testMaskedKmerSet() {
        final int K = 31;
        final int seqLen = 100;
        final RandomDNA randDNA = new RandomDNA(379483L);
        final Random rand = new Random(2738493L);
        int numKmerTrue = 0;
        int numMaskTrue = 0;
        for (int trial = 0; trial < 10000; trial++) {
            final byte[] baseBytes = randDNA.nextBases(seqLen);
            final Iterator<SVKmer> kmerIter = SVKmerizer.stream(baseBytes, K, 1, new SVKmerShort(K)).iterator();
            final HashSet<Long> kmerLib = new HashSet<>(seqLen);
            while (kmerIter.hasNext()) {
                kmerLib.add(kmerIter.next().mask(new byte[0], K).getLong());
            }

            final byte[] mask = {0, 15};
            final Iterator<SVKmer> kmerMaskIter = SVKmerizer.stream(baseBytes, K, 1, new SVKmerShort(K)).iterator();
            final HashSet<Long> kmerMaskLib = new HashSet<>(seqLen);
            while (kmerMaskIter.hasNext()) {
                kmerMaskLib.add(kmerMaskIter.next().mask(mask, K).getLong());
            }
            final byte[] baseBytesMut = Arrays.copyOf(baseBytes, baseBytes.length);
            for (int i = 0; i < baseBytesMut.length; i++) {
                if (rand.nextDouble() < 0.035) {
                    switch (rand.nextInt(4)) {
                        case 0:
                            baseBytesMut[i] = 'A';
                            break;
                        case 1:
                            baseBytesMut[i] = 'T';
                            break;
                        case 2:
                            baseBytesMut[i] = 'C';
                            break;
                        case 3:
                            baseBytesMut[i] = 'G';
                            break;
                        default:
                            throw new TestException("Unrecognized base");
                    }
                }
            }
            final Iterator<SVKmer> mutIter = SVKmerizer.stream(baseBytesMut, K, 1, new SVKmerShort(K)).iterator();
            boolean bFoundKmer = false;
            while (mutIter.hasNext()) {
                if (kmerLib.contains(mutIter.next().getLong())) {
                    bFoundKmer = true;
                    break;
                }
            }
            if (bFoundKmer) numKmerTrue++;

            final Iterator<SVKmer> mutMaskIter = SVKmerizer.stream(baseBytesMut, K, 1, new SVKmerShort(K)).iterator();
            boolean bFoundMask = false;
            while (mutMaskIter.hasNext()) {
                if (kmerMaskLib.contains(mutMaskIter.next().mask(mask, K).getLong())) {
                    bFoundMask = true;
                    break;
                }
            }
            if (bFoundMask) numMaskTrue++;
        }
        System.out.println(numKmerTrue);
        System.out.println(numMaskTrue);
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

    @Test
    public void testMask() {
        final RandomDNA randDNA = new RandomDNA(4738389L);
        final Random rand = new Random(72039493L);
        final int kSize = 31;
        final byte[] bases = randDNA.nextBases(1000);
        for (int index = 0; index < bases.length - kSize + 1; index++) {
            final String kmerString = new String(Arrays.copyOfRange(bases,index,index+kSize));
            final StringBuilder sb = new StringBuilder(kSize);
            sb.append(kmerString);
            int maskSize = 2 + rand.nextInt(kSize-1);
            maskSize = maskSize - (maskSize % 2);
            List<Byte> maskList = new ArrayList<>(kSize);
            for (int i = 0; i < kSize; i++) {
                maskList.add((byte)i);
            }
            Collections.shuffle(maskList);
            maskList = maskList.subList(0,maskSize);
            Collections.sort(maskList);
            final byte[] mask = new byte[maskSize];
            for (int i = 0; i < maskSize; i++) {
                mask[i] = maskList.get(i);
                sb.deleteCharAt(mask[i] - i);
            }
            final SVKmer kmer = SVKmerizer.toKmer(kmerString.getBytes(), new SVKmerShort(kSize));
            final SVKmer kmerMasked = SVKmerizer.toKmer(sb.toString(), new SVKmerShort(kSize - mask.length));
            Assert.assertTrue(kmer.mask(mask, kSize).equals(kmerMasked));
        }
    }
}