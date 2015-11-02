package org.broadinstitute.hellbender.tools.spark.sv;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

/**
 * Unit tests for SVKmer and SVKmerizer.
 */
public class SVKmerUnitTest extends BaseTest {
    @Test
    public void testDefaultConstruction() {
        Assert.assertEquals(new SVKmer(10).toString(10), "AAAAAAAAAA");
        Assert.assertEquals(new SVKmer(11).toString(11), "AAAAAAAAAAA");
    }

    @DataProvider(name = "sequenceStrings")
    public Object[][] getSequenceStrings() {
        return new Object[][] {
                {"ACGTACGTACGT"}, {"ACGTACGTACGTC"}, {"ACGTACGTACGTCC"}, {"ACGTACGTACGTCCC"}, {"ACGTACGTACGTCCCC"}
        };
    }

    @Test(dataProvider = "sequenceStrings")
    public void testConstructionAndToString( final String str ) {
        Assert.assertEquals(str, SVKmerizer.toKmer(str).toString(str.length()));
        Assert.assertEquals(str, SVKmerizer.toKmer(str.getBytes()).toString(str.length()));
    }

    @Test(dataProvider = "sequenceStrings")
    public void testSuccessor( final String str ) {
        final int K = str.length();
        final StringBuilder sb = new StringBuilder(str);
        sb.delete(0, 1);
        sb.append('A');
        final SVKmer kkk = SVKmerizer.toKmer(str);
        Assert.assertEquals(sb.toString(), kkk.successor(SVKmer.Base.A, K).toString(K));
        sb.setCharAt(K-1, 'C');
        Assert.assertEquals(sb.toString(), kkk.successor(SVKmer.Base.C, K).toString(K));
        sb.setCharAt(K-1, 'G');
        Assert.assertEquals(sb.toString(), kkk.successor(SVKmer.Base.G, K).toString(K));
        sb.setCharAt(K-1, 'T');
        Assert.assertEquals(sb.toString(), kkk.successor(SVKmer.Base.T, K).toString(K));
    }

    @Test(dataProvider = "sequenceStrings")
    public void testPredecessor( final String str ) {
        final int K = str.length();
        final StringBuilder sb = new StringBuilder(str);
        sb.insert(0, 'A');
        sb.setLength(K);
        final SVKmer kkk = SVKmerizer.toKmer(str);
        Assert.assertEquals(sb.toString(), kkk.predecessor(SVKmer.Base.A, K).toString(K));
        sb.setCharAt(0, 'C');
        Assert.assertEquals(sb.toString(), kkk.predecessor(SVKmer.Base.C, K).toString(K));
        sb.setCharAt(0, 'G');
        Assert.assertEquals(sb.toString(), kkk.predecessor(SVKmer.Base.G, K).toString(K));
        sb.setCharAt(0, 'T');
        Assert.assertEquals(sb.toString(), kkk.predecessor(SVKmer.Base.T, K).toString(K));
    }

    @Test
    public void testReverseComplementation() {
        Assert.assertEquals(SVKmerizer.toKmer("ACGTACGA").reverseComplement(8), SVKmerizer.toKmer("TCGTACGT"));
        Assert.assertEquals(SVKmerizer.toKmer("ACGTACGTC").reverseComplement(9), SVKmerizer.toKmer("GACGTACGT"));
        Assert.assertEquals(SVKmerizer.toKmer("ACGTCCGTC").reverseComplement(9), SVKmerizer.toKmer("GACGGACGT"));
        Assert.assertEquals(SVKmerizer.toKmer("ACGTGCGTC").reverseComplement(9), SVKmerizer.toKmer("GACGCACGT"));
        Assert.assertEquals(SVKmerizer.toKmer("ACGTTCGTC").reverseComplement(9), SVKmerizer.toKmer("GACGAACGT"));
    }

    @Test
    public void testCanonicalization() {
        Assert.assertEquals(SVKmerizer.toKmer("ACGTACGTC").canonical(9), SVKmerizer.toKmer("ACGTACGTC"));
        Assert.assertEquals(SVKmerizer.toKmer("GACGTACGT").canonical(9), SVKmerizer.toKmer("ACGTACGTC"));
    }

    @Test
    public void testComparison() {
        final SVKmer kkk1 = SVKmerizer.toKmer("ACGTA");
        final SVKmer kkk2 = SVKmerizer.toKmer("ACGTC");
        Assert.assertTrue(kkk1.compareTo(kkk1) == 0);
        Assert.assertTrue(kkk2.compareTo(kkk2) == 0);
        Assert.assertTrue(kkk1.compareTo(kkk2) < 0);
        Assert.assertTrue(kkk2.compareTo(kkk1) > 0);
    }

    @Test
    public void testHashCode() {
        Assert.assertNotEquals(SVKmerizer.toKmer("TAGCGTA").hashCode(), SVKmerizer.toKmer("TAGCGTC").hashCode());
        Assert.assertEquals(SVKmerizer.toKmer("TAGGGTC").hashCode(), SVKmerizer.toKmer("TAGGGTC").hashCode());
    }

    @Test
    public void testKmerization() {
        final SVKmerizer kmerizer = new SVKmerizer("AAAAATT", 5);
        Assert.assertTrue(kmerizer.hasNext());
        Assert.assertEquals(kmerizer.next(), SVKmerizer.toKmer("AAAAA"));
        Assert.assertTrue(kmerizer.hasNext());
        Assert.assertEquals(kmerizer.next(), SVKmerizer.toKmer("AAAAT"));
        Assert.assertTrue(kmerizer.hasNext());
        Assert.assertEquals(kmerizer.next(), SVKmerizer.toKmer("AAATT"));
        Assert.assertTrue(!kmerizer.hasNext());
    }

    @Test
    public void testKmerizationAcrossN() {
        final SVKmerizer kmerizer = new SVKmerizer("AAAAANTTTTT", 5);
        Assert.assertTrue(kmerizer.hasNext());
        Assert.assertEquals(kmerizer.next(), SVKmerizer.toKmer("AAAAA"));
        Assert.assertTrue(kmerizer.hasNext());
        Assert.assertEquals(kmerizer.next(), SVKmerizer.toKmer("TTTTT"));
        Assert.assertTrue(!kmerizer.hasNext());
    }
}
