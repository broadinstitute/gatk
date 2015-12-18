package org.broadinstitute.hellbender.tools.spark.sv;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

/**
 * Created by tsharpe on 1/12/16.
 */
public class KmerTest extends BaseTest
{
    @Test
    public void testDefaultConstruction()
    {
        Assert.assertEquals(new Kmer(10).toString(10),"AAAAAAAAAA");
        Assert.assertEquals(new Kmer(11).toString(11),"AAAAAAAAAAA");
    }

    @Test
    public void testConstructionFromStringsAndBytes()
    {
        testConstruction("ACGTACGTACGT");
        testConstruction("ACGTACGTACGTC");
        testConstruction("ACGTACGTACGTCC");
        testConstruction("ACGTACGTACGTCCC");
        testConstruction("ACGTACGTACGTCCCC");
    }

    private void testConstruction( final String str )
    {
        Assert.assertEquals(str,StringKmerizer.toKmer(str).toString(str.length()));
        Assert.assertEquals(str,BytesKmerizer.toKmer(str.getBytes()).toString(str.length()));
    }

    @Test
    public void testSuccessor()
    {
        testSuccessor("ACGTTGCA");
        testSuccessor("ACGTATGCA");
    }

    private void testSuccessor( final String str )
    {
        final int K = str.length();
        final StringBuilder sb = new StringBuilder(str);
        sb.delete(0,1);
        sb.append('A');
        Kmer kkk = StringKmerizer.toKmer(str);
        Assert.assertEquals(sb.toString(),kkk.successor(0,K).toString(K));
        sb.setCharAt(K-1,'C');
        Assert.assertEquals(sb.toString(),kkk.successor(1,K).toString(K));
        sb.setCharAt(K-1,'G');
        Assert.assertEquals(sb.toString(),kkk.successor(2,K).toString(K));
        sb.setCharAt(K-1,'T');
        Assert.assertEquals(sb.toString(),kkk.successor(3,K).toString(K));
    }

    @Test
    public void testPredecessor()
    {
        testPredecessor("ACGTTGCA");
        testPredecessor("ACGTATGCA");
    }

    private void testPredecessor( final String str )
    {
        final int K = str.length();
        final StringBuilder sb = new StringBuilder(str);
        sb.insert(0,'A');
        sb.setLength(K);
        Kmer kkk = StringKmerizer.toKmer(str);
        Assert.assertEquals(sb.toString(),kkk.predecessor(0,K).toString(K));
        sb.setCharAt(0,'C');
        Assert.assertEquals(sb.toString(),kkk.predecessor(1,K).toString(K));
        sb.setCharAt(0,'G');
        Assert.assertEquals(sb.toString(),kkk.predecessor(2,K).toString(K));
        sb.setCharAt(0,'T');
        Assert.assertEquals(sb.toString(),kkk.predecessor(3,K).toString(K));
    }

    @Test
    public void testReverseComplementation()
    {
        Assert.assertEquals(StringKmerizer.toKmer("ACGTACGA").rc(8),StringKmerizer.toKmer("TCGTACGT"));
        Assert.assertEquals(StringKmerizer.toKmer("ACGTACGTC").rc(9),StringKmerizer.toKmer("GACGTACGT"));
    }

    @Test
    public void testCanonicalization()
    {
        Assert.assertEquals(StringKmerizer.toKmer("ACGTACGTC").canonical(9),StringKmerizer.toKmer("ACGTACGTC"));
        Assert.assertEquals(StringKmerizer.toKmer("GACGTACGT").canonical(9),StringKmerizer.toKmer("ACGTACGTC"));
    }

    @Test
    public void testComparison()
    {
        Kmer kkk1 = StringKmerizer.toKmer("ACGTA");
        Kmer kkk2 = StringKmerizer.toKmer("ACGTC");
        Assert.assertTrue(kkk1.compareTo(kkk1) == 0);
        Assert.assertTrue(kkk2.compareTo(kkk2) == 0);
        Assert.assertTrue(kkk1.compareTo(kkk2) < 0);
        Assert.assertTrue(kkk2.compareTo(kkk1) > 0);
    }

    @Test
    public void testHashCode()
    {
        Assert.assertNotEquals(StringKmerizer.toKmer("TAGCGTA").hashCode(),StringKmerizer.toKmer("TAGCGTC").hashCode());
        Assert.assertEquals(StringKmerizer.toKmer("TAGGGTC").hashCode(),StringKmerizer.toKmer("TAGGGTC").hashCode());
    }

    @Test
    public void testKmerizationAcrossN()
    {
        StringKmerizer kmerizer = new StringKmerizer("AAAAANTTTTT",5);
        Assert.assertTrue(kmerizer.hasNext());
        Assert.assertEquals(kmerizer.next(),StringKmerizer.toKmer("AAAAA"));
        Assert.assertTrue(kmerizer.hasNext());
        Assert.assertEquals(kmerizer.next(),StringKmerizer.toKmer("TTTTT"));
        Assert.assertTrue(!kmerizer.hasNext());
    }
}
