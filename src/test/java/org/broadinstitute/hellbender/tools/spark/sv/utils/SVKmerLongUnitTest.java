package org.broadinstitute.hellbender.tools.spark.sv.utils;

import org.broadinstitute.hellbender.tools.spark.sv.SVConstants;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Unit tests for SVKmerLong and SVKmerizer.
 */
public class SVKmerLongUnitTest extends BaseTest {
    @Test
    public void testDefaultConstruction() {
        Assert.assertEquals(new SVKmerLong(10).toString(10), "AAAAAAAAAA");
        Assert.assertEquals(new SVKmerLong(11).toString(11), "AAAAAAAAAAA");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testDefaultConstructionWithTooLargeK() {
        final SVKmerLong tooBigK = new SVKmerLong(64);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testDefaultConstructionWithTooSmallK() {
        final SVKmerLong tooSmallK = new SVKmerLong(0);
    }

    @DataProvider(name = "sequenceStrings")
    public Object[][] getSequenceStrings() {
        return new Object[][] {
                {"ACGTACGTACGT"}, {"ACGTACGTACGTC"}, {"ACGTACGTACGTCC"}, {"ACGTACGTACGTCCC"}, {"ACGTACGTACGTCCCC"}
        };
    }

    @Test(dataProvider = "sequenceStrings")
    public void testConstructionAndToString( final String str ) {
        Assert.assertEquals(str, SVKmerizer.toKmer(str,new SVKmerLong(str.length())).toString(str.length()));
        Assert.assertEquals(str, SVKmerizer.toKmer(str.getBytes(),new SVKmerLong(str.length())).toString(str.length()));
    }

    @Test(dataProvider = "sequenceStrings")
    public void testSuccessor( final String str ) {
        final int K = str.length();
        final StringBuilder sb = new StringBuilder(str);
        sb.delete(0, 1);
        sb.append('A');
        final SVKmer kkk = SVKmerizer.toKmer(str,new SVKmerLong(str.length()));
        Assert.assertEquals(sb.toString(), kkk.successor(SVKmerLong.Base.A, K).toString(K));
        sb.setCharAt(K-1, 'C');
        Assert.assertEquals(sb.toString(), kkk.successor(SVKmerLong.Base.C, K).toString(K));
        sb.setCharAt(K-1, 'G');
        Assert.assertEquals(sb.toString(), kkk.successor(SVKmerLong.Base.G, K).toString(K));
        sb.setCharAt(K-1, 'T');
        Assert.assertEquals(sb.toString(), kkk.successor(SVKmerLong.Base.T, K).toString(K));
    }

    @Test(dataProvider = "sequenceStrings")
    public void testPredecessor( final String str ) {
        final int K = str.length();
        final StringBuilder sb = new StringBuilder(str);
        sb.insert(0, 'A');
        sb.setLength(K);
        final SVKmer kkk = SVKmerizer.toKmer(str,new SVKmerLong(str.length()));
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
        Assert.assertEquals(SVKmerizer.toKmer("ACGTACGA",new SVKmerLong(8)).reverseComplement(8), SVKmerizer.toKmer("TCGTACGT",new SVKmerLong(8)));
        Assert.assertEquals(SVKmerizer.toKmer("ACGTACGTC",new SVKmerLong(9)).reverseComplement(9), SVKmerizer.toKmer("GACGTACGT",new SVKmerLong(9)));
        Assert.assertEquals(SVKmerizer.toKmer("ACGTCCGTC",new SVKmerLong(9)).reverseComplement(9), SVKmerizer.toKmer("GACGGACGT",new SVKmerLong(9)));
        Assert.assertEquals(SVKmerizer.toKmer("ACGTGCGTC",new SVKmerLong(9)).reverseComplement(9), SVKmerizer.toKmer("GACGCACGT",new SVKmerLong(9)));
        Assert.assertEquals(SVKmerizer.toKmer("ACGTTCGTC",new SVKmerLong(9)).reverseComplement(9), SVKmerizer.toKmer("GACGAACGT",new SVKmerLong(9)));
    }

    @Test
    public void testCanonicalization() {
        Assert.assertEquals(SVKmerizer.toKmer("ACGTACGTC",new SVKmerLong(9)).canonical(9), SVKmerizer.toKmer("ACGTACGTC",new SVKmerLong(9)));
        Assert.assertEquals(SVKmerizer.toKmer("GACGTACGT",new SVKmerLong(9)).canonical(9), SVKmerizer.toKmer("ACGTACGTC",new SVKmerLong(9)));
    }

    @Test
    public void testComparison() {
        final SVKmerLong kkk1 = (SVKmerLong)SVKmerizer.toKmer("ACGTA",new SVKmerLong(5));
        final SVKmerLong kkk2 = (SVKmerLong)SVKmerizer.toKmer("ACGTC",new SVKmerLong(5));
        Assert.assertTrue(kkk1.compareTo(kkk1) == 0);
        Assert.assertTrue(kkk2.compareTo(kkk2) == 0);
        Assert.assertTrue(kkk1.compareTo(kkk2) < 0);
        Assert.assertTrue(kkk2.compareTo(kkk1) > 0);
    }

    @Test
    public void testHashCode() {
        Assert.assertNotEquals(SVKmerizer.toKmer("TAGCGTA",new SVKmerLong(7)).hashCode(), SVKmerizer.toKmer("TAGCGTC",new SVKmerLong(7)).hashCode());
        Assert.assertEquals(SVKmerizer.toKmer("TAGGGTC",new SVKmerLong(7)).hashCode(), SVKmerizer.toKmer("TAGGGTC",new SVKmerLong(7)).hashCode());
    }

    @Test
    public void testKmerization() {
        final SVKmerizer kmerizer = new SVKmerizer("AAAAATT", 5, 1, new SVKmerLong(7));
        Assert.assertTrue(kmerizer.hasNext());
        Assert.assertEquals(kmerizer.next(), SVKmerizer.toKmer("AAAAA",new SVKmerLong(5)));
        Assert.assertTrue(kmerizer.hasNext());
        Assert.assertEquals(kmerizer.next(), SVKmerizer.toKmer("AAAAT",new SVKmerLong(5)));
        Assert.assertTrue(kmerizer.hasNext());
        Assert.assertEquals(kmerizer.next(), SVKmerizer.toKmer("AAATT",new SVKmerLong(5)));
        Assert.assertTrue(!kmerizer.hasNext());
    }

    @Test
    public void testKmerizationAcrossN() {
        final SVKmerizer kmerizer = new SVKmerizer("AAAAANTTTTT", 5, 1, new SVKmerLong(11));
        Assert.assertTrue(kmerizer.hasNext());
        Assert.assertEquals(kmerizer.next(), SVKmerizer.toKmer("AAAAA",new SVKmerLong(5)));
        Assert.assertTrue(kmerizer.hasNext());
        Assert.assertEquals(kmerizer.next(), SVKmerizer.toKmer("TTTTT",new SVKmerLong(5)));
        Assert.assertTrue(!kmerizer.hasNext());
    }


    @Test
    public void testLowComplexityFilteredKmerization() {
        final byte[] seq =
                "CATAAAGCCTAAATAGCCCACACGTTCCCCTTAAATAAGACTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGATCACGATGGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCAT".getBytes();
        final SVKmer kmer = new SVKmerLong();
        final List<SVKmer> kmers =
                SVDUSTFilteredKmerizer.stream(seq, SVConstants.KMER_SIZE, SVConstants.MAX_DUST_SCORE, kmer)
                        .collect(SVUtils.arrayListCollector(seq.length-SVConstants.KMER_SIZE+1));
        final List<SVKmer> expectedKmers = Stream.of(
                "CATAAAGCCTAAATAGCCCACACGTTCCCCTTAAATAAGACTTAGGGTTAG",
                "ATAAAGCCTAAATAGCCCACACGTTCCCCTTAAATAAGACTTAGGGTTAGG",
                "TAAAGCCTAAATAGCCCACACGTTCCCCTTAAATAAGACTTAGGGTTAGGG",
                "AAAGCCTAAATAGCCCACACGTTCCCCTTAAATAAGACTTAGGGTTAGGGT",
                "AAGCCTAAATAGCCCACACGTTCCCCTTAAATAAGACTTAGGGTTAGGGTT",
                "AGCCTAAATAGCCCACACGTTCCCCTTAAATAAGACTTAGGGTTAGGGTTA",
                "GCCTAAATAGCCCACACGTTCCCCTTAAATAAGACTTAGGGTTAGGGTTAG",
                "CCTAAATAGCCCACACGTTCCCCTTAAATAAGACTTAGGGTTAGGGTTAGG",
                "CTAAATAGCCCACACGTTCCCCTTAAATAAGACTTAGGGTTAGGGTTAGGG",
                "TAAATAGCCCACACGTTCCCCTTAAATAAGACTTAGGGTTAGGGTTAGGGT",
                "AAATAGCCCACACGTTCCCCTTAAATAAGACTTAGGGTTAGGGTTAGGGTT",
                "AATAGCCCACACGTTCCCCTTAAATAAGACTTAGGGTTAGGGTTAGGGTTA",
                "ATAGCCCACACGTTCCCCTTAAATAAGACTTAGGGTTAGGGTTAGGGTTAG",
                "TAGCCCACACGTTCCCCTTAAATAAGACTTAGGGTTAGGGTTAGGGTTAGG",
                "AGCCCACACGTTCCCCTTAAATAAGACTTAGGGTTAGGGTTAGGGTTAGGG",
                "GCCCACACGTTCCCCTTAAATAAGACTTAGGGTTAGGGTTAGGGTTAGGGT",
                "TAGGGTTAGGGTTAGGGTTAGGGATCACGATGGATCACAGGTCTATCACCC",
                "AGGGTTAGGGTTAGGGTTAGGGATCACGATGGATCACAGGTCTATCACCCT",
                "GGGTTAGGGTTAGGGTTAGGGATCACGATGGATCACAGGTCTATCACCCTA",
                "GGTTAGGGTTAGGGTTAGGGATCACGATGGATCACAGGTCTATCACCCTAT",
                "GTTAGGGTTAGGGTTAGGGATCACGATGGATCACAGGTCTATCACCCTATT",
                "TTAGGGTTAGGGTTAGGGATCACGATGGATCACAGGTCTATCACCCTATTA",
                "TAGGGTTAGGGTTAGGGATCACGATGGATCACAGGTCTATCACCCTATTAA",
                "AGGGTTAGGGTTAGGGATCACGATGGATCACAGGTCTATCACCCTATTAAC",
                "GGGTTAGGGTTAGGGATCACGATGGATCACAGGTCTATCACCCTATTAACC",
                "GGTTAGGGTTAGGGATCACGATGGATCACAGGTCTATCACCCTATTAACCA",
                "GTTAGGGTTAGGGATCACGATGGATCACAGGTCTATCACCCTATTAACCAC",
                "TTAGGGTTAGGGATCACGATGGATCACAGGTCTATCACCCTATTAACCACT",
                "TAGGGTTAGGGATCACGATGGATCACAGGTCTATCACCCTATTAACCACTC",
                "AGGGTTAGGGATCACGATGGATCACAGGTCTATCACCCTATTAACCACTCA",
                "GGGTTAGGGATCACGATGGATCACAGGTCTATCACCCTATTAACCACTCAC",
                "GGTTAGGGATCACGATGGATCACAGGTCTATCACCCTATTAACCACTCACG",
                "GTTAGGGATCACGATGGATCACAGGTCTATCACCCTATTAACCACTCACGG",
                "TTAGGGATCACGATGGATCACAGGTCTATCACCCTATTAACCACTCACGGG",
                "TAGGGATCACGATGGATCACAGGTCTATCACCCTATTAACCACTCACGGGA",
                "AGGGATCACGATGGATCACAGGTCTATCACCCTATTAACCACTCACGGGAG",
                "GGGATCACGATGGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGC",
                "GGATCACGATGGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCT",
                "GATCACGATGGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTC",
                "ATCACGATGGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCT",
                "TCACGATGGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTC",
                "CACGATGGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCC",
                "ACGATGGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCA",
                "CGATGGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCAT",
                "GATGGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATG",
                "ATGGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGC",
                "TGGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCA",
                "GGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCAT"
        ).map(str -> SVKmerizer.toKmer(str,kmer)).collect(Collectors.toList());
        Assert.assertEquals(kmers,expectedKmers);
    }
}
