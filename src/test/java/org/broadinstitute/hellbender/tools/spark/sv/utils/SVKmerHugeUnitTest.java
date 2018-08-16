package org.broadinstitute.hellbender.tools.spark.sv.utils;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.Utils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection.KMER_SIZE;
import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection.MAX_DUST_SCORE;

/**
 * Unit tests for SVKmerHuge and SVKmerizer.
 */
public class SVKmerHugeUnitTest extends GATKBaseTest {
    @Test(groups = "sv")
    public void testDefaultConstruction() {
        Assert.assertEquals(new SVKmerHuge(10).toString(10), "AAAAAAAAAA");
        Assert.assertEquals(new SVKmerHuge(11).toString(11), "AAAAAAAAAAA");
    }

    @Test(expectedExceptions = IllegalArgumentException.class, groups = "sv")
    public void testDefaultConstructionWithTooSmallK() {
        final SVKmerHuge tooSmallK = new SVKmerHuge(0);
    }

    @DataProvider(name = "sequenceStrings")
    public Object[][] getSequenceStrings() {
        return new Object[][] {
                {"ACGTACGTACGTACGTACGTACGTACGTAC"},
                {"ACGTACGTACGTACGTACGTACGTACGTACG"},
                {"ACGTACGTACGTACGTACGTACGTACGTACGT"},
                {"ACGTACGTACGTACGTACGTACGTACGTACGTA"},
                {"ACGTACGTACGTACGTACGTACGTACGTACGTAC"},
                {"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"}
        };
    }

    @Test(dataProvider = "sequenceStrings", groups = "sv")
    public void testConstructionAndToString( final String str ) {
        Assert.assertEquals(toKmer(str).toString(str.length()), str);
        Assert.assertEquals(SVKmerizer.toKmer(str.getBytes(),new SVKmerHuge(str.length())).toString(str.length()), str);
    }

    @Test(dataProvider = "sequenceStrings", groups = "sv")
    public void testSuccessor( final String str ) {
        final int K = str.length();
        final StringBuilder sb = new StringBuilder(str);
        sb.delete(0, 1);
        sb.append('A');
        final SVKmer kkk = toKmer(str);
        Assert.assertEquals(kkk.successor(SVKmer.Base.A, K).toString(K), sb.toString());
        sb.setCharAt(K-1, 'C');
        Assert.assertEquals(kkk.successor(SVKmer.Base.C, K).toString(K), sb.toString());
        sb.setCharAt(K-1, 'G');
        Assert.assertEquals(kkk.successor(SVKmer.Base.G, K).toString(K), sb.toString());
        sb.setCharAt(K-1, 'T');
        Assert.assertEquals(kkk.successor(SVKmer.Base.T, K).toString(K), sb.toString());
    }

    @Test(dataProvider = "sequenceStrings", groups = "sv")
    public void testPredecessor( final String str ) {
        final int K = str.length();
        final StringBuilder sb = new StringBuilder(str);
        sb.insert(0, 'A');
        sb.setLength(K);
        final SVKmer kkk = toKmer(str);
        Assert.assertEquals(kkk.predecessor(SVKmer.Base.A, K).toString(K), sb.toString());
        sb.setCharAt(0, 'C');
        Assert.assertEquals(kkk.predecessor(SVKmer.Base.C, K).toString(K), sb.toString());
        sb.setCharAt(0, 'G');
        Assert.assertEquals(kkk.predecessor(SVKmer.Base.G, K).toString(K), sb.toString());
        sb.setCharAt(0, 'T');
        Assert.assertEquals(kkk.predecessor(SVKmer.Base.T, K).toString(K), sb.toString());
    }

    @Test(dataProvider = "sequenceStrings", groups = "sv")
    public void testFirstBase( final String str ) {
        Assert.assertEquals(toKmer(str).firstBase(str.length()),
                            SVKmer.Base.valueOf(str.substring(0,1)));
    }

    @Test(dataProvider = "sequenceStrings", groups = "sv")
    public void testFirstTrimer( final String str ) {
        final int testVal = (int)((SVKmer.Base.valueOf(str.substring(0,1)).value << 4) |
                                  (SVKmer.Base.valueOf(str.substring(1,2)).value << 2) |
                                   SVKmer.Base.valueOf(str.substring(2,3)).value);
        Assert.assertEquals(toKmer(str).firstTrimer(str.length()), testVal);
    }

    @Test(dataProvider = "sequenceStrings", groups = "sv")
    public void testLastBase( final String str ) {
        Assert.assertEquals(toKmer(str).lastBase(),
                            SVKmer.Base.valueOf(str.substring(str.length()-1)));
    }

    @Test(dataProvider = "sequenceStrings", groups = "sv")
    public void testLastTrimer( final String str ) {
        final int K = str.length();
        int testVal = (int)((SVKmer.Base.valueOf(str.substring(K-3, K-2)).value << 4) |
                            (SVKmer.Base.valueOf(str.substring(K-2, K-1)).value << 2) |
                             SVKmer.Base.valueOf(str.substring(K-1)).value);
        Assert.assertEquals(toKmer(str).lastTrimer(), testVal);
    }

    @Test(dataProvider = "sequenceStrings", groups = "sv")
    public void testReverseComplementation( final String str ) {
        final int K = str.length();
        Assert.assertEquals(toKmer(str).reverseComplement(K).toString(K), rcString(str));
    }

    @Test(dataProvider = "sequenceStrings", groups = "sv")
    public void testCanonicalization( final String str ) {
        final int K = str.length();
        if ( (K & 1) == 0 ) return;
        final boolean isCanonical = str.charAt(K/2) == 'A' || str.charAt(K/2) == 'C';
        Assert.assertEquals(toKmer(str).canonical(K).toString(K), isCanonical ? str : rcString(str));
    }

    @Test(groups = "sv")
    public void testComparison() {
        final SVKmerHuge kkk1 = toKmer("ACGTA");
        final SVKmerHuge kkk2 = toKmer("ACGTC");
        Assert.assertTrue(kkk1.compareTo(kkk1) == 0);
        Assert.assertTrue(kkk2.compareTo(kkk2) == 0);
        Assert.assertTrue(kkk1.compareTo(kkk2) < 0);
        Assert.assertTrue(kkk2.compareTo(kkk1) > 0);
    }

    private SVKmerHuge toKmer( final String str ) {
        return SVKmerizer.toKmer(str, new SVKmerHuge(str.length()));
    }

    private String rcString( final String str ) {
        final int K = str.length();
        final StringBuilder sb = new StringBuilder(K);
        for ( int idx = 0; idx != K; ++idx ) {
            sb.append(SVKmer.Base.valueOf(str.substring(idx, idx + 1)).complement().name());
        }
        return sb.reverse().toString();
    }
}
