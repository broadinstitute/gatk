package org.broadinstitute.hellbender.transformers;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public final class NDNCigarReadTransformerUnitTest {

    @DataProvider(name = "filteringIteratorTestData")
    public String[][] getFilteringIteratorTestData() {
        return new String[][] {
                {"1M1N1N1M","1M1N1N1M"},           // NN elements
                {"1M1N1D4M","1M1N1D4M"},           // ND
                {"1M1N3M","1M1N3M"},               // N
                {"1M1N2I1N3M","1M1N2I1N3M"},       // NIN
                {"1M1N3D2N1M","1M6N1M"},
                {"1M2N2D2N1M1D3N1D1N1M2H","1M6N1M1D5N1M2H"},
                {"1H2S1M1N3D2N1M","1H2S1M6N1M"},
                {"10M628N2D203N90M","10M833N90M"}
        };
    }

    @Test(dataProvider = "filteringIteratorTestData")
    public void testCigarRefactoring (final String originalCigarString, final String expectedString) {
        Cigar originalCigar = TextCigarCodec.decode(originalCigarString);
        String actualString = NDNCigarReadTransformer.refactorNDNtoN(originalCigar).toString();
        Assert.assertEquals(actualString, expectedString, "cigar string " + originalCigarString + " should become: " + expectedString + " but got: " + actualString);
    }

}
