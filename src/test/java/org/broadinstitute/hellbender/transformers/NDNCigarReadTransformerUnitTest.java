/*
* Copyright (c) 2012 The Broad Institute
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.hellbender.transformers;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class NDNCigarReadTransformerUnitTest {

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
