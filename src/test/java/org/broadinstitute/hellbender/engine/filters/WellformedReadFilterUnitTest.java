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

package org.broadinstitute.hellbender.engine.filters;


import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.utils.sam.ArtificialSAMUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.lang.reflect.Method;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

/**
 * Tests for the Wellformed read filter. .
 */
public class WellformedReadFilterUnitTest extends ReadFilterTest {

    @Test(enabled = true)
    public void testCheckSeqStored () {

        final SAMRecord goodRead = ArtificialSAMUtils.createArtificialRead(new byte[]{(byte) 'A'}, new byte[]{(byte) 'A'}, "1M");
        final SAMRecord badRead = ArtificialSAMUtils.createArtificialRead(new byte[]{}, new byte[]{}, "1M");
        badRead.setReadString("*");

        Assert.assertTrue(WellformedReads.seqIsStored(goodRead));
        Assert.assertFalse(WellformedReads.seqIsStored(badRead));
    }

    @Test(enabled = true, dataProvider= "UnsupportedCigarOperatorDataProvider")
    public void testCigarNOperatorFilter(String cigarString) {

        final ReadFilter filter = ReadFilterLibrary.WELLFORMED;
        final SAMRecord read = buildSAMRecord(cigarString);
        final boolean containsN = cigarString.contains("N");
        Assert.assertEquals(containsN, !filter.test(read), cigarString);
    }

    protected SAMRecord buildSAMRecord(final String cigarString) {
        final Cigar nContainingCigar = TextCigarCodec.decode(cigarString);
        return  this.createRead(nContainingCigar, 1, 0, 10);
    }

    /**
     * Cigar test data for unsupported operator test.
     * Each element of this array corresponds to a test case. In turn the first element of the test case array is the
     * Cigar string for that test case and the second indicates whether it should be filtered due to the presence of a
     * unsupported operator
     */
    private static final String[] TEST_CIGARS =  {
       "101M10D20I10M",
       "6M14N5M",
       "1N",
       "101M",
       "110N",
       "2N4M",
       "4M2N",
       "3M1I1M",
       "1M2I2M",
       "1M10N1I1M",
       "1M1I1D",
       "11N12M1I34M12N"
    };

    @DataProvider(name= "UnsupportedCigarOperatorDataProvider")
    public Iterator<Object[]> unsupportedOperatorDataProvider(final Method testMethod) {
        final List<Object[]> result = new LinkedList<>();
        for (final String cigarString : TEST_CIGARS) {
            result.add(new Object[] { cigarString });
        }
        return result.iterator();
    }
}
