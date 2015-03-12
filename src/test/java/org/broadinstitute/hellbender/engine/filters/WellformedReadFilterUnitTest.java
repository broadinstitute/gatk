package org.broadinstitute.hellbender.engine.filters;


import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.hellbender.utils.sam.ReadUtils;
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

    @Test
    public void testCheckSeqStored () {

        final SAMRecord goodRead = ArtificialSAMUtils.createArtificialRead(new byte[]{(byte) 'A'}, new byte[]{(byte) 'A'}, "1M");
        final SAMRecord badRead = ArtificialSAMUtils.createArtificialRead(new byte[]{}, new byte[]{}, "1M");
        badRead.setReadString("*");

        Assert.assertTrue(ReadUtils.seqIsStored(goodRead));
        Assert.assertFalse(ReadUtils.seqIsStored(badRead));
    }

    @Test(dataProvider= "UnsupportedCigarOperatorDataProvider")
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
