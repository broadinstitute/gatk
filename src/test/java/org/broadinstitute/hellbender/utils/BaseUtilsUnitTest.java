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


package org.broadinstitute.hellbender.utils;

import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Random;


public class BaseUtilsUnitTest extends BaseTest {
    @BeforeClass
    public void init() { }

    @Test
    public void testMostFrequentBaseFraction() {
        compareFrequentBaseFractionToExpected("AAAAA", 1.0);
        compareFrequentBaseFractionToExpected("ACCG", 0.5);
        compareFrequentBaseFractionToExpected("ACCCCTTTTG", 4.0/10.0);
    }

    private void compareFrequentBaseFractionToExpected(String sequence, double expected) {
        double fraction = BaseUtils.mostFrequentBaseFraction(sequence.getBytes());
        Assert.assertTrue(MathUtils.compareDoubles(fraction, expected) == 0);
    }

    @Test
    public void testConvertIUPACtoN() {

        checkBytesAreEqual(BaseUtils.convertIUPACtoN(new byte[]{'A', 'A', 'A'}, false, false), new byte[]{'A', 'A', 'A'});
        checkBytesAreEqual(BaseUtils.convertIUPACtoN(new byte[]{'W', 'A', 'A'}, false, false), new byte[]{'N', 'A', 'A'});
        checkBytesAreEqual(BaseUtils.convertIUPACtoN(new byte[]{'A', 'M', 'A'}, false, false), new byte[]{'A', 'N', 'A'});
        checkBytesAreEqual(BaseUtils.convertIUPACtoN(new byte[]{'A', 'A', 'K'}, false, false), new byte[]{'A', 'A', 'N'});
        checkBytesAreEqual(BaseUtils.convertIUPACtoN(new byte[]{'M', 'M', 'M'}, false, false), new byte[]{'N', 'N', 'N'});
    }

    private void checkBytesAreEqual(final byte[] b1, final byte[] b2) {
        for ( int i = 0; i < b1.length; i++ )
            Assert.assertEquals(b1[i], b2[i]);
    }

    /**
     * Converts a IUPAC nucleotide code to a pair of bases
     *
     * @param code
     * @return 0, 1, 2, 3, or -1 if the base can't be understood
     */
    @Deprecated
    static public char[] iupacToBases(char code) {
        char[] bases = new char[2];
        switch (code) {
            case '*':               // the wildcard character counts as an A
            case 'A':
            case 'a':
                bases[0] = bases[1] = 'A';
                break;
            case 'C':
            case 'c':
                bases[0] = bases[1] = 'C';
                break;
            case 'G':
            case 'g':
                bases[0] = bases[1] = 'G';
                break;
            case 'T':
            case 't':
                bases[0] = bases[1] = 'T';
                break;
            case 'R':
            case 'r':
                bases[0] = 'A';
                bases[1] = 'G';
                break;
            case 'Y':
            case 'y':
                bases[0] = 'C';
                bases[1] = 'T';
                break;
            case 'S':
            case 's':
                bases[0] = 'G';
                bases[1] = 'C';
                break;
            case 'W':
            case 'w':
                bases[0] = 'A';
                bases[1] = 'T';
                break;
            case 'K':
            case 'k':
                bases[0] = 'G';
                bases[1] = 'T';
                break;
            case 'M':
            case 'm':
                bases[0] = 'A';
                bases[1] = 'C';
                break;
            default:
                bases[0] = bases[1] = 'N';
        }
        return bases;
    }

    @Test
    public void testConvertBasesToIUPAC() {

        for ( final BaseUtils.Base b : BaseUtils.Base.values() ) {
            if ( BaseUtils.isRegularBase(b.base) )
                Assert.assertEquals(BaseUtils.basesToIUPAC(b.base, b.base), b.base, "testing same base");
        }

        Assert.assertEquals(BaseUtils.basesToIUPAC((byte) 'A', (byte) 'X'), 'N', "testing non-standard base");
        Assert.assertEquals(BaseUtils.basesToIUPAC((byte)'X', (byte)'A'), 'N', "testing non-standard base");
        Assert.assertEquals(BaseUtils.basesToIUPAC((byte)'X', (byte)'X'), 'N', "testing non-standard base");

        Assert.assertEquals(BaseUtils.basesToIUPAC((byte)'A', (byte)'T'), 'W', "testing A/T=W");
        Assert.assertEquals(BaseUtils.basesToIUPAC((byte)'T', (byte)'A'), 'W', "testing T/A=W");
        Assert.assertEquals(BaseUtils.basesToIUPAC((byte) 'G', (byte) 'T'), 'K', "testing G/T=K");
        Assert.assertEquals(BaseUtils.basesToIUPAC((byte) 'T', (byte) 'G'), 'K', "testing T/G=K");
    }

    @Test
    public void testTransitionTransversion() {
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'A', (byte)'T' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'A', (byte)'C' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'A', (byte)'G' ) == BaseUtils.BaseSubstitutionType.TRANSITION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'C', (byte)'A' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'C', (byte)'T' ) == BaseUtils.BaseSubstitutionType.TRANSITION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'C', (byte)'G' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'T', (byte)'A' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'T', (byte)'C' ) == BaseUtils.BaseSubstitutionType.TRANSITION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'T', (byte)'G' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'G', (byte)'A' ) == BaseUtils.BaseSubstitutionType.TRANSITION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'G', (byte)'T' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'G', (byte)'C' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );

        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'a', (byte)'T' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'a', (byte)'C' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'A', (byte)'T' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'A', (byte)'C' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'A', (byte)'t' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'A', (byte)'c' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'a', (byte)'t' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'a', (byte)'c' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
    }

    @Test
    public void testReverseComplementString() {
        compareRCStringToExpected("ACGGT", "ACCGT");
        compareRCStringToExpected("TCGTATATCTCGCTATATATATATAGCTCTAGTATA", "TATACTAGAGCTATATATATATAGCGAGATATACGA");
        compareRCStringToExpected("AAAN", "NTTT");
    }

    private void compareRCStringToExpected(String fw, String rcExp) {
        String rcObs = BaseUtils.simpleReverseComplement(fw);

        Assert.assertTrue(rcObs.equals(rcExp));
    }

    @Test(dataProvider="baseComparatorData")
    public void testBaseComparator(final Collection<byte[]> basesToSort) {
        final ArrayList<byte[]> sorted = new ArrayList<>(basesToSort);
        Collections.sort(sorted, BaseUtils.BASES_COMPARATOR);
        for (int i = 0; i < sorted.size(); i++)   {
            Assert.assertEquals(BaseUtils.BASES_COMPARATOR.compare(sorted.get(i),sorted.get(i)),0);
            final String iString = new String(sorted.get(i));
            for (int j = i; j < sorted.size(); j++) {
                final String jString = new String(sorted.get(j));
                if (iString.compareTo(jString) == 0)
                    Assert.assertEquals(BaseUtils.BASES_COMPARATOR.compare(sorted.get(i),sorted.get(j)),0);
                else
                    Assert.assertTrue(BaseUtils.BASES_COMPARATOR.compare(sorted.get(i),sorted.get(j)) * iString.compareTo(jString) > 0);
                Assert.assertTrue(BaseUtils.BASES_COMPARATOR.compare(sorted.get(i),sorted.get(j)) <= 0);
            }
        }
    }

    @DataProvider(name="baseComparatorData")
    public Object[][] baseComparatorData() {
        final int testCount = 10;
        final int testSizeAverage = 10;
        final int testSizeDeviation = 10;
        final int haplotypeSizeAverage = 100;
        final int haplotypeSizeDeviation = 100;

        final Object[][] result = new Object[testCount][];

        Utils.resetRandomGenerator();
        final Random rnd = Utils.getRandomGenerator();

        for (int i = 0; i < testCount; i++) {
            final int size = (int) Math.max(0,rnd.nextDouble() * testSizeDeviation + testSizeAverage);
            final ArrayList<byte[]> bases = new ArrayList<>(size);
            for (int j = 0; j < size; j++) {
                final int jSize = (int) Math.max(0,rnd.nextDouble() * haplotypeSizeDeviation + haplotypeSizeAverage);
                final byte[] b = new byte[jSize];
                for (int k = 0; k < jSize; k++)
                    b[k] = BaseUtils.baseIndexToSimpleBase(rnd.nextInt(4));
                bases.add(b);
            }
            result[i] = new Object[] { bases };
        }
        return result;
    }
}
