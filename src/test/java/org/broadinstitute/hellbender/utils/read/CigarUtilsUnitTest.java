package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.utils.clipping.ReadClipperTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.lang.reflect.Method;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

public final class CigarUtilsUnitTest {

    @DataProvider(name = "testData_testCombineAdjacentCigarElements")
    public Iterator<Object[]> testData_testCombineAdjacentCigarElements(final Method testMethod) {
        final String[][] TEST_CIGARS = {
                {"10M", "10M"},
                {"10M10M", "20M"},
                {"10M10D", "10M10D"},
                {"2I2D", "4I"},//Ds and Is get combined
                {"2D2I", "4D"},//Ds and Is get combined
        };
        final List<Object[]> result = new LinkedList<>();
        Collections.addAll(result, TEST_CIGARS);
        return result.iterator();
    }

    @Test(dataProvider = "testData_testCombineAdjacentCigarElements")
    public void testCombineAdjacentCigarElements(final String cigarStrIn, final String expectedCigarStrOut){
        final Cigar cigarIn = TextCigarCodec.decode(cigarStrIn);
        final Cigar cigarOut = CigarUtils.combineAdjacentCigarElements(cigarIn);
        final String actualCigarStrOut = TextCigarCodec.encode(cigarOut);
        Assert.assertEquals(actualCigarStrOut, expectedCigarStrOut);
    }

    @DataProvider(name = "testData_ReadHasNonClippedBases")
    public Iterator<Object[]> testData_ReadHasNonClippedBases(final Method testMethod) {
        final Object[][] TEST_CIGARS = {
                {"10M", true},
                {"10M10M", true},
                {"10S10H", false},
        };
        final List<Object[]> result = new LinkedList<>();
        Collections.addAll(result, TEST_CIGARS);
        return result.iterator();
    }

    @Test(dataProvider = "testData_ReadHasNonClippedBases")
    public void testReadHasNonClippedBases(final String cigarStrIn, final boolean expected){
        final Cigar cigarIn = TextCigarCodec.decode(cigarStrIn);
        final boolean actual = CigarUtils.hasNonClippedBases(cigarIn);
        Assert.assertEquals(actual, expected);
    }

    @DataProvider(name = "testData_reclipCigar")
    public Iterator<Object[]> testData_reclipCigar(final Method testMethod) {
        final String[][] TEST_CIGARS = {
                {"10M", "10M", "10M"},
                {"10M", "1H1M1H", "1H10M1H"},
                {"10M", "1S1M1S", "1S10M1S"},
                {"10M", "1S1M1H", "1S10M1H"},
                {"10M", "1H1M1S", "1H10M1S"},
                {"10M", "1H1S1M1S1H", "1H1S10M1S1H"},
                {"10M", "1H1S1M1I1M1S1H", "1H1S10M1S1H"},
        };
        final List<Object[]> result = new LinkedList<>();
        Collections.addAll(result, TEST_CIGARS);
        return result.iterator();
    }

    @Test(dataProvider = "testData_reclipCigar")
    public void testReclipCigar(final String cigarStrIn1, final String cigarStrIn2, final String expectedCigarStrOut){
        final Cigar cigarIn = TextCigarCodec.decode(cigarStrIn1);
        final SAMRecord read = ReadClipperTestUtils.makeReadFromCigar(cigarStrIn2);
        final Cigar cigarOut = CigarUtils.reclipCigar(cigarIn, read);
        final String actualCigarStrOut = TextCigarCodec.encode(cigarOut);
        Assert.assertEquals(actualCigarStrOut, expectedCigarStrOut);
    }


    @DataProvider(name = "testData_invertCigar")
    public Iterator<Object[]> testData_invertCigar(final Method testMethod) {
        final String[][] TEST_CIGARS = {
                {"10M", "10M"},
                {"1M2M", "2M1M"},
                {"10M10D", "10D10M"},
                {"2I2D", "2D2I"},
                {"2D2I", "2I2D"},
        };
        final List<Object[]> result = new LinkedList<>();
        Collections.addAll(result, TEST_CIGARS);
        return result.iterator();
    }

    @Test(dataProvider = "testData_invertCigar")
    public void testInvertCigar(final String cigarStrIn, final String expectedCigarStrOut){
        final Cigar cigarIn = TextCigarCodec.decode(cigarStrIn);
        final Cigar cigarOut = CigarUtils.invertCigar(cigarIn);
        final String actualCigarStrOut = TextCigarCodec.encode(cigarOut);
        Assert.assertEquals(actualCigarStrOut, expectedCigarStrOut);
    }

    @DataProvider(name = "testData_unclipCigar")
    public Iterator<Object[]> testData_unclipCigar(final Method testMethod) {
        final String[][] TEST_CIGARS = {
                {"10M", "10M"},
                {"1M1D1M1I1M", "1M1D1M1I1M"},
                {"1H10M", "10M"},
                {"1H1S10M", "10M"},
                {"10M1S1H", "10M"},
                {"1H1S10M1S1H", "10M"},
                {"1H1S10M1D10M1S1H", "10M1D10M"},
        };
        final List<Object[]> result = new LinkedList<>();
        Collections.addAll(result, TEST_CIGARS);
        return result.iterator();
    }

    @Test(dataProvider = "testData_unclipCigar")
    public void testUnclipCigar(final String cigarStrIn, final String expectedCigarStrOut){
        final Cigar cigarIn = TextCigarCodec.decode(cigarStrIn);
        final Cigar cigarOut = CigarUtils.trimReadToUnclippedBases(cigarIn);
        final String actualCigarStrOut = TextCigarCodec.encode(cigarOut);
        Assert.assertEquals(actualCigarStrOut, expectedCigarStrOut);
    }

    @DataProvider(name = "testData_containsNOperator")
    public Iterator<Object[]> testData_containsNOperator(final Method testMethod) {
        final Object[][] TEST_CIGARS = {
                {"10M", false},
                {"10M1N", true},
                {"10M1N1M", true},
        };
        final List<Object[]> result = new LinkedList<>();
        Collections.addAll(result, TEST_CIGARS);
        return result.iterator();
    }

    @Test(dataProvider = "testData_containsNOperator")
    public void testContainsNOperator(final String cigarStrIn, final boolean expected){
        final Cigar cigarIn = TextCigarCodec.decode(cigarStrIn);
        final boolean actual = CigarUtils.containsNOperator(cigarIn);
        Assert.assertEquals(actual, expected, cigarStrIn);
    }

    @DataProvider(name = "testData_countRefBasesBasedOnCigar")
    public Iterator<Object[]> testData_countRefBasesBasedOnCigar(final Method testMethod) {
        final Object[][] TEST_CIGARS = {
                {"10M", 0, 1, 10},
                {"10M1D", 0, 1, 10},
                {"10M1D1S", 0, 1, 10},
                {"10M1D1N1S", 0, 1, 10},
                {"10M1D1N1S1H", 0, 1, 10},
                {"10M1D", 0, 2, 11},
                {"10M1D1S", 0, 2, 11},
                {"10M1D1N1S", 0, 2, 11},
                {"10M1D1N1S1H", 0, 2, 11},
                {"10M1=1X1S", 0, 2, 11},
                {"10M1=1X1S1H", 0, 2, 11},
                {"10M1D2N4S", 1, 3, 3},
                {"10M1D2N4S8H", 1, 4, 7},
                {"10M1I2N4S", 1, 3, 2},
                {"10M1I2N4S8H", 1, 4, 6},

                //https://github.com/broadinstitute/hellbender/issues/450
                {"1M1S1H1H", 1, 2, 1},
                {"1M1I1H1H", 1, 2, 0},//I is skipped
                {"1M1P1H1H", 1, 2, 0},//P is skipped
                {"1M1H1H1H", 1, 2, 1},
        };
        final List<Object[]> result = new LinkedList<>();
        Collections.addAll(result, TEST_CIGARS);
        return result.iterator();
    }

    @Test(dataProvider = "testData_countRefBasesBasedOnCigar")
    public void testCountRefBasesBasedOnCigar(final String cigarStrIn, final int start, final int end, final int expected){
        final SAMRecord read = ReadClipperTestUtils.makeReadFromCigar(cigarStrIn);
        final int actual = CigarUtils.countRefBasesBasedOnCigar(read, start, end);
        Assert.assertEquals(actual, expected, cigarStrIn);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testCountRefBasesBasedOnCigarNull(){
        CigarUtils.countRefBasesBasedOnCigar(null, 1, 2);
    }
    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testCountRefBasesBasedOnCigarStart1(){
        final String cigarStrIn = "1M1=1X";
        final SAMRecord read = ReadClipperTestUtils.makeReadFromCigar(cigarStrIn);
        CigarUtils.countRefBasesBasedOnCigar(read, -1, 1);
    }
    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testCountRefBasesBasedOnCigarStart2(){
        final String cigarStrIn = "1M1=1X";
        final SAMRecord read = ReadClipperTestUtils.makeReadFromCigar(cigarStrIn);
        CigarUtils.countRefBasesBasedOnCigar(read, 2, 1);
    }
    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testCountRefBasesBasedOnCigarEnd2(){
        final String cigarStrIn = "1M1=1X";
        final SAMRecord read = ReadClipperTestUtils.makeReadFromCigar(cigarStrIn);
        CigarUtils.countRefBasesBasedOnCigar(read, 1, 6);
    }

}
