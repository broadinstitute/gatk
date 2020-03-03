package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Stream;


public class SvCigarUtilsUnitTest extends GATKBaseTest {

    @Test(expectedExceptions=IllegalArgumentException.class, groups = "sv")
    public void testCigarChecker_emptyCigarElementList(){
        final List<CigarElement> emptyList = Collections.emptyList();
        SvCigarUtils.validateCigar(emptyList);
    }

    @Test(expectedExceptions=IllegalArgumentException.class, groups = "sv")
    public void testCigarChecker_deletionNeighboringClipping(){
        SvCigarUtils.validateCigar(TextCigarCodec.decode("10S10D10M").getCigarElements());
        SvCigarUtils.validateCigar(TextCigarCodec.decode("10H10S10D10M").getCigarElements());
        SvCigarUtils.validateCigar(TextCigarCodec.decode("10M10D10S").getCigarElements());
        SvCigarUtils.validateCigar(TextCigarCodec.decode("10M10D10S10H").getCigarElements());
    }

    @Test(expectedExceptions=IllegalArgumentException.class, groups = "sv")
    public void testCigarChecker_only1NonAlignment(){
        SvCigarUtils.validateCigar(TextCigarCodec.decode("10S").getCigarElements());
    }

    @Test(expectedExceptions=IllegalArgumentException.class, groups = "sv")
    public void testCigarChecker_noAlignment(){
        SvCigarUtils.validateCigar(TextCigarCodec.decode("10H10S10I10S10H").getCigarElements());
    }

    @Test(groups = "sv")
    public void testGetIndexOfFirstNonClippingBase(){
        Assert.assertEquals(SvCigarUtils.findIndexOfFirstNonClippingOperation(TextCigarCodec.decode("151M").getCigarElements(), true), 0);
        Assert.assertEquals(SvCigarUtils.findIndexOfFirstNonClippingOperation(TextCigarCodec.decode("151M").getCigarElements(), false), 0);
        Assert.assertEquals(SvCigarUtils.findIndexOfFirstNonClippingOperation(TextCigarCodec.decode("10S10D10M").getCigarElements(), true), 1);
        Assert.assertEquals(SvCigarUtils.findIndexOfFirstNonClippingOperation(TextCigarCodec.decode("10H10S10D10M").getCigarElements(), true), 2);
        Assert.assertEquals(SvCigarUtils.findIndexOfFirstNonClippingOperation(TextCigarCodec.decode("10M10D10S").getCigarElements(), false), 1);
        Assert.assertEquals(SvCigarUtils.findIndexOfFirstNonClippingOperation(TextCigarCodec.decode("10M10D10S10H").getCigarElements(), false), 1);
    }

    @DataProvider(name = "readWalkDistanceTestDataException")
    private Object[][] createReadWalkDistanceTestDataException() {
        final List<Object[]> data = new ArrayList<>(20);

        data.add(new Object[]{TextCigarCodec.decode("50M10N101M"), 41, 10, false, 0});
        data.add(new Object[]{TextCigarCodec.decode("50M10P101M"), 41, 10, false, 0});

        final Cigar cigar = TextCigarCodec.decode("35H40S10M20I25M30D50M55S60H");
        data.add(new Object[]{cigar, -1, 10, false, 0});
        data.add(new Object[]{cigar, 0, 10, false, 0});
        data.add(new Object[]{cigar, 41, -1, false, 0});
        data.add(new Object[]{cigar, 41, 0, false, 0});

        data.add(new Object[]{cigar, 1, 116, false, 0});
        data.add(new Object[]{cigar, 200, 116, false, 0});

        data.add(new Object[]{cigar, 96, 51, false, 0});
        data.add(new Object[]{cigar, 145, 2, false, 0});
        data.add(new Object[]{cigar, 146, 1, false, 0});

        data.add(new Object[]{cigar, 96, 67, true, 0});
        data.add(new Object[]{cigar, 40, 1, true, 0});
        data.add(new Object[]{cigar, 41, 2, true, 0});

        return data.toArray(new Object[data.size()][]);
    }

    @Test(dataProvider = "readWalkDistanceTestDataException", groups = "sv", expectedExceptions = IllegalArgumentException.class)
    public void testReadWalkDistanceTestDataException(final Cigar cigar, final int startInclusive, final int refWalkDist,
                                                      final boolean walkBackwards, final int expectedReadWalkDist) {
        Assert.assertEquals(SvCigarUtils.computeAssociatedDistOnRead(cigar, startInclusive, refWalkDist, walkBackwards), expectedReadWalkDist);
    }

    @DataProvider(name = "readWalkDistanceTestData")
    private Object[][] createReadWalkDistanceTestData() {

        final List<Object[]> data = new ArrayList<>(20);
        final Cigar cigar = TextCigarCodec.decode("35H40S10M20I25M30D50M55S60H");

        data.add(new Object[]{cigar, 1, 5, false, 45});
        data.add(new Object[]{cigar, 1, 10, false, 50});
        data.add(new Object[]{cigar, 1, 16, false, 76});
        data.add(new Object[]{cigar, 1, 64, false, 95});
        data.add(new Object[]{cigar, 1, 66, false, 96});
        data.add(new Object[]{cigar, 11, 5, false, 35});
        data.add(new Object[]{cigar, 41, 64, false, 55});

        data.add(new Object[]{cigar, 146, 1, true, 2});

        data.add(new Object[]{cigar, 181, 10, true, 46});
        data.add(new Object[]{cigar, 181, 50, true, 86});
        data.add(new Object[]{cigar, 181, 51, true, 86});
        data.add(new Object[]{cigar, 181, 80, true, 86});
        data.add(new Object[]{cigar, 181, 105, true, 111});
        data.add(new Object[]{cigar, 181, 106, true, 132});
        data.add(new Object[]{cigar, 181, 115, true, 141});

        return data.toArray(new Object[data.size()][]);
    }

    @Test(dataProvider = "readWalkDistanceTestData", groups = "sv")
    public void testReadWalkDistanceTestData(final Cigar cigar, final int startInclusiveOnRead, final int refWalkDist,
                                             final boolean walkBackwards, final int expectedReadWalkDist) {
        Assert.assertEquals(SvCigarUtils.computeAssociatedDistOnRead(cigar, startInclusiveOnRead, refWalkDist, walkBackwards), expectedReadWalkDist);
    }
}
