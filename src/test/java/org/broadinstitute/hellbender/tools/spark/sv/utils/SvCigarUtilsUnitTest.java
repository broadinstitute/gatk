package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SvCigarUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Collections;
import java.util.List;
import java.util.stream.Stream;


public class SvCigarUtilsUnitTest {

    @Test(groups = "sv")
    public void testClippingArithmetic() {
        Cigar cigar = TextCigarCodec.decode("100M51S");
        Assert.assertEquals(SvCigarUtils.getTotalHardClipping(cigar), 0);
        Assert.assertEquals(SvCigarUtils.getNumClippedBases(true, cigar), 0);
        Assert.assertEquals(SvCigarUtils.getNumClippedBases(false, cigar), 51);

        cigar = TextCigarCodec.decode("51S100M");
        Assert.assertEquals(SvCigarUtils.getTotalHardClipping(cigar), 0);
        Assert.assertEquals(SvCigarUtils.getNumClippedBases(true, cigar), 51);
        Assert.assertEquals(SvCigarUtils.getNumClippedBases(false, cigar), 0);

        cigar = TextCigarCodec.decode("100M51H");
        Assert.assertEquals(SvCigarUtils.getTotalHardClipping(cigar), 51);
        Assert.assertEquals(SvCigarUtils.getNumClippedBases(true, cigar), 0);
        Assert.assertEquals(SvCigarUtils.getNumClippedBases(false, cigar), 51);

        cigar = TextCigarCodec.decode("51H100M");
        Assert.assertEquals(SvCigarUtils.getTotalHardClipping(cigar), 51);
        Assert.assertEquals(SvCigarUtils.getNumClippedBases(true, cigar), 51);
        Assert.assertEquals(SvCigarUtils.getNumClippedBases(false, cigar), 0);

        cigar = TextCigarCodec.decode("12H12S101M13S13H");
        Assert.assertEquals(SvCigarUtils.getTotalHardClipping(cigar), 25);
        Assert.assertEquals(SvCigarUtils.getNumClippedBases(true, cigar), 24);
        Assert.assertEquals(SvCigarUtils.getNumClippedBases(false, cigar), 26);
    }

    @Test(expectedExceptions=IllegalArgumentException.class, groups = "sv")
    public void testCigarChecker_emptyCigarElementList(){
        @SuppressWarnings("unchecked")
        final List<CigarElement> emptyList = Collections.EMPTY_LIST;
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
    public void testGetNumClippingBases_hardAndSoftSeparately() {
        List<CigarElement> cigarElements = TextCigarCodec.decode("10H20S30M40D50M60S70H").getCigarElements();
        Assert.assertEquals(SvCigarUtils.getNumHardClippingBases(true, cigarElements), 10);
        Assert.assertEquals(SvCigarUtils.getNumHardClippingBases(false, cigarElements), 70);
        Assert.assertEquals(SvCigarUtils.getNumSoftClippingBases(true, cigarElements), 20);
        Assert.assertEquals(SvCigarUtils.getNumSoftClippingBases(false, cigarElements), 60);
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

    @DataProvider(name = "CigarTestData")
    private Object[][] createCigarTestData() {
        final Cigar[] cigars = Stream.of("5M5H", "5H5M", "14S6M", "14M6S",  "12H12S20M51I30M13S13H", "10H20S30M40D50M60S70H").map(TextCigarCodec::decode).toArray(Cigar[]::new);
        final int[] lengths = {10, 10, 20, 20, 151, 240};
        final Object[][] data = new Object[cigars.length][];
        for (int i = 0; i < cigars.length; ++i) {
            data[i] = new Object[]{cigars[i], lengths[i]};
        }
        return data;
    }

    @Test(dataProvider = "CigarTestData", groups = "sv")
    public void testGetUnclippedReadLengthFromCigar(final Cigar cigar, final int expectedReadLength) {
        Assert.assertEquals(SvCigarUtils.getUnclippedReadLength(cigar), expectedReadLength);
    }
}
