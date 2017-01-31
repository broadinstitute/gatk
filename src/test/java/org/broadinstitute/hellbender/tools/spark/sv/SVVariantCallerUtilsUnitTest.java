package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Collections;
import java.util.List;


public class SVVariantCallerUtilsUnitTest {

    @Test
    public void testAlignmentRegionOverlap() throws Exception {

        //overlap by 1
        final AlignmentRegion ar1 = new AlignmentRegion("1","1", new SimpleInterval("1",1,5), TextCigarCodec.decode("5M5H"),true, 60, 0, 1,5);
        final AlignmentRegion ar2 = new AlignmentRegion("1","1", new SimpleInterval("1",10,16), TextCigarCodec.decode("5S5M"),true, 60, 0, 5,10);
        Assert.assertEquals(SVVariantCallerUtils.overlapOnContig(ar1, ar2), 1);

        // don't overlap
        final AlignmentRegion ar3 = new AlignmentRegion("1","1", new SimpleInterval("1",1,5), TextCigarCodec.decode("5M5H"),true, 60, 0, 1,5);
        final AlignmentRegion ar4 = new AlignmentRegion("1","1", new SimpleInterval("1",11,16), TextCigarCodec.decode("5S5M"),true, 60, 0, 6,10);
        Assert.assertEquals(SVVariantCallerUtils.overlapOnContig(ar3, ar4), 0);
    }

    @Test
    public void testClippingArithmetic() {
        Cigar cigar = TextCigarCodec.decode("100M51S");
        Assert.assertEquals(SVVariantCallerUtils.getTotalHardClipping(cigar), 0);
        Assert.assertEquals(SVVariantCallerUtils.getNumClippedBases(true, cigar), 0);
        Assert.assertEquals(SVVariantCallerUtils.getNumClippedBases(false, cigar), 51);

        cigar = TextCigarCodec.decode("51S100M");
        Assert.assertEquals(SVVariantCallerUtils.getTotalHardClipping(cigar), 0);
        Assert.assertEquals(SVVariantCallerUtils.getNumClippedBases(true, cigar), 51);
        Assert.assertEquals(SVVariantCallerUtils.getNumClippedBases(false, cigar), 0);

        cigar = TextCigarCodec.decode("100M51H");
        Assert.assertEquals(SVVariantCallerUtils.getTotalHardClipping(cigar), 51);
        Assert.assertEquals(SVVariantCallerUtils.getNumClippedBases(true, cigar), 0);
        Assert.assertEquals(SVVariantCallerUtils.getNumClippedBases(false, cigar), 51);

        cigar = TextCigarCodec.decode("51H100M");
        Assert.assertEquals(SVVariantCallerUtils.getTotalHardClipping(cigar), 51);
        Assert.assertEquals(SVVariantCallerUtils.getNumClippedBases(true, cigar), 51);
        Assert.assertEquals(SVVariantCallerUtils.getNumClippedBases(false, cigar), 0);

        cigar = TextCigarCodec.decode("12H12S101M13S13H");
        Assert.assertEquals(SVVariantCallerUtils.getTotalHardClipping(cigar), 25);
        Assert.assertEquals(SVVariantCallerUtils.getNumClippedBases(true, cigar), 24);
        Assert.assertEquals(SVVariantCallerUtils.getNumClippedBases(false, cigar), 26);
    }

    @Test(expectedExceptions=IllegalArgumentException.class)
    public void testCigarChecker_emptyCigarElementList(){
        @SuppressWarnings("unchecked")
        final List<CigarElement> emptyList = Collections.EMPTY_LIST;
        SVVariantCallerUtils.validateCigar(emptyList);
    }

    @Test(expectedExceptions=IllegalArgumentException.class)
    public void testCigarChecker_deletionNeighboringClipping(){
        SVVariantCallerUtils.validateCigar(TextCigarCodec.decode("10S10D10M").getCigarElements());
        SVVariantCallerUtils.validateCigar(TextCigarCodec.decode("10H10S10D10M").getCigarElements());
        SVVariantCallerUtils.validateCigar(TextCigarCodec.decode("10M10D10S").getCigarElements());
        SVVariantCallerUtils.validateCigar(TextCigarCodec.decode("10M10D10S10H").getCigarElements());
    }

    @Test(expectedExceptions=IllegalArgumentException.class)
    public void testCigarChecker_only1NonAlignment(){
        SVVariantCallerUtils.validateCigar(TextCigarCodec.decode("10S").getCigarElements());
    }

    @Test(expectedExceptions=IllegalArgumentException.class)
    public void testCigarChecker_noAlignment(){
        SVVariantCallerUtils.validateCigar(TextCigarCodec.decode("10H10S10I10S10H").getCigarElements());
    }

    @Test
    public void testGetNumClippingBases_hardAndSoftSeparately() {
        List<CigarElement> cigarElements = TextCigarCodec.decode("10H20S30M40D50M60S70H").getCigarElements();
        Assert.assertEquals(SVVariantCallerUtils.getNumHardClippingBases(true, cigarElements), 10);
        Assert.assertEquals(SVVariantCallerUtils.getNumHardClippingBases(false, cigarElements), 70);
        Assert.assertEquals(SVVariantCallerUtils.getNumSoftClippingBases(true, cigarElements), 20);
        Assert.assertEquals(SVVariantCallerUtils.getNumSoftClippingBases(false, cigarElements), 60);
    }

    @Test
    public void testGetIndexOfFirstNonClippingBase(){
        Assert.assertEquals(SVVariantCallerUtils.findIndexOfFirstNonClippingOperation(TextCigarCodec.decode("151M").getCigarElements(), true), 0);
        Assert.assertEquals(SVVariantCallerUtils.findIndexOfFirstNonClippingOperation(TextCigarCodec.decode("151M").getCigarElements(), false), 0);
        Assert.assertEquals(SVVariantCallerUtils.findIndexOfFirstNonClippingOperation(TextCigarCodec.decode("10S10D10M").getCigarElements(), true), 1);
        Assert.assertEquals(SVVariantCallerUtils.findIndexOfFirstNonClippingOperation(TextCigarCodec.decode("10H10S10D10M").getCigarElements(), true), 2);
        Assert.assertEquals(SVVariantCallerUtils.findIndexOfFirstNonClippingOperation(TextCigarCodec.decode("10M10D10S").getCigarElements(), false), 1);
        Assert.assertEquals(SVVariantCallerUtils.findIndexOfFirstNonClippingOperation(TextCigarCodec.decode("10M10D10S10H").getCigarElements(), false), 1);
    }
}
