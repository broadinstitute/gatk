package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Collections;
import java.util.List;


public class SVVariantDiscoveryUtilsUnitTest {

    @Test(groups = "sv")
    public void testAlignmentIntervalOverlap() throws Exception {

        final AlignedAssembly.AlignmentInterval ar1 = new AlignedAssembly.AlignmentInterval(new SimpleInterval("1",1,5), 1,5, TextCigarCodec.decode("5M5H"),true, 60, 0);
        final AlignedAssembly.AlignmentInterval ar2 = new AlignedAssembly.AlignmentInterval(new SimpleInterval("1",10,16), 5,10, TextCigarCodec.decode("4S6M"),true, 60, 0);
        Assert.assertEquals(SVVariantDiscoveryUtils.overlapOnContig(ar1, ar2), 1);

        final AlignedAssembly.AlignmentInterval ar3 = new AlignedAssembly.AlignmentInterval(new SimpleInterval("1",1,5), 1,5, TextCigarCodec.decode("5M5H"),true, 60, 0);
        final AlignedAssembly.AlignmentInterval ar4 = new AlignedAssembly.AlignmentInterval(new SimpleInterval("1",11,16), 6,10, TextCigarCodec.decode("5S5M"),true, 60, 0);
        Assert.assertEquals(SVVariantDiscoveryUtils.overlapOnContig(ar3, ar4), 0);
    }

    @Test(groups = "sv")
    public void testClippingArithmetic() {
        Cigar cigar = TextCigarCodec.decode("100M51S");
        Assert.assertEquals(SVVariantDiscoveryUtils.getTotalHardClipping(cigar), 0);
        Assert.assertEquals(SVVariantDiscoveryUtils.getNumClippedBases(true, cigar), 0);
        Assert.assertEquals(SVVariantDiscoveryUtils.getNumClippedBases(false, cigar), 51);

        cigar = TextCigarCodec.decode("51S100M");
        Assert.assertEquals(SVVariantDiscoveryUtils.getTotalHardClipping(cigar), 0);
        Assert.assertEquals(SVVariantDiscoveryUtils.getNumClippedBases(true, cigar), 51);
        Assert.assertEquals(SVVariantDiscoveryUtils.getNumClippedBases(false, cigar), 0);

        cigar = TextCigarCodec.decode("100M51H");
        Assert.assertEquals(SVVariantDiscoveryUtils.getTotalHardClipping(cigar), 51);
        Assert.assertEquals(SVVariantDiscoveryUtils.getNumClippedBases(true, cigar), 0);
        Assert.assertEquals(SVVariantDiscoveryUtils.getNumClippedBases(false, cigar), 51);

        cigar = TextCigarCodec.decode("51H100M");
        Assert.assertEquals(SVVariantDiscoveryUtils.getTotalHardClipping(cigar), 51);
        Assert.assertEquals(SVVariantDiscoveryUtils.getNumClippedBases(true, cigar), 51);
        Assert.assertEquals(SVVariantDiscoveryUtils.getNumClippedBases(false, cigar), 0);

        cigar = TextCigarCodec.decode("12H12S101M13S13H");
        Assert.assertEquals(SVVariantDiscoveryUtils.getTotalHardClipping(cigar), 25);
        Assert.assertEquals(SVVariantDiscoveryUtils.getNumClippedBases(true, cigar), 24);
        Assert.assertEquals(SVVariantDiscoveryUtils.getNumClippedBases(false, cigar), 26);
    }

    @Test(expectedExceptions=IllegalArgumentException.class, groups = "sv")
    public void testCigarChecker_emptyCigarElementList(){
        @SuppressWarnings("unchecked")
        final List<CigarElement> emptyList = Collections.EMPTY_LIST;
        SVVariantDiscoveryUtils.validateCigar(emptyList);
    }

    @Test(expectedExceptions=IllegalArgumentException.class, groups = "sv")
    public void testCigarChecker_deletionNeighboringClipping(){
        SVVariantDiscoveryUtils.validateCigar(TextCigarCodec.decode("10S10D10M").getCigarElements());
        SVVariantDiscoveryUtils.validateCigar(TextCigarCodec.decode("10H10S10D10M").getCigarElements());
        SVVariantDiscoveryUtils.validateCigar(TextCigarCodec.decode("10M10D10S").getCigarElements());
        SVVariantDiscoveryUtils.validateCigar(TextCigarCodec.decode("10M10D10S10H").getCigarElements());
    }

    @Test(expectedExceptions=IllegalArgumentException.class, groups = "sv")
    public void testCigarChecker_only1NonAlignment(){
        SVVariantDiscoveryUtils.validateCigar(TextCigarCodec.decode("10S").getCigarElements());
    }

    @Test(expectedExceptions=IllegalArgumentException.class, groups = "sv")
    public void testCigarChecker_noAlignment(){
        SVVariantDiscoveryUtils.validateCigar(TextCigarCodec.decode("10H10S10I10S10H").getCigarElements());
    }

    @Test(groups = "sv")
    public void testGetNumClippingBases_hardAndSoftSeparately() {
        List<CigarElement> cigarElements = TextCigarCodec.decode("10H20S30M40D50M60S70H").getCigarElements();
        Assert.assertEquals(SVVariantDiscoveryUtils.getNumHardClippingBases(true, cigarElements), 10);
        Assert.assertEquals(SVVariantDiscoveryUtils.getNumHardClippingBases(false, cigarElements), 70);
        Assert.assertEquals(SVVariantDiscoveryUtils.getNumSoftClippingBases(true, cigarElements), 20);
        Assert.assertEquals(SVVariantDiscoveryUtils.getNumSoftClippingBases(false, cigarElements), 60);
    }

    @Test(groups = "sv")
    public void testGetIndexOfFirstNonClippingBase(){
        Assert.assertEquals(SVVariantDiscoveryUtils.findIndexOfFirstNonClippingOperation(TextCigarCodec.decode("151M").getCigarElements(), true), 0);
        Assert.assertEquals(SVVariantDiscoveryUtils.findIndexOfFirstNonClippingOperation(TextCigarCodec.decode("151M").getCigarElements(), false), 0);
        Assert.assertEquals(SVVariantDiscoveryUtils.findIndexOfFirstNonClippingOperation(TextCigarCodec.decode("10S10D10M").getCigarElements(), true), 1);
        Assert.assertEquals(SVVariantDiscoveryUtils.findIndexOfFirstNonClippingOperation(TextCigarCodec.decode("10H10S10D10M").getCigarElements(), true), 2);
        Assert.assertEquals(SVVariantDiscoveryUtils.findIndexOfFirstNonClippingOperation(TextCigarCodec.decode("10M10D10S").getCigarElements(), false), 1);
        Assert.assertEquals(SVVariantDiscoveryUtils.findIndexOfFirstNonClippingOperation(TextCigarCodec.decode("10M10D10S10H").getCigarElements(), false), 1);
    }
}
