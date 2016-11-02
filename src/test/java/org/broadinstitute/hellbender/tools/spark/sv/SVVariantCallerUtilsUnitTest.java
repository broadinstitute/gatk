package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;


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
}
