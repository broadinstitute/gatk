package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;

import static org.broadinstitute.hellbender.tools.spark.sv.BreakpointAllele.BreakpointAlleleInversion;


public class SVVariantCallerUtilsUnitTest {

    @Test
    public void testAlignmentRegionOverlap() throws Exception {

        //overlap by 1
        final AlignmentRegion ar1 = new AlignmentRegion("1","1", TextCigarCodec.decode("5M5H"),true,new SimpleInterval("1",1,5),60,1,5,0);
        final AlignmentRegion ar2 = new AlignmentRegion("1","1", TextCigarCodec.decode("5S5M"),true,new SimpleInterval("1",10,16),60,5,10,0);
        Assert.assertEquals(SVVariantCallerUtils.overlapOnContig(ar1, ar2), 1);

        // don't overlap
        final AlignmentRegion ar3 = new AlignmentRegion("1","1", TextCigarCodec.decode("5M5H"),true,new SimpleInterval("1",1,5),60,1,5,0);
        final AlignmentRegion ar4 = new AlignmentRegion("1","1", TextCigarCodec.decode("5S5M"),true,new SimpleInterval("1",11,16),60,6,10,0);
        Assert.assertEquals(SVVariantCallerUtils.overlapOnContig(ar3, ar4), 0);
    }

    @Test
    public void testIsInversion() {

        final AlignmentRegion region1 = new AlignmentRegion("1", "contig-1", TextCigarCodec.decode("100M"), true, new SimpleInterval("1", 10000, 10100), 60, 1, 100, 0);
        final AlignmentRegion region2 = new AlignmentRegion("1", "contig-1", TextCigarCodec.decode("100M"), false, new SimpleInterval("1", 20100, 20200), 60, 101, 200, 0);
        final ChimericAlignment breakpoint1 = new ChimericAlignment(region1, region2, "", "", new ArrayList<>());

        Assert.assertTrue(breakpoint1.involvesStrandSwitch());
        Assert.assertTrue(SVVariantCallerUtils.isInversion(new BreakpointAlleleInversion(breakpoint1)));

        final AlignmentRegion region3 = new AlignmentRegion("4", "contig-7", TextCigarCodec.decode("137M141S"), true, new SimpleInterval("19", 38343346, 38343483), 60, 1, 137, 0);
        final AlignmentRegion region4 = new AlignmentRegion("4", "contig-7", TextCigarCodec.decode("137S141M"), false, new SimpleInterval("10", 38342908, 38343049), 60, 138, 278, 0);
        final ChimericAlignment breakpoint2 = new ChimericAlignment(region3, region4, "", "", new ArrayList<>());

        Assert.assertTrue(breakpoint2.involvesStrandSwitch());
        Assert.assertFalse(SVVariantCallerUtils.isInversion(new BreakpointAlleleInversion(breakpoint2)));

        final AlignmentRegion region5 = new AlignmentRegion("3", "contig-7", TextCigarCodec.decode("137M141S"), true, new SimpleInterval("19", 38343346, 38343483), 60, 1, 137, 0);
        final AlignmentRegion region6 = new AlignmentRegion("3", "contig-7", TextCigarCodec.decode("137S141M"), false, new SimpleInterval("19", 38342908, 38343049), 60, 138, 278, 0);
        final ChimericAlignment breakpoint3 = new ChimericAlignment(region5, region6, "", "", new ArrayList<>());

        Assert.assertTrue(breakpoint3.involvesStrandSwitch());
        Assert.assertTrue(SVVariantCallerUtils.isInversion(new BreakpointAlleleInversion(breakpoint3)));
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
