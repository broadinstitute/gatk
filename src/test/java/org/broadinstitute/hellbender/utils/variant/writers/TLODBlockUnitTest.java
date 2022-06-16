package org.broadinstitute.hellbender.utils.variant.writers;

import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;

public class TLODBlockUnitTest {
    private static final String SAMPLE_NAME = "foo";
    private static final Allele REF = Allele.create("A", true);
    private static final double TLOD_TO_DP_RATIO = -2.3;  //totally made up

    private static List<Allele> getAlleles() {
        return Arrays.asList(REF, Allele.create("C"));
    }

    private static Genotype getGenotype() {
        return new GenotypeBuilder(SAMPLE_NAME, getAlleles()).make();
    }

    private static VariantContext getVariantContext() {
        return new VariantContextBuilder(SAMPLE_NAME, "20", 1, 1, getAlleles()).genotypes(getGenotype()).make();
    }

    private static TLODBlock getTLODBlock(VariantContext vc) {
        return new TLODBlock(vc, -1650, -450, 0);
    }

    @Test
    public void testAddAndCreateGenotype() {
        final VariantContext vc = getVariantContext();
        final TLODBlock band = getTLODBlock(vc);
        final GenotypeBuilder gb = new GenotypeBuilder("TUMOR");
        gb.alleles(vc.getAlleles());

        final int[] DPs = new int[]{400, 450, 500, 550, 600, 700, 425, 320, 200};
        final double[] TLODs = new double[DPs.length];
        final double[] expectedTLODs = new double[DPs.length];
        final int[] expectedMinDPs = new int[DPs.length];
        final int[] expectedMedianDPs = new int[DPs.length];
        for (int i = 0; i < DPs.length; i++) {
            TLODs[i] = TLOD_TO_DP_RATIO * DPs[i];
            expectedTLODs[i] = Arrays.stream(Arrays.copyOfRange(TLODs, 0, i + 1)).min().getAsDouble();
            expectedMinDPs[i] = MathUtils.arrayMin(Arrays.copyOfRange(DPs, 0, i + 1));
            expectedMedianDPs[i] = MathUtils.median(Arrays.copyOfRange(DPs, 0, i + 1));
        }

        int pos = band.getStart();
        for (int i = 0; i < DPs.length; i++) {
            band.add(++pos, gb.DP(DPs[i]).attribute(GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY, TLODs[i]).make());
            Assert.assertEquals(band.getEnd(), pos);
            Assert.assertEquals(band.getMinDP(), expectedMinDPs[i]);
            Assert.assertEquals(band.getMedianDP(), expectedMedianDPs[i]);
            Assert.assertEquals(band.getMinBlockLOD(), expectedTLODs[i]);
        }

        final Genotype g = band.createHomRefGenotype("TUMOR", false);
        Assert.assertTrue(g.hasDP() && g.getDP() == expectedMedianDPs[DPs.length-1]);
        Assert.assertTrue(g.hasExtendedAttribute(GATKVCFConstants.MIN_DP_FORMAT_KEY) && ((int) g.getExtendedAttribute(GATKVCFConstants.MIN_DP_FORMAT_KEY)) == expectedMinDPs[DPs.length-1]);
        Assert.assertTrue(g.hasExtendedAttribute(GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY) && ((double) g.getExtendedAttribute(GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY)) == expectedTLODs[DPs.length-1]);
    }


    @Test
    public void testWithinBounds() {
        final VariantContext vc = new VariantContextBuilder("TUMOR", "20", 1, 1, getAlleles()).genotypes(getGenotype()).make();
        //Note that the default precision is 1, so the bounds on this block are actually [-0.5, 0)
        TLODBlock band = new TLODBlock(vc, -5, 0, 1);
        Assert.assertFalse(band.withinBounds(0.0));
        Assert.assertTrue(band.withinBounds(-.10));
        Assert.assertTrue(band.withinBounds(-.20));
        Assert.assertTrue(band.withinBounds(-.50));
        Assert.assertFalse(band.withinBounds(-.51));
        Assert.assertFalse(band.withinBounds(-.501));

        //at precision 2 the bounds here are [-0.51, -0.50)
        band = new TLODBlock(vc, -51, -50, 2);
        Assert.assertFalse(band.withinBounds(-.50));
        Assert.assertTrue(band.withinBounds(-.51));
        Assert.assertTrue(band.withinBounds(-.501));
        Assert.assertFalse(band.withinBounds(-.511));
        Assert.assertFalse(band.withinBounds(-.51001));
    }
}