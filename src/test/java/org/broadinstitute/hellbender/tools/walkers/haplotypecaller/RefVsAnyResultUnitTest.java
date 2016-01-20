package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.testng.Assert;
import org.testng.annotations.Test;

public final class RefVsAnyResultUnitTest {

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNotNegative() throws Exception {
        new RefVsAnyResult(-1);
    }
    @Test
    public void testConstructor() throws Exception {
        final RefVsAnyResult res = new RefVsAnyResult(3);
        Assert.assertEquals(res.getAD().length, 2);
        Assert.assertEquals(res.getGenotypeLikelihoodsCappedByHomRefLikelihood().length, 3);
        Assert.assertEquals(res.getDP(), 0);
        Assert.assertEquals(res.getAD(), new int[]{0, 0});
        Assert.assertEquals(res.getGenotypeLikelihoodsCappedByHomRefLikelihood(), new double[]{0.0, 0.0, 0.0});
    }

    @Test
    public void testADInc() throws Exception {
        final RefVsAnyResult res = new RefVsAnyResult(3);
        res.incrementRefAD(2);
        res.incrementNonRefAD(3);
        Assert.assertEquals(res.getAD(), new int[]{2, 3});
    }

    @Test
    public void testCapByHomRefLikelihood() throws Exception {
        final RefVsAnyResult res = new RefVsAnyResult(3);
        res.addGenotypeLikelihood(0, 100);
        res.addGenotypeLikelihood(1, 200);
        res.addGenotypeLikelihood(2, 60);
        Assert.assertEquals(res.getGenotypeLikelihoodsCappedByHomRefLikelihood(), new double[]{100.0, 100.0, 60.0});
    }

    @Test
    public void testArrays() throws Exception {
        final RefVsAnyResult res = new RefVsAnyResult(3);
        res.incrementRefAD(2);
        res.incrementNonRefAD(3);
        final int[] adArray = res.getAD();
        Assert.assertEquals(adArray, new int[]{2, 3});

        adArray[0] = 17;
        Assert.assertEquals(res.getAD(), new int[]{2, 3}); //verify that the ad array is a copy

        Assert.assertEquals(res.getGenotypeLikelihoodsCappedByHomRefLikelihood(), new double[]{0, 0, 0});
        double[] liks = res.getGenotypeLikelihoodsCappedByHomRefLikelihood();
        liks[0] = 19;
        Assert.assertEquals(res.getGenotypeLikelihoodsCappedByHomRefLikelihood(), new double[]{0, 0, 0}); //verify that the GL array is a copy
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testArgumentCheckRefAD() throws Exception {
        final RefVsAnyResult res = new RefVsAnyResult(3);
        res.incrementRefAD(-2);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testArgumentCheckNonRefAD() throws Exception {
        final RefVsAnyResult res = new RefVsAnyResult(3);
        res.incrementNonRefAD(-2);
    }

}
