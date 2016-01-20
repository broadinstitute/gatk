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
        Assert.assertEquals(res.get_AD_Ref_Any().length, 2);
        Assert.assertEquals(res.getGenotypeLikelihoods().length, 3);
        Assert.assertEquals(res.getDP(), 0);
        Assert.assertEquals(res.get_AD_Ref_Any(), new int[]{0, 0});
        Assert.assertEquals(res.getGenotypeLikelihoods(), new double[]{0.0, 0.0, 0.0});
    }

    @Test
    public void testADInc() throws Exception {
        final RefVsAnyResult res = new RefVsAnyResult(3);
        res.incrementRefAD(2);
        res.incrementNonRefAD(3);
        Assert.assertEquals(res.get_AD_Ref_Any(), new int[]{2, 3});
    }

    @Test
    public void testCapByHomRefLikelihood() throws Exception {
        final RefVsAnyResult res = new RefVsAnyResult(3);
        res.increaseGenotypeLikelihood(0, 100);
        res.increaseGenotypeLikelihood(1, 200);
        res.increaseGenotypeLikelihood(2, 60);
        Assert.assertEquals(res.getGenotypeLikelihoods(), new double[]{100.0, 200.0, 60.0});
        res.capByHomRefLikelihood();
        Assert.assertEquals(res.getGenotypeLikelihoods(), new double[]{100.0, 100.0, 60.0});
    }

}
