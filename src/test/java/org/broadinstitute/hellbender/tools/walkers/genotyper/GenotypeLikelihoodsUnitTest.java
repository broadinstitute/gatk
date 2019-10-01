package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;

public class GenotypeLikelihoodsUnitTest {

    @Test
    public void testCalcNumLikelihoods() {
        int val4_2 =  GenotypeLikelihoods.numLikelihoods(6, 2);
        Assert.assertEquals(val4_2,21 );
        int val20_2 =  GenotypeLikelihoods.numLikelihoods(20, 2);
        Assert.assertEquals(val20_2,210 );
        int val10_7 =  GenotypeLikelihoods.numLikelihoods(10, 7);
        Assert.assertEquals(val10_7,11440 );
    }

}
