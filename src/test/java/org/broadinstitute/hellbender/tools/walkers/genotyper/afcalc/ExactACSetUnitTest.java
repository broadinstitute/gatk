package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import org.broadinstitute.hellbender.utils.MathUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

public final class ExactACSetUnitTest {
    @Test
    public void test1() throws Exception {
        final int[] c1 = {1,2,3};
        final int size1 = 7;
        ExactACcounts ec1= new ExactACcounts(c1);
        ExactACset acs1 = new ExactACset(size1, ec1);
        Assert.assertEquals(acs1.getACsum(), MathUtils.sum(c1));
        Assert.assertEquals(acs1.getACcounts(), ec1);
        Assert.assertEquals(acs1.getLog10Likelihoods().length, size1);

        final int[] c2 = {1,2,3};
        final int size2 = 7;
        ExactACcounts ec2= new ExactACcounts(c2);
        ExactACset acs2 = new ExactACset(size2, ec2);
        Assert.assertEquals(acs1, acs2);
        Assert.assertEquals(acs1.hashCode(), acs2.hashCode());
        Assert.assertEquals(acs2, acs1);

    }
}
