package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import org.testng.Assert;
import org.testng.annotations.Test;

public final class ExactACcountsUnitTest {
    @Test
    public void test1() throws Exception {
        final int[] c1 = {1, 2, 3};
        final int[] c1equal = {1, 2, 3};
        final int[] c2 = {5,6,7};
        final ExactACcounts ec1 = new ExactACcounts(c1);
        final ExactACcounts ec1Same = new ExactACcounts(c1);//same array
        final ExactACcounts ec1Equal = new ExactACcounts(c1equal);
        final ExactACcounts ec2 = new ExactACcounts(c2);

        Assert.assertEquals(ec1, ec1);
        Assert.assertEquals(ec1Equal, ec1Equal);
        Assert.assertEquals(ec2, ec2);

        Assert.assertTrue(ec1.getCounts() == ec1Same.getCounts());
        Assert.assertFalse(ec1.getCounts() == ec1Equal.getCounts());
        Assert.assertEquals(ec1.getCounts(), ec1Equal.getCounts());

        Assert.assertEquals(ec1, ec1Equal);
        Assert.assertEquals(ec1Equal, ec1);
        Assert.assertFalse(ec1 == ec1Equal);
        Assert.assertEquals(ec1.hashCode(), ec1Equal.hashCode());

        Assert.assertNotEquals(ec1, ec2);
        Assert.assertNotEquals(ec2, ec1);
        Assert.assertNotEquals(ec1Equal, ec2);
        Assert.assertNotEquals(ec2, ec1Equal);

        Assert.assertNotNull(ec1.toString());
        Assert.assertNotNull(ec1Equal.toString());
        Assert.assertNotNull(ec2.toString());
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testEmpty() throws Exception {
        new ExactACcounts(new int[0]);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNull() throws Exception {
        new ExactACcounts(null);
    }

}
