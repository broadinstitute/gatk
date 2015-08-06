package org.broadinstitute.hellbender.utils.genotyper;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

public final class MostLikelyAlleleUnitTest extends BaseTest {
    final Allele a = Allele.create("A");
    final Allele b = Allele.create("C");

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testErrorMostLikely() throws Exception {
        final double second = -1 - MostLikelyAllele.INFORMATIVE_LIKELIHOOD_THRESHOLD - 1;
        new MostLikelyAllele(a, b, 0.5, second);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testErrorSecondMostLikely() throws Exception {
        new MostLikelyAllele(a, b, -1.0, 0.5);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testErrorMostLikelyLessLikelyThanSecondMostLikely() throws Exception {
        new MostLikelyAllele(a, b, -2.0, -1.0);
    }

    @Test
    public void testBasicCreation() {
        final double second = -1 - MostLikelyAllele.INFORMATIVE_LIKELIHOOD_THRESHOLD - 1;
        MostLikelyAllele mla = new MostLikelyAllele(a, b, -1.0, second);
        Assert.assertEquals(mla.getMostLikelyAllele(), a);
        Assert.assertEquals(mla.getSecondMostLikelyAllele(), b);
        Assert.assertEquals(mla.getLog10LikelihoodOfMostLikely(), -1.0);
        Assert.assertEquals(mla.getLog10LikelihoodOfSecondBest(), second);

        Assert.assertEquals(mla.isInformative(), true);
        Assert.assertEquals(mla.isInformative(10), false);
        Assert.assertEquals(mla.isInformative(0), true);
        Assert.assertEquals(mla.getAlleleIfInformative(), a);
        Assert.assertEquals(mla.getAlleleIfInformative(10), Allele.NO_CALL);
        Assert.assertEquals(mla.getAlleleIfInformative(0), a);
    }

    @Test
    public void testNotDefaultInformative() {
        final double second = -1.0 - (MostLikelyAllele.INFORMATIVE_LIKELIHOOD_THRESHOLD - 1e-2);
        MostLikelyAllele mla = new MostLikelyAllele(a, b, -1.0, second);
        Assert.assertEquals(mla.isInformative(), false);
        Assert.assertEquals(mla.isInformative(10), false);
        Assert.assertEquals(mla.isInformative(0), true);
        Assert.assertEquals(mla.getAlleleIfInformative(), Allele.NO_CALL);
        Assert.assertEquals(mla.getAlleleIfInformative(10), Allele.NO_CALL);
        Assert.assertEquals(mla.getAlleleIfInformative(0), a);
    }

    @Test
    public void testCreationNoGoodSecond() {
        MostLikelyAllele mla = new MostLikelyAllele(a, null, -1.0, Double.NEGATIVE_INFINITY);
        Assert.assertEquals(mla.getMostLikelyAllele(), a);
        Assert.assertEquals(mla.getSecondMostLikelyAllele(), null);
        Assert.assertEquals(mla.getLog10LikelihoodOfMostLikely(), -1.0);
        Assert.assertEquals(mla.getLog10LikelihoodOfSecondBest(), Double.NEGATIVE_INFINITY);

        Assert.assertEquals(mla.isInformative(), true);
        Assert.assertEquals(mla.isInformative(10), true);
        Assert.assertEquals(mla.isInformative(0), true);
        Assert.assertEquals(mla.getAlleleIfInformative(), a);
        Assert.assertEquals(mla.getAlleleIfInformative(10), a);
        Assert.assertEquals(mla.getAlleleIfInformative(0), a);
    }

    @Test
    public void testCreationNoAllele() {
        MostLikelyAllele mla = new MostLikelyAllele(Allele.NO_CALL, null, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY);
        Assert.assertEquals(mla.getMostLikelyAllele(), Allele.NO_CALL);
        Assert.assertEquals(mla.getLog10LikelihoodOfMostLikely(), Double.NEGATIVE_INFINITY);
        Assert.assertEquals(mla.getLog10LikelihoodOfSecondBest(), Double.NEGATIVE_INFINITY);

        Assert.assertEquals(mla.isInformative(), false);
        Assert.assertEquals(mla.isInformative(10), false);
        Assert.assertEquals(mla.isInformative(0), false);
        Assert.assertEquals(mla.getAlleleIfInformative(), Allele.NO_CALL);
        Assert.assertEquals(mla.getAlleleIfInformative(10), Allele.NO_CALL);
        Assert.assertEquals(mla.getAlleleIfInformative(0), Allele.NO_CALL);
    }
}
