package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.apache.commons.math3.special.Beta;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Collections;
import java.util.List;

public class Mutect2EngineUnitTest extends GATKBaseTest {


    /**
     * Test the active region log likelihood ratio to lowest order.  From the Mutect docs the no-variation likelihood is
     * prod_alt(eps_j), where eps_j is the error probability of read j and alt means "over all alt reads".
     *
     * The likelihood for a variant with allele fraction f is (1-f)^N_ref * prod_alt [f(1-eps_j) + (1-f)eps_j].
     *
     * To leading order in epsilon the log likelihood ratio is
     *
     * log int_{f 0 to 1} (1-f)^N_ref f^N_alt - log prod (eps_j)
     * = log Beta(N_ref + 1, N_alt + 1) - sum_j log eps_j
     *
     * Thus we should get approximately this result when the error rate is very low.
     */
    @Test(dataProvider = "leadingOrderData")
    public void testLeadingOrderInErrorRate(final int numRef, final int numAlt, final double errorRate) {
        final byte qual = QualityUtils.errorProbToQual(errorRate);
        final List<Byte> altQuals = Collections.nCopies(numAlt, qual);

        final double calculated = Mutect2Engine.logLikelihoodRatio(numRef, altQuals, 1);
        final double expected = Beta.logBeta(numRef + 1, numAlt + 1) - numAlt * Math.log(errorRate);
        Assert.assertEquals(calculated, expected, 0.07);

    }

    // test that the repeats option does the same thing as the explicit calculation
    @Test(dataProvider = "leadingOrderData")
    public void testNumRepeats(final int numRef, final int numAlt, final double errorRate) {
        final byte qual = QualityUtils.errorProbToQual(errorRate);
        for (final int numRepeats : new int[] {2, 5}) {
            final List<Byte> altQuals = Collections.nCopies(numAlt, qual);
            final List<Byte> repeatedAltQuals = Collections.nCopies(numRepeats * numAlt, qual);
            final double calculated = Mutect2Engine.logLikelihoodRatio(numRef, altQuals, numRepeats);
            final double explicit = Mutect2Engine.logLikelihoodRatio(numRef, repeatedAltQuals, 1);
            Assert.assertEquals(calculated, explicit, 0.000001);
        }
    }

    @DataProvider(name = "leadingOrderData")
    public Object[][] getLeadingOrderData() {
        return new Object[][] {
                { 100, 5, 0.0001},
                { 100, 20, 0.0001},
                { 100, 200, 0.0001},
                { 10, 2, 0.0001},
                { 1000, 2, 0.00001},
                { 10000, 2, 0.000001}
        };
    }

}