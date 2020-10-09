package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.apache.commons.math3.special.Beta;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
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
     * Note that, as in the docs, this assumes the approximation that ref reads have infinite quality -- that is, we
     * don't try to squeeze out the last bit of extra variant likelihood by accounting for the chance that some alt reads
     * were actually misread as ref reads.
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

    /**
     * We can also (again using the perfect ref reads approximation) exactly calculate the integral when the alt count is small.
     * Assuming a constant error rate eps for the alt reads we have a variant likelihood
     * int_(0 to 1) (1-f)^N_ref * [f(1-eps_j) + (1-f)eps_j]^N_alt df
     *
     * We can expand the binomial raised to the N_alt power explicitly for small N_alt to obtain the sum of N_alt terms
     * in the integrand, each of which is the normalization constant of a Beta distribution with integer shape parameters.
     * @param numRef
     * @param errorRate
     */
    @Test(dataProvider = "fewAltData")
    public void testSmallNumAltExact(final int numRef, final double errorRate) {

        for (final int numAlt : new int[] {0, 1, 2}) {
            final double calculated = Mutect2Engine.logLikelihoodRatio(numRef, numAlt, errorRate);
            final double expected;
            switch (numAlt) {
                case 0:
                    expected = Math.log(1.0 / (numRef + 1));;
                    break;
                case 1:
                    expected = Math.log((errorRate / (numRef + 2)) + (1 - errorRate) /((numRef + 1) * (numRef + 2)))
                            - Math.log(errorRate);
                    break;
                case 2:
                    expected = Math.log( (errorRate*errorRate/(numRef + 3)) + (2*errorRate*(1-errorRate)/((numRef+2)*(numRef+3)))
                            + (2*(1-errorRate)*(1-errorRate)/((numRef+1)*(numRef+2)*(numRef+3)))) - 2*Math.log(errorRate);
                    break;
                default:
                    throw new GATKException.ShouldNeverReachHereException("Didn't write this test case");

            }

            // we don't really care about high accuracy if things are obvious:
            if (expected < - 3.0) {
                Assert.assertTrue(calculated < -1.0);
                continue;
            }

            final double precision;
            if (expected < -2) {
                precision = 2.0;
            } else if (expected < 0) {
                precision = 1.0;
            } else if (expected < 1) {
                precision = 0.5;
            } else {
                precision = 0.25;
            }

            Assert.assertEquals(calculated, expected, precision);
        }



    }

    @DataProvider(name = "fewAltData")
    public Object[][] getFewAltData() {
        return new Object[][] {
                { 1, 0.0001},
                { 5, 0.0001},
                { 10, 0.0001},
                { 100, 0.0001},
                { 1000, 0.0001},
                { 1, 0.001},
                { 5, 0.001},
                { 10, 0.001},
                { 100, 0.001},
                { 1000, 0.001},
                { 1, 0.01},
                { 5, 0.01},
                { 10, 0.01},
                { 100, 0.01},
                { 1000, 0.01},
        };
    }

}