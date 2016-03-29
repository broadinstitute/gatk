package org.broadinstitute.hellbender.tools.exome.allelefraction;

import org.apache.commons.math3.special.Gamma;
import org.broadinstitute.hellbender.tools.exome.AllelicCount;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;

import static java.lang.Math.log;

/**
 * Tests the likelihood methods of the allele-fraction model. <p>
 *
 * We subject the likelihood (with bias marginalized out) to several tests.  Recall that the exact likelihood is:
 * <p>
 * if indicator == ALT_MINOR:
 * <p>
 * log { [beta^alpha / Gamma(alpha)][(1-pi)/2] * int_{0 to infty} f^a * (1-f)^r * lambda^(alpha + r - 1) * exp(-beta*lambda)/(f + (1-f)*lambda)^n d lambda }
 * ,<p> if indicator == REF_MINOR same as ALT_MINOR but with f <--> 1 - f </p>
 * ,<p>if indicator == OUTLIER, log {pi * a!r!/(n+1)!} </p>
 * where alpha = mu*beta and n = a + r.
 * @author David Benjamin
 */
public class AlleleFractionLikelihoodsUnitTest {
    private static final SimpleInterval DUMMY = new SimpleInterval("dummy", 1, 2);
    private static final double EPSILON = 1e-10;

    //if f is very close to 0 we have an analytic result for comparison
    @Test
    public void testHetLogLikelihoodMinorFractionNearZero() {
        final double pi = 0.01; //pi is just a prefactor so we don't need to test it thoroughly here
        for (final double f : Arrays.asList(1e-6, 1e-7, 1e-8)) {
            for (final double mean : Arrays.asList(0.9, 1.0, 1.1)) {
                for (final double variance : Arrays.asList(0.01, 0.005, 0.001)) {
                    final double alpha = mean * mean / variance;
                    final double beta = mean / variance;
                    final AlleleFractionState state = new AlleleFractionState(mean, variance, pi, f);
                    for (final int a : Arrays.asList(1, 2, 3)) {  //alt count
                        for (final int r : Arrays.asList(50, 100, 200)) { //ref count
                            final AllelicCount count = new AllelicCount(DUMMY, r, a);
                            final double actual = AlleleFractionLikelihoods.hetLogLikelihood(state, 0, count, AlleleFractionIndicator.ALT_MINOR);
                            final double expected = a * log(beta) + Gamma.logGamma(alpha - a) - Gamma.logGamma(alpha)
                                    + log((1 - pi) / 2) + a * log(f / (1 - f));
                            Assert.assertEquals(actual, expected, 1e-3);
                        }
                    }
                }
            }
        }
    }

    // if f is very close to 1 we have an analytic result for comparison
    @Test
    public void testHetLogLikelihoodMinorFractionNearOne() {
        final double pi = 0.01; //pi is just a prefactor so we don't need to test it thoroughly here
        for (final double f : Arrays.asList(1 - 1e-6, 1 - 1e-7, 1 - 1e-8)) {
            for (final double mean : Arrays.asList(0.9, 1.0, 1.1)) {
                for (final double variance : Arrays.asList(0.01, 0.005, 0.001)) {
                    final double alpha = mean * mean / variance;
                    final double beta = mean / variance;
                    final AlleleFractionState state = new AlleleFractionState(mean, variance, pi, f);
                    for (final int a : Arrays.asList(1, 10, 20)) {  //alt count
                        for (final int r : Arrays.asList(1, 10, 20)) { //ref count
                            final AllelicCount count = new AllelicCount(DUMMY, r, a);
                            final double actual = AlleleFractionLikelihoods.hetLogLikelihood(state, 0, count, AlleleFractionIndicator.ALT_MINOR);
                            final double expected = -r * log(beta) + Gamma.logGamma(alpha + r) - Gamma.logGamma(alpha)
                                    + log((1 - pi) / 2) - r * log(f / (1 - f));
                            Assert.assertEquals(actual, expected,1e-4);
                        }
                    }
                }
            }
        }
    }

    //if variance is tiny we can approximate lambda ~ mu in the tricky part of the integral to get an analytic result
    @Test
    public void testHetLogLikelihoodTightDistribution() {
        final double pi = 0.01; //pi is just a prefactor so we don't need to test it thoroughly here
        for (final double f : Arrays.asList(0.1, 0.2, 0.3)) {
            for (final double mean : Arrays.asList(0.9, 1.0, 1.1)) {
                for (final double variance : Arrays.asList(1e-6, 1e-7, 1e-8)) {
                    final AlleleFractionState state = new AlleleFractionState(mean, variance, pi, f);
                    for (final int a : Arrays.asList(1, 10, 20)) {  //alt count
                        for (final int r : Arrays.asList(1, 10, 20)) { //ref count
                            final AllelicCount count = new AllelicCount(DUMMY, r, a);
                            final double actual = AlleleFractionLikelihoods.hetLogLikelihood(state, 0, count, AlleleFractionIndicator.ALT_MINOR);
                            final double expected = log((1 - pi) / 2) + a * log(f) + r * log(1-f) + r * log(mean) - (a+r) * log(f + (1-f)*mean);
                            Assert.assertEquals(actual, expected, 1e-3);
                        }
                    }
                }
            }
        }
    }

    //ALT_MINOR <--> REF_MINOR is equivalent to f <--> 1 - f
    @Test
    public void testRefMinor() {
        final double pi = 0.01; //pi is just a prefactor so we don't need to test it thoroughly here
        for (final double f : Arrays.asList(0.1, 0.2, 0.3)) {
            for (final double mean : Arrays.asList(0.9, 1.0, 1.1)) {
                for (final double variance : Arrays.asList(0.02, 0.01)) {
                    final AlleleFractionState altMinorState = new AlleleFractionState(mean, variance, pi, f);
                    final AlleleFractionState refMinorState = new AlleleFractionState(mean, variance, pi, 1 - f);
                    for (final int a : Arrays.asList(1, 10, 20)) {  //alt count
                        for (final int r : Arrays.asList(1, 10, 20)) { //ref count
                            final AllelicCount count = new AllelicCount(DUMMY, r, a);
                            final double altMinorLk = AlleleFractionLikelihoods.hetLogLikelihood(altMinorState, 0, count, AlleleFractionIndicator.ALT_MINOR);
                            final double refMinorLk = AlleleFractionLikelihoods.hetLogLikelihood(refMinorState, 0, count, AlleleFractionIndicator.REF_MINOR);
                            Assert.assertEquals(altMinorLk, refMinorLk, 1e-10);
                        }
                    }
                }
            }
        }
    }

    @Test
    public void testHetLogLikelihoodOutlierProbabilityDependence() {
        final AllelicCount count = new AllelicCount(DUMMY, 11, 37);
        final double f = 0.25;
        final double mean = 1.0;
        final double variance = 0.01;
        final double pi1 = 0.1;
        final double pi2 = 0.2;
        final double pi3 = 0.3;

        final AlleleFractionState state1 = new AlleleFractionState(mean, variance, pi1, f);
        final AlleleFractionState state2 = new AlleleFractionState(mean, variance, pi2, f);
        final AlleleFractionState state3 = new AlleleFractionState(mean, variance, pi3, f);

        final double lk1 = AlleleFractionLikelihoods.hetLogLikelihood(state1, 0, count, AlleleFractionIndicator.ALT_MINOR);
        final double lk2 = AlleleFractionLikelihoods.hetLogLikelihood(state2, 0, count, AlleleFractionIndicator.ALT_MINOR);
        final double lk3 = AlleleFractionLikelihoods.hetLogLikelihood(state3, 0, count, AlleleFractionIndicator.ALT_MINOR);

        Assert.assertEquals(lk2 - lk1, log(1 - pi2) - log(1 - pi1), EPSILON);
        Assert.assertEquals(lk3 - lk2, log(1 - pi3) - log(1 - pi2), EPSILON);
    }
}