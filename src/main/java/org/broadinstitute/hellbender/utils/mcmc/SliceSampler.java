package org.broadinstitute.hellbender.utils.mcmc;


import org.apache.commons.math3.random.RandomGenerator;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.function.Function;

/**
 * Implements slice sampling of a continuous, univariate, unnormalized probability density function (PDF),
 * which is assumed to be unimodal.  See Neal 2003 at https://projecteuclid.org/euclid.aos/1056562461 for details.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class SliceSampler extends AbstractSliceSampler {
    private final Function<Double, Double> logPDF;

    private Double xSampleCache = null;
    private Double logPDFCache = null;

    /**
     * Creates a new sampler for a bounded univariate random variable, given a random number generator,
     * a continuous, univariate, unimodal, unnormalized log probability density function,
     * hard limits on the random variable, and a step width.
     * @param rng      random number generator, never {@code null}
     * @param logPDF   continuous, univariate, unimodal log probability density function (up to additive constant), never {@code null}
     * @param xMin     minimum allowed value of the random variable
     * @param xMax     maximum allowed value of the random variable
     * @param width    step width for slice expansion
     */
    public SliceSampler(final RandomGenerator rng,
                        final Function<Double, Double> logPDF,
                        final double xMin,
                        final double xMax,
                        final double width) {
        super(rng, xMin, xMax, width);
        Utils.nonNull(logPDF);
        this.logPDF = logPDF;
    }

    /**
     * Creates a new sampler for an unbounded univariate random variable, given a random number generator,
     * a continuous, univariate, unimodal, unnormalized log probability density function,
     * and a step width.
     * @param rng      random number generator, never {@code null}
     * @param logPDF   continuous, univariate, unimodal log probability density function (up to additive constant), never {@code null}
     * @param width    step width for slice expansion
     */
    public SliceSampler(final RandomGenerator rng,
                        final Function<Double, Double> logPDF,
                        final double width) {
        this(rng, logPDF, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, width);
    }

    @Override
    boolean isGreaterThanSliceHeight(final double xProposed,
                                     final double xSample,
                                     final double z) {
        if (xProposed < xMin || xMax < xProposed) {
            return false;
        }

        //we cache logPDF(xSample), since this method is called multiple times for the same value of xSample
        //when expanding slice interval and proposing samples
        if (xSampleCache == null || xSampleCache != xSample) {
            xSampleCache = xSample;
            logPDFCache = logPDF.apply(xSample);
        }
        if (!((xSampleCache == null && logPDFCache == null) ||
                (xSampleCache != null && logPDFCache != null))) {
            throw new GATKException.ShouldNeverReachHereException("Cache for xSample is in an invalid state.");
        }

        return logPDF.apply(xProposed) > logPDFCache - z;
    }
}