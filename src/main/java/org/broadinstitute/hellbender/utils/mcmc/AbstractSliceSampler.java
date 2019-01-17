package org.broadinstitute.hellbender.utils.mcmc;


import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;

/**
 * Abstract class for slice sampling of a continuous, univariate, unnormalized probability density function (PDF),
 * which is assumed to be unimodal.  This class extracts the logic for initializing the slice height and interval
 * (via doubling) and shrinking the interval during proposal of samples.  The abstract method
 * {@link AbstractSliceSampler#isGreaterThanSliceHeight} should be implemented and determines whether the PDF evaluated
 * at proposed interval endpoints or samples is greater than the slice height.
 * See Neal 2003 at https://projecteuclid.org/euclid.aos/1056562461 for details.
 *
 * The standard batch implementation can be found in {@link SliceSampler}.
 *
 * A minibatch implementation can be be found in {@link MinibatchSliceSampler}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
abstract class AbstractSliceSampler {
    private static final int MAXIMUM_NUMBER_OF_DOUBLINGS = 16;
    private static final int MAXIMUM_NUMBER_OF_SLICE_SAMPLINGS = 100;
    private static final double EPSILON = 1E-10;

    final RandomGenerator rng;
    final double xMin;
    final double xMax;
    private final double width;
    private final ExponentialDistribution exponentialDistribution;

    /**
     * Creates a new sampler for a bounded univariate random variable, given a random number generator, hard limits
     * on the random variable, and a step width.
     * @param rng      random number generator, never {@code null}
     * @param xMin     minimum allowed value of the random variable
     * @param xMax     maximum allowed value of the random variable
     * @param width    step width for slice expansion
     */
    AbstractSliceSampler(final RandomGenerator rng,
                         final double xMin,
                         final double xMax,
                         final double width) {
        Utils.nonNull(rng);
        Utils.validateArg(xMin < xMax, "Maximum bound must be greater than minimum bound.");
        ParamUtils.isPositive(width, "Slice-sampling width must be positive.");
        this.rng = rng;
        this.xMin = xMin;
        this.xMax = xMax;
        this.width = width;
        exponentialDistribution = new ExponentialDistribution(rng, 1.);
    }

    /**
     * Generate a single sample, given an initial value to use in slice construction.
     * @param xInitial      initial value to use in slice construction; must be in [xMin, xMax]
     * @return              sample drawn from the probability density function
     */
    public double sample(final double xInitial) {
        Utils.validateArg(xMin <= xInitial && xInitial <= xMax, "Initial point in slice sampler is not within specified range.");

        //adjust xInitial if on boundary
        final double xSample = Math.min(Math.max(xInitial, xMin + EPSILON), xMax - EPSILON);

        //follow Neal 2003 procedure to slice sample with doubling of interval (assuming unimodal distribution)

        //sample the variable used to randomly pick height of slice from uniform distribution under PDF(xSample);
        //slice height = u * PDF(xSample), where u ~ Uniform(0, 1)
        //however, since we are working with logPDF, we instead sample z = -log u ~ Exponential(1)
        final double z = exponentialDistribution.sample();

        //randomly position slice with given width so that it brackets xSample; position is uniformly distributed
        double xLeft = xSample - width * rng.nextDouble();
        double xRight = xLeft + width;

        final Function<Double, Boolean> isGreaterThanSliceHeight =
                xProposed -> isGreaterThanSliceHeight(xProposed, xSample, z);

        int k = MAXIMUM_NUMBER_OF_DOUBLINGS;
        //expand slice interval by doubling until it brackets PDF
        //(i.e., PDF at both ends is less than the slice height)
        while (k > 0 && (isGreaterThanSliceHeight.apply(xLeft) || isGreaterThanSliceHeight.apply(xRight))) {
            if (rng.nextBoolean()) {
                xLeft = xLeft - (xRight - xLeft);
            } else {
                xRight = xRight + (xRight - xLeft);
            }
            k--;
        }

        //sample uniformly from slice interval until sample over slice height found, shrink slice on each iteration if not found
        //limited to MAXIMUM_NUMBER_OF_SLICE_SAMPLINGS, after which last proposed sample is returned
        //(shouldn't happen if width is chosen appropriately)
        int numIterations = 1;
        double xProposed = rng.nextDouble() * (xRight - xLeft) + xLeft;
        while (numIterations <= MAXIMUM_NUMBER_OF_SLICE_SAMPLINGS) {
            if (isGreaterThanSliceHeight.apply(xProposed)) {
                break;
            }
            if (xProposed < xSample) {
                xLeft = xProposed;
            } else {
                xRight = xProposed;
            }
            xProposed = rng.nextDouble() * (xRight - xLeft) + xLeft;
            numIterations++;
        }
        return Math.min(Math.max(xProposed, xMin + EPSILON), xMax - EPSILON);
    }

    /**
     * Generate multiple samples from the probability density function, given an initial value to use in slice construction.
     * @param xInitial      initial value to use in slice construction; if outside [xMin, xMax], forced to be within
     * @param numSamples    number of samples to generate
     * @return              samples drawn from the probability density function
     */
    public List<Double> sample(final double xInitial,
                               final int numSamples) {
        ParamUtils.isPositive(numSamples, "Number of samples must be positive.");
        final List<Double> samples = new ArrayList<>(numSamples);
        double xSample = xInitial;
        for (int i = 0; i < numSamples; i++) {
            xSample = sample(xSample);
            samples.add(xSample);
        }
        return samples;
    }

    /**
     * Returns true if PDF(xProposed) > slice height = u * PDF(xSample), where u = exp(-z);
     * this is equivalent to returning true if logPDF(xProposed) > logPDF(xSample) - z.
     */
    abstract boolean isGreaterThanSliceHeight(final double xProposed,
                                              final double xSample,
                                              final double z);
}