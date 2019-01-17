package org.broadinstitute.hellbender.utils.mcmc;


import org.apache.commons.math3.distribution.TDistribution;
import org.apache.commons.math3.primes.Primes;
import org.apache.commons.math3.random.RandomGenerator;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.*;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Implements slice sampling of a continuous, univariate, unnormalized probability density function (PDF),
 * which is assumed to be unimodal.  See Neal 2003 at https://projecteuclid.org/euclid.aos/1056562461 for details.
 * Minibatching is implemented as in Dubois et al. 2014 at http://proceedings.mlr.press/v33/dubois14.pdf and requires
 * that the PDF, which is assumed to be a posterior function of a parameter value and the data, is specified in terms
 * of a prior, a likelihood, and the data.
 *
 * Note that likelihoods are cached and retrieved based on the hash of {@code DATA}, i.e., we assume it is
 * cheaper to compute this hash rather than the likelihood.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class MinibatchSliceSampler<DATA> extends AbstractSliceSampler {
    private final List<DATA> data;
    private final Function<Double, Double> logPrior;
    private final BiFunction<DATA, Double, Double> logLikelihood;
    private final Integer minibatchSize;
    private final Double approxThreshold;

    private final int numDataPoints;

    private Double xSampleCache = null;
    private Double logPriorCache = null;
    private Map<DATA, Double> logLikelihoodsCache = null;    //data -> log likelihood

    /**
     * Creates a new sampler for a bounded univariate random variable, given a random number generator, a list of data,
     * a continuous, univariate, unimodal, unnormalized log probability density function
     * (assumed to be a posterior and specified by a prior and a likelihood),
     * hard limits on the random variable, a step width, a minibatch size, and a minibatch approximation threshold.
     * @param rng                       random number generator, never {@code null}
     * @param data                      list of data, never {@code null}
     * @param logPrior                  log prior component of continuous, univariate, unimodal log posterior (up to additive constant), never {@code null}
     * @param logLikelihood             log likelihood component of continuous, univariate, unimodal log posterior (up to additive constant), never {@code null}
     * @param xMin                      minimum allowed value of the random variable
     * @param xMax                      maximum allowed value of the random variable
     * @param width                     step width for slice expansion
     * @param minibatchSize             minibatch size
     * @param approxThreshold           threshold for approximation used in {@link MinibatchSliceSampler#isGreaterThanSliceHeight};
     *                                  approximation is exact when this threshold is zero
     */
    public MinibatchSliceSampler(final RandomGenerator rng,
                                 final List<DATA> data,
                                 final Function<Double, Double> logPrior,
                                 final BiFunction<DATA, Double, Double> logLikelihood,
                                 final double xMin,
                                 final double xMax,
                                 final double width,
                                 final int minibatchSize,
                                 final double approxThreshold) {
        super(rng, xMin, xMax, width);
        Utils.nonNull(data);
        Utils.nonNull(logPrior);
        Utils.nonNull(logLikelihood);
        Utils.validateArg(minibatchSize > 1, "Minibatch size must be greater than 1.");
        ParamUtils.isPositiveOrZero(approxThreshold, "Minibatch approximation threshold must be non-negative.");
        this.data = Collections.unmodifiableList(new ArrayList<>(data));
        this.logPrior = logPrior;
        this.logLikelihood = logLikelihood;
        this.minibatchSize = minibatchSize;
        this.approxThreshold = approxThreshold;
        numDataPoints = data.size();
    }

    /**
     * Creates a new sampler for an unbounded univariate random variable, given a random number generator, a list of data,
     * a continuous, univariate, unimodal, unnormalized log probability density function
     * (assumed to be a posterior and specified by a prior and a likelihood),
     * a step width, a minibatch size, and a minibatch approximation threshold.
     * @param rng                       random number generator, never {@code null}
     * @param data                      list of data, never {@code null}
     * @param logPrior                  log prior component of continuous, univariate, unimodal log posterior (up to additive constant), never {@code null}
     * @param logLikelihood             log likelihood component of continuous, univariate, unimodal log posterior (up to additive constant), never {@code null}
     * @param width                     step width for slice expansion
     * @param minibatchSize             minibatch size
     * @param approxThreshold           threshold for approximation used in {@link MinibatchSliceSampler#isGreaterThanSliceHeight};
     *                                  approximation is exact when this threshold is zero
     */
    public MinibatchSliceSampler(final RandomGenerator rng,
                                 final List<DATA> data,
                                 final Function<Double, Double> logPrior,
                                 final BiFunction<DATA, Double, Double> logLikelihood,
                                 final double width,
                                 final int minibatchSize,
                                 final double approxThreshold) {
        this(rng, data, logPrior, logLikelihood, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, width, minibatchSize, approxThreshold);
    }

    /**
     * Implements the OnSlice procedure from Dubois et al. 2014.
     */
    @Override
    boolean isGreaterThanSliceHeight(final double xProposed,
                                     final double xSample,
                                     final double z) {
        if (xProposed < xMin || xMax < xProposed) {
            return false;
        }

        //we cache values calculated from xSample, since this method is called multiple times for the same value
        //of xSample when expanding slice interval and proposing samples
        if (xSampleCache == null || xSampleCache != xSample) {
            xSampleCache = xSample;
            logPriorCache = logPrior.apply(xSample);
            logLikelihoodsCache = new HashMap<>(numDataPoints);
        }
        if (!(xSampleCache != null && logPriorCache != null && logLikelihoodsCache != null)) {
            throw new GATKException.ShouldNeverReachHereException("Cache for xSample is in an invalid state.");
        }

        //if no data points are provided, sample from the prior
        if (numDataPoints == 0) {
            return logPrior.apply(xProposed) > logPriorCache - z;
        }

        //see Eq. 1 of Dubois et al. 2014
        final double mu0 = (logPriorCache - logPrior.apply(xProposed) - z) / numDataPoints;

        //initialize the lazy data iterator (or just use the standard iterator if only a single batch is required)
        final int numMinibatches = Math.max(numDataPoints / minibatchSize, 1);
        final Iterator<DATA> shuffledDataIterator = numMinibatches > 1
                ? lazyShuffleIterator(rng, data)
                : data.iterator();

        //initialize running quantities needed for statistical test
        int numDataIndicesSeen = 0;
        double logLikelihoodDifferencesMean = 0.;
        double logLikelihoodDifferencesSquaredMean = 0.;

        for (int minibatchIndex = 0; minibatchIndex < numMinibatches; minibatchIndex++) {
            //get a minibatch of data
            final int dataIndexStart = minibatchIndex * minibatchSize;
            final int dataIndexEnd = Math.min((minibatchIndex + 1) * minibatchSize, numDataPoints);
            final int actualMinibatchSize = dataIndexEnd - dataIndexStart;  //equals minibatchSize except perhaps for last minibatch
            final List<DATA> dataMinibatch = IntStream.range(0, actualMinibatchSize).boxed()
                    .map(i -> shuffledDataIterator.next())
                    .collect(Collectors.toList());

            //calculate quantities for this minibatch
            double logLikelihoodDifferencesMinibatchSum = 0.;
            double logLikelihoodDifferencesSquaredMinibatchSum = 0.;
            for (final DATA dataPoint : dataMinibatch) {
                final double logLikelihoodxSample = logLikelihoodsCache.computeIfAbsent(
                        dataPoint, d -> logLikelihood.apply(d, xSample));
                final double logLikelihoodxProposed = logLikelihood.apply(dataPoint, xProposed);
                final double logLikelihoodDifference = logLikelihoodxProposed - logLikelihoodxSample;
                logLikelihoodDifferencesMinibatchSum += logLikelihoodDifference;
                logLikelihoodDifferencesSquaredMinibatchSum += logLikelihoodDifference * logLikelihoodDifference;
            }

            //update running quantities
            logLikelihoodDifferencesMean =
                    (numDataIndicesSeen * logLikelihoodDifferencesMean + logLikelihoodDifferencesMinibatchSum) /
                            (numDataIndicesSeen + actualMinibatchSize);
            logLikelihoodDifferencesSquaredMean =
                    (numDataIndicesSeen * logLikelihoodDifferencesSquaredMean + logLikelihoodDifferencesSquaredMinibatchSum) /
                            (numDataIndicesSeen + actualMinibatchSize);
            numDataIndicesSeen += actualMinibatchSize;

            //if there is only a single batch, then there is no need to perform the statistical test
            //(as all running quantities will be exact)
            if (numMinibatches == 1) {
                break;
            }

            //perform the statistical test and terminate if appropriate
            //(i.e., if we can determine if we are OnSlice to within the approximation threshold)
            final double s = Math.sqrt(1. - (double) numDataIndicesSeen / numDataPoints) *
                    Math.sqrt((logLikelihoodDifferencesSquaredMean - Math.pow(logLikelihoodDifferencesMean, 2)) / (numDataIndicesSeen - 1));
            final double delta = 1. - new TDistribution(null, numDataIndicesSeen - 1)
                    .cumulativeProbability(Math.abs((logLikelihoodDifferencesMean - mu0) / s));
            if (delta < approxThreshold) {
                break;
            }
        }
        return logLikelihoodDifferencesMean > mu0;
    }

    /**
     * To efficiently sample without replacement with the possibility of early stopping when creating minibatches,
     * we lazily shuffle to avoid unnecessarily shuffling all data.  Uses the properties of relative primes and is
     * random enough for our purposes.  Adapted from https://stackoverflow.com/questions/16165128/lazy-shuffle-algorithms.
     */
    private static <T> Iterator<T> lazyShuffleIterator(final RandomGenerator rng,
                                                       final List<T> data) {
        final int numDataPoints = data.size();

        //find first prime greater than or equal to numDataPoints
        final int nextPrime = Primes.nextPrime(numDataPoints);

        return new Iterator<T>() {
            int numSeen = 0;
            int index = rng.nextInt(numDataPoints) + 1;
            final int increment = index;

            public boolean hasNext() {
                return numSeen < data.size();
            }

            @Override
            public T next() {
                while (true) {
                    index = (index + increment) % nextPrime;
                    if (index < numDataPoints) {
                        numSeen++;
                        return data.get(index);
                    }
                }
            }
        };
    }
}