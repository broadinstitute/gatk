package org.broadinstitute.hellbender.utils.mcmc;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class MultivariateMetropolisSampler {
    private final int size;
    private final Random random = new Random();
    private final Function<double[], Double> logLikelihood;
    private final Supplier<double[]> stepSampler;
    private final double[] lowerBound;
    private final double[] upperBound;
    private final Logger logger = LogManager.getLogger(this.getClass());

    public MultivariateMetropolisSampler(final int size,
                                    final Function<double[], Double> logLikelihood,
                                    final Supplier<double[]> stepSampler,
                                    final double[] lowerBound,
                                    final double[] upperBound) {
        Utils.nonNull(logLikelihood, "Likelihood function cannot be null");
        Utils.nonNull(stepSampler, "Step sampler cannot be null");
        Utils.nonNull(lowerBound, "Lower bound array cannot be null");
        Utils.nonNull(upperBound, "Upper bound array cannot be null");
        Utils.validateArg(lowerBound.length == size, "Lower bound array must have length equal to size parameter");
        Utils.validateArg(upperBound.length == size, "Upper bound array must have length equal to size parameter");
        this.size = size;
        this.logLikelihood = logLikelihood;
        this.stepSampler = stepSampler;
        this.lowerBound = lowerBound;
        this.upperBound = upperBound;
        seed(0L);
    }

    public void seed(final long seed) {
        random.setSeed(seed);
    }

    private void updateNeighbor(final double[] x, final double[] neighborX) {
        final double[] steps = stepSampler.get();
        IntStream.range(0, size).forEach(i -> neighborX[i] = x[i] + steps[i]);
    }

    private double computeLikelihood(final double[] x) {
        return logLikelihood.apply(x);
    }

    private void recordSample(final double[] x, final ArrayList<double[]> samples) {
        samples.add(Arrays.copyOf(x, x.length));
    }

    final String stateString(final double[] x) {
        return "[" + String.join(", ", Arrays.stream(x).boxed().map(String::valueOf).collect(Collectors.toList())) + "]";
    }

    public ArrayList<double[]> run(final double[] x0, final int numSamples, final int burnIn, final int loggingInterval) {
        double[] x = x0;
        double[] neighborX = new double[size];
        final ArrayList<double[]> samples = new ArrayList<>(numSamples);
        double currentLikelihood = computeLikelihood(x);
        if (loggingInterval > 0) {
            logger.info("Sampling " + x0.length + " parameters over " + numSamples + " steps; f0 = " + currentLikelihood);
        }
        for (int step = 0; step < numSamples + burnIn; step++) {
            if (loggingInterval > 0 && step % loggingInterval == 0) {
                logger.info("Iteration " + step + " : f = " + currentLikelihood);
                logger.info("\tx = " + stateString(x));
            }
            updateNeighbor(x, neighborX);
            final double[] finalNeighborX = neighborX;
            if (!IntStream.range(0, size).anyMatch(i -> finalNeighborX[i] > upperBound[i] || finalNeighborX[i] < lowerBound[i])) {
                final double neighborLikelihood = computeLikelihood(neighborX);
                final double acceptanceProbability = neighborLikelihood >= currentLikelihood ? 1 : Math.exp(neighborLikelihood - currentLikelihood);
                if (acceptanceProbability == 1 || random.nextDouble() <= acceptanceProbability) {
                    final double[] tempX = x;
                    x = neighborX;
                    neighborX = tempX;
                    currentLikelihood = neighborLikelihood;
                }
            }
            if (step >= burnIn) {
                recordSample(x, samples);
            }
        }
        if (loggingInterval > 0) {
            logger.info("Finished after " + numSamples + " : f = " + currentLikelihood);
        }
        return samples;
    }

    public double[] getMAP(final List<double[]> samples) {
        if (samples == null || samples.isEmpty()) {
            throw new IllegalArgumentException("Samples was null or empty");
        }
        double maxLikelihood = computeLikelihood(samples.get(0));
        double[] maxLikelihoodSample = samples.get(0);
        for (int i = 1; i < samples.size(); i++) {
            final double[] sample = samples.get(i);
            final double sampleLikelihood = computeLikelihood(sample);
            if (sampleLikelihood > maxLikelihood) {
                maxLikelihood = sampleLikelihood;
                maxLikelihoodSample = sample;
            }
        }
        return maxLikelihoodSample;
    }

    public double[] getStateProbability(final double state, final List<double[]> samples) {
        if (samples == null) {
            throw new IllegalArgumentException("Samples was null");
        }
        final double[] p = new double[size];
        for (final double[] sample : samples) {
            for (int i = 0; i < size; i++) {
                if (sample[i] == state) {
                    p[i]++;
                }
            }
        }
        final int numSamples = samples.size();
        IntStream.range(0, size).forEach(i -> p[i] = p[i] / numSamples);
        return p;
    }
}
