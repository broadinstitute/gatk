package org.broadinstitute.hellbender.utils.solver;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Arrays;
import java.util.Random;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class SimulatedAnnealingSolver {
    private final int size;
    private final Random random = new Random();
    private final Function<double[], Double> energyFunction;
    private final Supplier<double[]> stepSampler;
    private final double[] lowerBound;
    private final double[] upperBound;
    private double[] lastSolution;
    private final Logger logger = LogManager.getLogger(this.getClass());

    public SimulatedAnnealingSolver(final int size,
                                    final Function<double[], Double> energyFunction,
                                    final Supplier<double[]> stepSampler,
                                    final double[] lowerBound,
                                    final double[] upperBound) {
        Utils.nonNull(energyFunction, "Energy function cannot be null");
        Utils.nonNull(stepSampler, "Step sampler cannot be null");
        Utils.nonNull(lowerBound, "Lower bound array cannot be null");
        Utils.nonNull(upperBound, "Upper bound array cannot be null");
        Utils.validateArg(lowerBound.length == size, "Lower bound array must have length equal to size parameter");
        Utils.validateArg(upperBound.length == size, "Upper bound array must have length equal to size parameter");
        this.size = size;
        this.energyFunction = energyFunction;
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
        IntStream.range(0, size).forEach(i -> neighborX[i] = Math.max(lowerBound[i], Math.min(upperBound[i], x[i] + steps[i])));
    }

    private double computeEnergy(final double[] x) {
        return energyFunction.apply(x);
    }

    private void copyStates(final double[] fromX, final double[] toX) {
        IntStream.range(0, size).forEach(i -> toX[i] = fromX[i]);
    }

    public double[] getSolution() {
        return lastSolution;
    }

    public double solve(final double[] x0, final int numSteps, final int loggingInterval) {
        final double T0 = 2.0;
        final Function<Integer,Double> temperatureSchedule = step -> T0 * (1 - (step / (double) numSteps));
        return solve(x0, numSteps, temperatureSchedule, loggingInterval);
    }

    final String stateString(final double[] x) {
        return "[" + String.join(", ", Arrays.stream(x).boxed().map(String::valueOf).collect(Collectors.toList())) + "]";
    }

    public double solve(final double[] x0, final int numSteps, final Function<Integer,Double> temperatureSchedule, final int loggingInterval) {
        double[] x = x0;
        double[] neighborX = new double[size];
        final double[] minX = new double[size];
        double currentEnergy = computeEnergy(x);
        if (loggingInterval > 0) {
            logger.info("Solving " + x0.length + " parameters over " + numSteps + " steps; f0 = " + currentEnergy);
        }
        double minEnergy = currentEnergy;
        int tMin = 0;
        copyStates(x, minX);
        for (int step = 0; step < numSteps; step++) {
            final double T = temperatureSchedule.apply(step);
            if (loggingInterval > 0 && step % loggingInterval == 0) {
                logger.info("Iteration " + step + " : T = " + T + ", f = " + currentEnergy + ", fmin = " + minEnergy);
                logger.info("\tx = " + stateString(x));
            }
            updateNeighbor(x, neighborX);
            final double neighborEnergy = computeEnergy(neighborX);
            final double acceptanceProbability;
            if (neighborEnergy <= currentEnergy) {
                acceptanceProbability = 1;
            } else {
                acceptanceProbability = Math.exp((currentEnergy - neighborEnergy) / T);
            }
            if (acceptanceProbability == 1 || random.nextDouble() <= acceptanceProbability) {
                final double[] tempX = x;
                x = neighborX;
                neighborX = tempX;
                currentEnergy = neighborEnergy;
                if (currentEnergy < minEnergy) {
                    minEnergy = currentEnergy;
                    copyStates(x, minX);
                    tMin = step;
                }
            }
        }
        if (loggingInterval > 0) {
            logger.info("Solved after " + numSteps + " : f = " + currentEnergy + ", fmin = " + minEnergy + ", tmin = " + tMin);
            logger.info("\tx_min = " + stateString(minX));
        }
        lastSolution = minX;
        return minEnergy;
    }
}
