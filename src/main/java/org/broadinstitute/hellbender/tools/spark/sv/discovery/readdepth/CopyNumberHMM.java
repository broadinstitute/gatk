package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import com.netflix.servo.util.VisibleForTesting;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.hmm.HMM;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Simple Hidden Markov Model for predicting copy number from copy ratios. Each copy number has a Gaussian emission
 * distribution with means in multiples of 0.5 (0:0, 1:0.5, 2:1, etc.) and fixed variance. The probability of transitioning
 * between any two different hidden states is a constant.
 */
public final class CopyNumberHMM implements HMM<Double, Integer, Integer> {

    @VisibleForTesting
    static final double COPY_RATIO_STDEV = 0.25;

    private static final double PRIOR_CHECK_TOLERANCE = 1e-6;
    private final RealVector logPriors;
    private final int numStates;
    private final double logStateTransitionProbability;
    private final double logNoStateTransitionProbability;
    private final List<NormalDistribution> copyRatioEmissionDistributions;

    public CopyNumberHMM(final RealVector priors, final double totalStateTransitionProbability) {
        double priorSum = 0;
        for (final double val : priors.toArray()) {
            Utils.validateArg(val >= 0 && val <= 1, "Prior values must be probabilities on the interval [0,1]");
            priorSum += val;
        }
        Utils.validateArg(priorSum >= 1 - PRIOR_CHECK_TOLERANCE && priorSum <= 1 + PRIOR_CHECK_TOLERANCE, "Prior vector must sum to 1");
        Utils.validateArg(totalStateTransitionProbability >= 0 && totalStateTransitionProbability <= 1, "State change probability must be on the interval [0,1]");
        this.numStates = priors.getDimension();
        this.logPriors = priors.map(Math::log);
        this.logStateTransitionProbability = Math.log(totalStateTransitionProbability / (numStates - 1));
        this.logNoStateTransitionProbability = Math.log(1 - totalStateTransitionProbability);
        copyRatioEmissionDistributions = new ArrayList<>(numStates);
        for (int i = 0; i < numStates; i++) {
            copyRatioEmissionDistributions.add(new NormalDistribution(i * 0.5, COPY_RATIO_STDEV));
        }
    }

    public static RealVector uniformPrior(final int numStates) {
        Utils.validateArg(numStates >= 0, "Number of states must be non-negative");
        final RealVector prior = new ArrayRealVector(new double[numStates]);
        prior.set(1.0 / numStates);
        return prior;
    }

    public static List<Integer> positionsList(final int length) {
        Utils.validateArg(length >= 0, "Length must be non-negative");
        return IntStream.range(0, length).boxed().collect(Collectors.toList());
    }

    private void checkState(final int state) {
        if (state < 0 || state >= numStates) throw new IllegalArgumentException("Invalid state " + state);
    }

    @Override
    public List<Integer> hiddenStates() {
        return IntStream.range(0, numStates).boxed().collect(Collectors.toList());
    }

    @Override
    public double logPriorProbability(final Integer state, final Integer position) {
        checkState(state);
        return logPriors.getEntry(state);
    }

    @Override
    public double logTransitionProbability(final Integer currentState,
                                           final Integer currentPosition,
                                           final Integer nextState,
                                           final Integer nextPosition) {
        checkState(currentState);
        checkState(nextState);
        return currentState == nextState ? logNoStateTransitionProbability : logStateTransitionProbability;
    }

    @Override
    public double logEmissionProbability(final Double copyRatio, final Integer state,
                                         final Integer position) {
        checkState(state);
        return copyRatioEmissionDistributions.get(state).logDensity(copyRatio);
    }

}
