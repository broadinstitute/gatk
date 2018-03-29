package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;
import org.broadinstitute.hellbender.utils.hmm.HMM;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public final class CopyRatioHMM implements HMM<Double, Integer, Integer> {

    private final RealVector logPriors;
    private final int numStates;
    private final double logStateChangeProbability;
    private final double logStateSameProbability;
    private final List<NormalDistribution> copyRatioEmissionDistributions;

    private static final double COPY_RATIO_STDEV = 0.25;

    public CopyRatioHMM(final RealVector priors, final double stateChangeProbability) {
        this.numStates = priors.getDimension();
        this.logPriors = priors.map(Math::log);
        this.logStateChangeProbability = Math.log(stateChangeProbability);
        this.logStateSameProbability = Math.log(1 - stateChangeProbability);
        copyRatioEmissionDistributions = new ArrayList<>(numStates);
        for (int i = 0; i < numStates; i++) {
            copyRatioEmissionDistributions.add(new NormalDistribution(i * 0.5, COPY_RATIO_STDEV));
        }
    }

    @Override
    public List<Integer> hiddenStates() {
        return IntStream.range(0, numStates).boxed().collect(Collectors.toList());
    }

    @Override
    public double logPriorProbability(final Integer state, final Integer position) {
        return logPriors.getEntry(state);
    }

    @Override
    public double logTransitionProbability(final Integer currentState,
                                           final Integer currentPosition,
                                           final Integer nextState,
                                           final Integer nextPosition) {
        return currentState == nextState ? logStateSameProbability : logStateChangeProbability;
    }

    @Override
    public double logEmissionProbability(final Double copyRatio, final Integer state,
                                         final Integer position) {
        return copyRatioEmissionDistributions.get(state).logDensity(copyRatio);
    }

    public static RealVector getCopyNumberPrior(final int numStates) {
        final RealVector prior = new ArrayRealVector(new double[numStates]);
        prior.set(1.0 / numStates);
        return prior;
    }

    public static List<Integer> getCopyNumberHMMPositions(final int length) {
        return IntStream.range(0, length).boxed().collect(Collectors.toList());
    }

}
