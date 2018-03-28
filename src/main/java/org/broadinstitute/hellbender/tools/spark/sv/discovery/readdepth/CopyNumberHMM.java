package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;
import org.broadinstitute.hellbender.utils.hmm.HMM;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public final class CopyNumberHMM implements HMM<Double, Integer, Integer> {

    private final RealVector priors;
    private final int numStates;
    private final double alpha;
    private final List<NormalDistribution> copyNumberRatioEmissionDistributions;

    public CopyNumberHMM(final RealVector priors, final double alpha) {
        this.priors = priors;
        numStates = priors.getDimension();
        this.alpha = alpha;
        copyNumberRatioEmissionDistributions = new ArrayList<>(numStates);
        for (int i = 0; i < numStates; i++) {
            copyNumberRatioEmissionDistributions.add(new NormalDistribution(i * 0.5, 0.25));
        }
    }

    @Override
    public List<Integer> hiddenStates() {
        return IntStream.range(0, numStates).boxed().collect(Collectors.toList());
    }

    @Override
    public double logPriorProbability(final Integer state, final Integer position) {
        return Math.log(priors.getEntry(state));
    }

    @Override
    public double logTransitionProbability(final Integer currentState,
                                           final Integer currentPosition,
                                           final Integer nextState,
                                           final Integer nextPosition) {
        return currentState == nextState ? Math.log(1 - alpha) : Math.log(alpha);
    }

    @Override
    public double logEmissionProbability(final Double log2CopyRatio, final Integer state,
                                         final Integer position) {
        return copyNumberRatioEmissionDistributions.get(state).logDensity(Math.pow(2, log2CopyRatio));
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
