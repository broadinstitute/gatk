package org.broadinstitute.hellbender.tools.exome.segmentation;

import org.apache.commons.collections4.ListUtils;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.OptimizationUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * HMM segmenter for scalar hidden data, such as allele fraction and copy ratio, but not the joint
 * allele fraction / copy ratio segmenter
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public abstract class ScalarHMMSegmenter<DATA> extends ClusteringGenomicHMMSegmenter<DATA, Double> {

    private static final double DEFAULT_INITIAL_CONCENTRATION = 1;

    private final int numConstantStates;
    private static final double DEFAULT_MEMORY_LENGTH = 5e7;

    public ScalarHMMSegmenter(final List<SimpleInterval> positions, final List<DATA> data,
                              final List<Double> constantHiddenStates, final List<Double> initialNonConstantHiddenStates) {
        super(positions, data, ListUtils.union(constantHiddenStates, initialNonConstantHiddenStates),
                uniformWeights(constantHiddenStates.size() + initialNonConstantHiddenStates.size()),
                DEFAULT_INITIAL_CONCENTRATION, DEFAULT_MEMORY_LENGTH);
        numConstantStates = constantHiddenStates.size();
    }

    @Override
    protected boolean hiddenStateValuesHaveConverged(final List<Double> oldHiddenStateValues) {
        return oldHiddenStateValues.size() == numStates() && GATKProtectedMathUtils.maxDifference(oldHiddenStateValues, getStates()) < CONVERGENCE_THRESHOLD;
    }

    // filter out components that have low weight and are too close to another component -- these will
    // die out eventually in EM, but very slowly, so we hasten their demise for quicker convergence
    @Override
    protected void pruneUnusedComponents() {
        final Set<Integer> componentsToPrune = new TreeSet<>();
        for (final int state : nonConstantStateIndices()) {
            final int closestOtherState = MathUtils.maxElementIndex(IntStream.range(0, numStates())
                    .mapToDouble(n -> n == state ? Double.NEGATIVE_INFINITY : -Math.abs(getState(n) - getState(state)))
                    .toArray());
            final boolean hasLowWeight = getWeight(state) < MAX_WEIGHT_CONSIDERED_FOR_PRUNING;
            final boolean isCloseToNeighbor = Math.abs(getState(state) - getState(closestOtherState)) < DISTANCE_TO_NEIGHBOR_TO_BE_CONSIDERED_SPURIOUS;
            final boolean isVeryCloseToNeighbor = Math.abs(getState(state) - getState(closestOtherState)) < DISTANCE_TO_NEIGHBOR_TO_BE_CONSIDERED_DEFINITELY_SPURIOUS;
            final boolean hasLessWeightThanNeighbor = getWeight(state) < getWeight(closestOtherState);
            final boolean isTinyWeight = getWeight(state) < AUTOMATICALLY_PRUNED_WEIGHT;
            final boolean isParasite = hasLowWeight && isCloseToNeighbor && hasLessWeightThanNeighbor;
            final boolean isClone = isVeryCloseToNeighbor && hasLessWeightThanNeighbor;
            if ( isTinyWeight || isParasite || isClone) {
                componentsToPrune.add(state);
            }
        }

        removeStates(componentsToPrune);
    }

    private static List<Double> uniformWeights(final int numStates) {
        return Collections.nCopies(numStates, 1.0/numStates);
    }

    @Override
    protected void relearnHiddenStateValues(final ExpectationStep eStep) {
        final ClusteringGenomicHMM<DATA, Double> model = makeModel();
        for (final int state : nonConstantStateIndices()) {
            final Function<Double, Double> objective = f -> IntStream.range(0, data.size())
                    .filter(n -> eStep.pStateAtPosition(state, n) > NEGLIGIBLE_POSTERIOR_FOR_M_STEP)
                    .mapToDouble(n -> eStep.pStateAtPosition(state, n) * model.logEmissionProbability(data.get(n), f))
                    .sum();
            setState(state, OptimizationUtils.singleNewtonArgmaxUpdate(objective, minHiddenStateValue(),
                    maxHiddenStateValue(), getState(state)));
        }
    }

    protected abstract double minHiddenStateValue();
    protected abstract double maxHiddenStateValue();

    // constant states are always first -- see constructor
    protected List<Integer> nonConstantStateIndices() { return IntStream.range(numConstantStates, numStates()).boxed().collect(Collectors.toList());}
}
