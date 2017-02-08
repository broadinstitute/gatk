package org.broadinstitute.hellbender.tools.exome.segmentation;

import org.broadinstitute.hellbender.utils.OptimizationUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.List;
import java.util.function.Function;
import java.util.stream.IntStream;

/**
 * HMM segmenter for scalar hidden data, such as allele fraction and copy ratio, but not the joint
 * allele fraction / copy ratio segmenter
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public abstract class ScalarHMMSegmenter<DATA> extends ClusteringGenomicHMMSegmenter<DATA, Double> {

    private static final double DEFAULT_MEMORY_LENGTH = 5e7;
    private static final double EPSILON = 1E-10;

    public ScalarHMMSegmenter(final List<SimpleInterval> positions, final List<DATA> data, final List<Double> initialNonConstantHiddenStates) {
        super(positions, data, initialNonConstantHiddenStates, DEFAULT_MEMORY_LENGTH);
    }

    @Override
    protected boolean hiddenStateValuesHaveConverged(final List<Double> oldHiddenStateValues) {
        return oldHiddenStateValues.size() == numStates() && maxRelativeDifference(oldHiddenStateValues, getStates()) < RELATIVE_CONVERGENCE_THRESHOLD;
    }

    @Override
    protected void relearnHiddenStateValues(final ExpectationStep eStep) {
        final ClusteringGenomicHMM<DATA, Double> model = makeModel();
        IntStream.range(0, numStates()).forEach(state -> {
            final Function<Double, Double> objective = f -> IntStream.range(0, data.size())
                    .filter(n -> eStep.pStateAtPosition(state, n) > NEGLIGIBLE_POSTERIOR_FOR_M_STEP)
                    .mapToDouble(n -> eStep.pStateAtPosition(state, n) * model.logEmissionProbability(data.get(n), f))
                    .sum();
            setState(state, OptimizationUtils.singleNewtonArgmaxUpdate(objective, minHiddenStateValue(),
                    maxHiddenStateValue(), getState(state)));
        });
    }

    protected abstract double minHiddenStateValue();
    protected abstract double maxHiddenStateValue();

    private static double maxRelativeDifference(final List<Double> array1, final List<Double> array2) {
        Utils.validateArg(array1.size() == array2.size(), "arrays must have same length.");
        Utils.validateArg(array1.size() > 0, "arrays must be non-empty");
        return IntStream.range(0, array1.size()).mapToDouble(n -> Math.abs((array1.get(n) - array2.get(n)) / (array1.get(n) + EPSILON))).max().getAsDouble();
    }
}
