package org.broadinstitute.hellbender.tools.exome.segmentation;

import org.broadinstitute.hellbender.utils.*;

import java.util.*;
import java.util.function.Function;

/**
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public final class CopyRatioSegmenter extends ClusteringGenomicHMMSegmenter<Double> {
    private double logCoverageCauchyWidth;

    private static final double DEFAULT_INITIAL_CAUCHY_WIDTH = 0.1;
    private static final double NEUTRAL_LOG_2_COPY_RATIO = 0.0;
    private static final double MAX_REASONABLE_CAUCHY_WIDTH = 1.0;
    private static final double MIN_LOG_2_COPY_RATIO = -5.0;
    private static final double MAX_LOG_2_COPY_RATIO = 5.0;

    /**
     * @param initialNumStates  A liberal estimate of the number of hidden minor allele fraction values to
     *                          include in the model.  Hidden states are pruned as the model is learned.
     */
    public CopyRatioSegmenter(final int initialNumStates, final List<SimpleInterval> positions, final List<Double> data) {
        super(initialNumStates, positions, data);
    }

    /**
     * evenly-spaced log-2 copy ratios
     * @param K the initial number of hidden states
     */
    @Override
    protected void initializeHiddenStateValues(final int K) {
        hiddenStateValues = GATKProtectedMathUtils.createEvenlySpacedPoints(-3, 2, K);
        hiddenStateValues[NEUTRAL_VALUE_INDEX] = NEUTRAL_LOG_2_COPY_RATIO;
    }

    @Override
    protected void initializeAdditionalParameters() {
        logCoverageCauchyWidth = DEFAULT_INITIAL_CAUCHY_WIDTH;
    }

    @Override
    protected ClusteringGenomicHMM<Double> makeModel() {
        return new CopyRatioHiddenMarkovModel(hiddenStateValues, weights, memoryLength, logCoverageCauchyWidth);
    }

    @Override
    protected void relearnAdditionalParameters(final ExpectationStep eStep) {
        //relearn the Cauchy width of the emission distribution
        final Function<Double, Double> emissionLogLikelihood = width -> {
            double logLikelihood = 0.0;
            for (int position = 0; position < positions.size(); position++) {
                for (int state = 0; state < weights.length; state++) {
                    final double eStepPosterior = eStep.pStateAtPosition(state, position);
                    logLikelihood += eStepPosterior < NEGLIGIBLE_POSTERIOR_FOR_M_STEP ? 0 :eStepPosterior
                            * CopyRatioHiddenMarkovModel.logEmissionProbability(data.get(position), hiddenStateValues[state], width);
                }
            }
            return logLikelihood;
        };

        logCoverageCauchyWidth = OptimizationUtils.argmax(emissionLogLikelihood, 0, MAX_REASONABLE_CAUCHY_WIDTH, logCoverageCauchyWidth,
                RELATIVE_TOLERANCE_FOR_OPTIMIZATION, ABSOLUTE_TOLERANCE_FOR_OPTIMIZATION, MAX_EVALUATIONS_FOR_OPTIMIZATION);

        logger.info("New coverage standard deviation learned: " + logCoverageCauchyWidth);
    }

    @Override
    protected double minHiddenStateValue() { return MIN_LOG_2_COPY_RATIO; }

    @Override
    protected double maxHiddenStateValue() { return  MAX_LOG_2_COPY_RATIO; }
}
