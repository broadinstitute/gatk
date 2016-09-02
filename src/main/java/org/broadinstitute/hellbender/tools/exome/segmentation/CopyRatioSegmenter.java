package org.broadinstitute.hellbender.tools.exome.segmentation;

import com.google.common.primitives.Doubles;
import org.apache.commons.collections4.ListUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.tools.exome.HashedListTargetCollection;
import org.broadinstitute.hellbender.tools.exome.ModeledSegment;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.TargetCollection;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.OptimizationUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.Arrays;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public final class CopyRatioSegmenter extends ScalarHMMSegmenter<Double> {
    private double logCoverageCauchyWidth;

    private static final double DEFAULT_INITIAL_CAUCHY_WIDTH = 0.1;
    private static final double MAX_REASONABLE_CAUCHY_WIDTH = 1.0;
    private static final double MIN_LOG_2_COPY_RATIO = -5.0;
    private static final double MAX_LOG_2_COPY_RATIO = 5.0;
    private static final List<Double> CONSTANT_LOG_2_COPY_RATIO_OF_1 = Arrays.asList(0.0);
    private static final double MIN_INITIAL_LOG_2_COPY_RATIO = -3.0;
    private static final double MAX_INITIAL_LOG_2_COPY_RATIO = 2.0;

    /**
     * @param initialNumStates  A liberal estimate of the number of hidden minor allele fraction values to
     *                          include in the model.  Hidden states are pruned as the model is learned.
     */
    public CopyRatioSegmenter(final int initialNumStates, final ReadCountCollection rcc) {
        super(rcc.targets().stream().map(SimpleInterval::new).collect(Collectors.toList()), Doubles.asList(rcc.getColumn(0)),
                CONSTANT_LOG_2_COPY_RATIO_OF_1, initialNonConstantLog2CopyRatios(initialNumStates - 1));
        Utils.validateArg(rcc.columnNames().size() == 1, "Only single-sample ReadCountCollection is supported.");
        logCoverageCauchyWidth = DEFAULT_INITIAL_CAUCHY_WIDTH;
    }

    /**
     * evenly-spaced log-2 copy ratios
     * @param K the initial number of hidden states
     */
    private static List<Double> initialNonConstantLog2CopyRatios(final int K) {
        ParamUtils.isPositive(K, "must have at least one non-constant state");
        final double spacing = (MAX_INITIAL_LOG_2_COPY_RATIO  - MIN_INITIAL_LOG_2_COPY_RATIO) / (K + 1);
        final int numNegativeStates = K / 2;
        final int numPositiveStates = K - numNegativeStates;
        final List<Double> negativeStates = Doubles.asList(GATKProtectedMathUtils.createEvenlySpacedPoints(MIN_INITIAL_LOG_2_COPY_RATIO, spacing, numNegativeStates));
        final List<Double> positiveStates = Doubles.asList(GATKProtectedMathUtils.createEvenlySpacedPoints(spacing, MAX_INITIAL_LOG_2_COPY_RATIO, numPositiveStates));
        return ListUtils.union(negativeStates, positiveStates);
    }

    public List<ModeledSegment> getModeledSegments() {
        final List<Pair<SimpleInterval, Double>> segmentation = findSegments();
        final TargetCollection<SimpleInterval> tc = new HashedListTargetCollection<>(positions);
        return segmentation.stream().map(pair ->
                new ModeledSegment(pair.getLeft(), tc.targetCount(pair.getLeft()), pair.getRight())).collect(Collectors.toList());
    }

    @Override
    protected ClusteringGenomicHMM<Double, Double> makeModel() {
        return new CopyRatioHiddenMarkovModel(getStates(), getWeights(), getMemoryLength(), logCoverageCauchyWidth);
    }

    @Override
    protected void relearnAdditionalParameters(final ExpectationStep eStep) {
        //relearn the Cauchy width of the emission distribution
        final Function<Double, Double> emissionLogLikelihood = width -> {
            double logLikelihood = 0.0;
            for (int position = 0; position < numPositions(); position++) {
                for (int state = 0; state < numStates(); state++) {
                    final double eStepPosterior = eStep.pStateAtPosition(state, position);
                    logLikelihood += eStepPosterior < NEGLIGIBLE_POSTERIOR_FOR_M_STEP ? 0 : eStepPosterior
                            * CopyRatioHiddenMarkovModel.logEmissionProbability(data.get(position), getState(state), width);
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

    public double getLogCoverageCauchyWidth() {
        return logCoverageCauchyWidth;
    }

}
