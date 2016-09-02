package org.broadinstitute.hellbender.tools.exome.segmentation;

import com.google.cloud.dataflow.sdk.repackaged.com.google.common.primitives.Doubles;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.tools.exome.HashedListTargetCollection;
import org.broadinstitute.hellbender.tools.exome.ModeledSegment;
import org.broadinstitute.hellbender.tools.exome.TargetCollection;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionGlobalParameters;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionInitializer;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionState;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.pon.allelic.AllelicPanelOfNormals;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.OptimizationUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public final class AlleleFractionSegmenter extends ScalarHMMSegmenter<AllelicCount> {
    private final AllelicPanelOfNormals allelicPoN;
    private AlleleFractionGlobalParameters globalParameters;
    private static final double INITIAL_MEAN_ALLELIC_BIAS = 1.0;
    private static final double INITIAL_ALLELIC_BIAS_VARIANCE = 1e-2;
    private static final double INITIAL_OUTLIER_PROBABILITY = 3e-2;
    private static final List<Double> NEUTRAL_MINOR_ALLELE_FRACTION = Arrays.asList(0.5);

    /**
     * Initialize the segmenter with its data and panel of normals, giving equal weight to a set of evenly-spaced
     * hidden minor allele fraction values.
     *
     * @param initialNumStates  A liberal estimate of the number of hidden minor allele fraction values to
     *                          include in the model.  Hidden states are pruned as the model is learned.
     * @param acc               The {@link AllelicCountCollection} data attached to this segmenter
     * @param allelicPoN        The {@link AllelicPanelOfNormals} attached to this segmenter
     */
    public AlleleFractionSegmenter(final int initialNumStates, final AllelicCountCollection acc,
                                   final AllelicPanelOfNormals allelicPoN) {
        super(acc.getCounts().stream().map(AllelicCount::getInterval).collect(Collectors.toList()),
                acc.getCounts(), NEUTRAL_MINOR_ALLELE_FRACTION, initialNonConstantMinorFractions(initialNumStates - 1));
        this.allelicPoN = Utils.nonNull(allelicPoN);
        globalParameters = new AlleleFractionGlobalParameters(INITIAL_MEAN_ALLELIC_BIAS, INITIAL_ALLELIC_BIAS_VARIANCE, INITIAL_OUTLIER_PROBABILITY);
    }

    /**
     * evenly-spaced minor allele fractions going from 0 (inclusive) to 1/2 (exclusive)
     * @param K the initial number of hidden states
     */
    private static List<Double> initialNonConstantMinorFractions(final int K) {
        ParamUtils.isPositive(K, "must have at least one non-constant state");
        final double spacing = 0.5 / (K + 1);
        return Doubles.asList(GATKProtectedMathUtils.createEvenlySpacedPoints(0.001, 0.5 - spacing, K));
    }

    public List<ModeledSegment> getModeledSegments() {
        final TargetCollection<SimpleInterval> tc = new HashedListTargetCollection<>(positions);
        final List<Pair<SimpleInterval, Double>> segmentation = findSegments();
        return segmentation.stream()
                .map(pair -> new ModeledSegment(pair.getLeft(), tc.targetCount(pair.getLeft()), pair.getRight()))
                .collect(Collectors.toList());
    }

    @Override
    protected ClusteringGenomicHMM<AllelicCount, Double> makeModel() {
        return new AlleleFractionHiddenMarkovModel(getStates(), getWeights(), getMemoryLength(), allelicPoN, globalParameters);
    }

    @Override
    protected void relearnAdditionalParameters(final ExpectationStep eStep) {
        final Function<AlleleFractionGlobalParameters, Double> emissionLogLikelihood = params -> {
            double logLikelihood = 0.0;
            for (int position = 0; position < numPositions(); position++) {
                for (int state = 0; state < numStates(); state++) {
                    final double eStepPosterior = eStep.pStateAtPosition(state, position);
                    logLikelihood += eStepPosterior < NEGLIGIBLE_POSTERIOR_FOR_M_STEP ? 0 :eStepPosterior
                            * AlleleFractionHiddenMarkovModel.logEmissionProbability(data.get(position), getState(state), params, allelicPoN);
                }
            }
            return logLikelihood;
        };

        final Function<Double, Double> meanBiasObjective = mean -> emissionLogLikelihood.apply(globalParameters.copyWithNewMeanBias(mean));
        final double newMeanBias = OptimizationUtils.argmax(meanBiasObjective, 0, AlleleFractionInitializer.MAX_REASONABLE_MEAN_BIAS, globalParameters.getMeanBias(),
                RELATIVE_TOLERANCE_FOR_OPTIMIZATION, ABSOLUTE_TOLERANCE_FOR_OPTIMIZATION, MAX_EVALUATIONS_FOR_OPTIMIZATION);

        final Function<Double, Double> biasVarianceObjective = variance -> emissionLogLikelihood.apply(globalParameters.copyWithNewBiasVariance(variance));
        final double newBiasVariance = OptimizationUtils.argmax(biasVarianceObjective, 0, AlleleFractionInitializer.MAX_REASONABLE_BIAS_VARIANCE, globalParameters.getBiasVariance(),
                RELATIVE_TOLERANCE_FOR_OPTIMIZATION, ABSOLUTE_TOLERANCE_FOR_OPTIMIZATION, MAX_EVALUATIONS_FOR_OPTIMIZATION);

        final Function<Double, Double> outlierProbabilityObjective = pOutlier -> emissionLogLikelihood.apply(globalParameters.copyWithNewOutlierProbability(pOutlier));
        final double newOutlierProbability = OptimizationUtils.argmax(outlierProbabilityObjective, 0, AlleleFractionInitializer.MAX_REASONABLE_OUTLIER_PROBABILITY, globalParameters.getOutlierProbability(),
                RELATIVE_TOLERANCE_FOR_OPTIMIZATION, ABSOLUTE_TOLERANCE_FOR_OPTIMIZATION, MAX_EVALUATIONS_FOR_OPTIMIZATION);

        globalParameters = new AlleleFractionGlobalParameters(newMeanBias, newBiasVariance, newOutlierProbability);

        logger.info(String.format("Global allelic bias parameters learned.  Mean allelic bias: %f, variance of allelic bias: %f, outlier probability: %f.",
                newMeanBias, newBiasVariance, newOutlierProbability));
    }

    @Override
    protected double minHiddenStateValue() { return AlleleFractionState.MIN_MINOR_FRACTION; }

    @Override
    protected double maxHiddenStateValue() { return  AlleleFractionState.MAX_MINOR_FRACTION; }

    public AllelicPanelOfNormals getAllelicPoN() {
        return allelicPoN;
    }

    public AlleleFractionGlobalParameters getGlobalParameters() {
        return globalParameters;
    }
}
