package org.broadinstitute.hellbender.tools.exome.segmentation;

import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionInitializer;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionState;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionGlobalParameters;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AllelicPanelOfNormals;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.utils.OptimizationUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public final class AlleleFractionSegmenter extends ClusteringGenomicHMMSegmenter<AllelicCount> {
    private final AllelicPanelOfNormals allelicPoN;
    private AlleleFractionGlobalParameters biasParameters;
    private static final double INITIAL_MEAN_ALLELIC_BIAS = 1.0;
    private static final double INITIAL_ALLELIC_BIAS_VARIANCE = 1e-2;
    private static final double INITIAL_OUTLIER_PROBABILITY = 3e-2;


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
        super(initialNumStates, acc.getCounts().stream().map(AllelicCount::getInterval).collect(Collectors.toList()), acc.getCounts());
        this.allelicPoN = Utils.nonNull(allelicPoN);
    }

    /**
     * evenly-spaced minor allele fractions going from 1/2 to 0
     * @param K the initial number of hidden states
     */
    @Override
    protected void initializeHiddenStateValues(final int K) {
        hiddenStateValues = IntStream.range(0, K).mapToDouble(n ->  ((double) K - n) / (2*K)).toArray();
    }

    @Override
    protected void initializeAdditionalParameters() {
        biasParameters = new AlleleFractionGlobalParameters(INITIAL_MEAN_ALLELIC_BIAS, INITIAL_ALLELIC_BIAS_VARIANCE, INITIAL_OUTLIER_PROBABILITY);
    }

    @Override
    protected ClusteringGenomicHMM<AllelicCount> makeModel() {
        return new AlleleFractionHiddenMarkovModel(hiddenStateValues, weights, memoryLength, allelicPoN, biasParameters);
    }

    @Override
    protected void relearnAdditionalParameters(final ExpectationStep eStep) {
        final Function<AlleleFractionGlobalParameters, Double> emissionLogLikelihood = params -> {
            double logLikelihood = 0.0;
            for (int position = 0; position < positions.size(); position++) {
                for (int state = 0; state < weights.length; state++) {
                    final double eStepPosterior = eStep.pStateAtPosition(state, position);
                    logLikelihood += eStepPosterior < NEGLIGIBLE_POSTERIOR_FOR_M_STEP ? 0 :eStepPosterior
                            * AlleleFractionHiddenMarkovModel.logEmissionProbability(data.get(position), hiddenStateValues[state], params, allelicPoN);
                }
            }
            return logLikelihood;
        };

        final Function<Double, Double> meanBiasObjective = mean -> emissionLogLikelihood.apply(biasParameters.copyWithNewMeanBias(mean));
        final double newMeanBias = OptimizationUtils.argmax(meanBiasObjective, 0, AlleleFractionInitializer.MAX_REASONABLE_MEAN_BIAS, biasParameters.getMeanBias(),
                RELATIVE_TOLERANCE_FOR_OPTIMIZATION, ABSOLUTE_TOLERANCE_FOR_OPTIMIZATION, MAX_EVALUATIONS_FOR_OPTIMIZATION);

        final Function<Double, Double> biasVarianceObjective = variance -> emissionLogLikelihood.apply(biasParameters.copyWithNewBiasVariance(variance));
        final double newBiasVariance = OptimizationUtils.argmax(biasVarianceObjective, 0, AlleleFractionInitializer.MAX_REASONABLE_BIAS_VARIANCE, biasParameters.getBiasVariance(),
                RELATIVE_TOLERANCE_FOR_OPTIMIZATION, ABSOLUTE_TOLERANCE_FOR_OPTIMIZATION, MAX_EVALUATIONS_FOR_OPTIMIZATION);

        final Function<Double, Double> outlierProbabilityObjective = pOutlier -> emissionLogLikelihood.apply(biasParameters.copyWithNewOutlierProbability(pOutlier));
        final double newOutlierProbability = OptimizationUtils.argmax(outlierProbabilityObjective, 0, AlleleFractionInitializer.MAX_REASONABLE_OUTLIER_PROBABILITY, biasParameters.getOutlierProbability(),
                RELATIVE_TOLERANCE_FOR_OPTIMIZATION, ABSOLUTE_TOLERANCE_FOR_OPTIMIZATION, MAX_EVALUATIONS_FOR_OPTIMIZATION);

        biasParameters = new AlleleFractionGlobalParameters(newMeanBias, newBiasVariance, newOutlierProbability);

        logger.info(String.format("Global allelic bias parameters learned.  Mean allelic bias: %f, variance of allelic bias: %f, outlier probability: %f.",
                newMeanBias, newBiasVariance, newOutlierProbability));
    }

    @Override
    protected double minHiddenStateValue() { return AlleleFractionState.MIN_MINOR_FRACTION; }

    @Override
    protected double maxHiddenStateValue() { return  AlleleFractionState.MAX_MINOR_FRACTION; }
}
