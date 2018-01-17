package org.broadinstitute.hellbender.tools.copynumber.models;

import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.special.Beta;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.copynumber.formats.CopyNumberFormatsUtils;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.AllelicCount;
import org.broadinstitute.hellbender.utils.OptimizationUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * The allele-fraction model (after marginalizing latent parameters as described in docs/CNVs/CNV-methods.pdf)
 * contains the following parameters:
 *      1.  minor-allele fractions for each segment
 *      2.  a global outlier probability
 *      3.  the mean allelic bias
 *      4.  the rate (mean / variance) of the allelic bias
 *
 * Note that 3 and 4 are hyperparameters specifying a gamma distribution prior on allelic bias -- the latent variables
 * for bias at each het site have been marginalized but the hyperparameters have not.
 *
 * The allele-fraction model samples the distribution of these parameters using Markov chain Monte Carlo and in principle
 * an initialization step is not necessary.  However, in practice this initialization finds the mode of the posterior
 * distributions in only a few iterations, whereas sampling would require many more.  Thus we greatly reduce the
 * number of burn-in samples that we must discard.
 *
 * The initialization is straightforward: first we set the minor fractions to reasonable guesses based on alt and ref
 * counts, assuming no allelic bias.  Then we numerically maximize the likelihood with respect to each parameter until
 * the likelihood converges to a maximum.  In practice this is the unique global maximum.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
final class AlleleFractionInitializer {
    private static final Logger logger = LogManager.getLogger(AlleleFractionInitializer.class);

    private static final double INITIAL_OUTLIER_PROBABILITY = 0.01;
    private static final double INITIAL_MEAN_BIAS = 1.0;
    private static final double INITIAL_BIAS_VARIANCE = 0.05;    //this is an overestimate, but starting small makes it slow for
                                                                 //mean bias to escape a bad initial guess
    private static final AlleleFractionGlobalParameters INITIAL_GLOBAL_PARAMETERS =
        new AlleleFractionGlobalParameters(INITIAL_MEAN_BIAS, INITIAL_BIAS_VARIANCE, INITIAL_OUTLIER_PROBABILITY);

    private static final double LOG_LIKELIHOOD_CONVERGENCE_THRESHOLD = 0.5;
    private static final int MAX_ITERATIONS = 50;

    //define maxima of search intervals for maximum likelihood -- parameter values above these would be ridiculous
    static final double MAX_REASONABLE_OUTLIER_PROBABILITY = 0.15;
    static final double MAX_REASONABLE_MEAN_BIAS = 5.0;
    static final double MAX_REASONABLE_BIAS_VARIANCE = 0.5;
    private static final double EPSILON_FOR_NEAR_MAX_WARNING = 1E-2;

    //the minor-allele fraction of a segment must be less than one half by definition
    private static final double MAX_MINOR_ALLELE_FRACTION = 0.5;

    private final AlleleFractionSegmentedData data;
    private AlleleFractionGlobalParameters globalParameters;
    private AlleleFractionState.MinorFractions minorFractions;

    /**
     * This constructor performs the initialization.
     */
    AlleleFractionInitializer(final AlleleFractionSegmentedData data) {
        this.data = Utils.nonNull(data);
        globalParameters = INITIAL_GLOBAL_PARAMETERS;
        minorFractions = calculateInitialMinorFractions(data);
        double previousIterationLogLikelihood;
        double nextIterationLogLikelihood = Double.NEGATIVE_INFINITY;
        logger.info(String.format("Initializing allele-fraction model, iterating until log likelihood converges to within %s...",
                CopyNumberFormatsUtils.formatDouble(LOG_LIKELIHOOD_CONVERGENCE_THRESHOLD)));
        int iteration = 1;
        do {
            previousIterationLogLikelihood = nextIterationLogLikelihood;
            globalParameters = new AlleleFractionGlobalParameters(
                    estimateMeanBias(), estimateBiasVariance(), estimateOutlierProbability());
            minorFractions = estimateMinorFractions();

            nextIterationLogLikelihood = AlleleFractionLikelihoods.logLikelihood(globalParameters, minorFractions, data);
            logger.info(String.format("Iteration %d, model log likelihood = %s...", iteration,
                    CopyNumberFormatsUtils.formatDouble(nextIterationLogLikelihood)));
            logger.info(globalParameters);
            iteration++;
        } while (iteration < MAX_ITERATIONS &&
                nextIterationLogLikelihood - previousIterationLogLikelihood > LOG_LIKELIHOOD_CONVERGENCE_THRESHOLD);
        warnIfNearMax(AlleleFractionParameter.MEAN_BIAS.name, globalParameters.getMeanBias(), MAX_REASONABLE_MEAN_BIAS, EPSILON_FOR_NEAR_MAX_WARNING);
        warnIfNearMax(AlleleFractionParameter.BIAS_VARIANCE.name, globalParameters.getBiasVariance(), MAX_REASONABLE_BIAS_VARIANCE, EPSILON_FOR_NEAR_MAX_WARNING);
        warnIfNearMax(AlleleFractionParameter.OUTLIER_PROBABILITY.name, globalParameters.getOutlierProbability(), MAX_REASONABLE_OUTLIER_PROBABILITY, EPSILON_FOR_NEAR_MAX_WARNING);
    }

    private static void warnIfNearMax(final String parameterName,
                                      final double value,
                                      final double maxValue,
                                      final double epsilon) {
        if (maxValue - value < epsilon) {
            logger.warn(String.format("The maximum-likelihood estimate for the global parameter %s (%s) was near its boundary (%s), " +
                            "the model is likely not a good fit to the data!  Consider changing parameters for filtering homozygous sites.",
                    parameterName,
                    CopyNumberFormatsUtils.formatDouble(value),
                    CopyNumberFormatsUtils.formatDouble(maxValue)));
        }
    }

    AlleleFractionState getInitializedState() {
        return new AlleleFractionState(
                globalParameters.getMeanBias(), globalParameters.getBiasVariance(), globalParameters.getOutlierProbability(), minorFractions);
    }

    /**
     * <p>
     *     Initialize minor fractions assuming no allelic bias.
     * </p>
     *
     * <p>
     *     We integrate over f to get posterior probabilities (responsibilities) of alt / ref minor
     *     that is, responsibility of alt minor is int_{0 to 1/2} f^a (1 - f)^r df
     *              responsibility of ref minor is int_{0 to 1/2} f^r (1 - f)^a df
     *     These are proportional to I(1/2, a + 1, r + 1) and I(1/2, r + 1, a + 1),
     *     respectively, where I is the (incomplete) regularized Beta function.
     *     By definition, these likelihoods sum to 1, i.e., they are already normalized.
     * </p>
     *
     * <p>
     *     Finally, we set each minor fraction to the responsibility-weighted total count of
     *     reads in minor allele divided by total reads, ignoring outliers.
     * </p>
     */
    private AlleleFractionState.MinorFractions calculateInitialMinorFractions(final AlleleFractionSegmentedData data) {
        final int numSegments = data.getNumSegments();
        final AlleleFractionState.MinorFractions result = new AlleleFractionState.MinorFractions(numSegments);
        for (int segment = 0; segment < numSegments; segment++) {
            double responsibilityWeightedMinorAlleleReadCount = 0.0;
            double responsibilityWeightedTotalReadCount = 0.0;
            for (final AllelicCount count : data.getIndexedAllelicCountsInSegment(segment)) {
                final int a = count.getAltReadCount();
                final int r = count.getRefReadCount();
                double altMinorResponsibility;
                try {
                    altMinorResponsibility = Beta.regularizedBeta(0.5, a + 1, r + 1);
                } catch (final MaxCountExceededException e) {
                    altMinorResponsibility = a < r ? 1.0 : 0.0; //if the special function can't be computed, give an all-or-nothing responsibility
                }
                responsibilityWeightedMinorAlleleReadCount += altMinorResponsibility * a + (1 - altMinorResponsibility) * r;
                responsibilityWeightedTotalReadCount += a + r;
            }

            // we achieve a flat prior via a single pseudocount for minor and non-minor reads, hence the  +1 and +2
            result.add((responsibilityWeightedMinorAlleleReadCount + 1)/(responsibilityWeightedTotalReadCount + 2));
        }
        return result;
    }

    private double estimateOutlierProbability() {
        final Function<Double, Double> objective = outlierProbability ->
                AlleleFractionLikelihoods.logLikelihood(globalParameters.copyWithNewOutlierProbability(outlierProbability), minorFractions, data);
        return OptimizationUtils.argmax(objective, 0.0, MAX_REASONABLE_OUTLIER_PROBABILITY, globalParameters.getOutlierProbability());
    }

    private double estimateMeanBias() {
        final Function<Double, Double> objective = meanBias ->
                AlleleFractionLikelihoods.logLikelihood(globalParameters.copyWithNewMeanBias(meanBias), minorFractions, data);
        return OptimizationUtils.argmax(objective, 0.0, MAX_REASONABLE_MEAN_BIAS, globalParameters.getMeanBias());
    }

    private double estimateBiasVariance() {
        final Function<Double, Double> objective = biasVariance ->
                AlleleFractionLikelihoods.logLikelihood(globalParameters.copyWithNewBiasVariance(biasVariance), minorFractions, data);
        return OptimizationUtils.argmax(objective, 0.0, MAX_REASONABLE_BIAS_VARIANCE, globalParameters.getBiasVariance());
    }

    private double estimateMinorFraction(final int segment) {
        final Function<Double, Double> objective = minorFraction ->
            AlleleFractionLikelihoods.segmentLogLikelihood(globalParameters, minorFraction, data.getIndexedAllelicCountsInSegment(segment));
        return OptimizationUtils.argmax(objective, 0.0, MAX_MINOR_ALLELE_FRACTION, minorFractions.get(segment));
    }

    private AlleleFractionState.MinorFractions estimateMinorFractions() {
        return new AlleleFractionState.MinorFractions(
                IntStream.range(0, data.getNumSegments()).boxed()
                        .map(this::estimateMinorFraction)
                        .collect(Collectors.toList()));
    }
}
