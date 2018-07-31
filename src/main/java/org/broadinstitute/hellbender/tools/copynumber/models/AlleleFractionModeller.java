package org.broadinstitute.hellbender.tools.copynumber.models;

import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.ParameterDecileCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SimpleIntervalCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleSampleMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.ModeledSegment;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.mcmc.DecileCollection;
import org.broadinstitute.hellbender.utils.mcmc.GibbsSampler;
import org.broadinstitute.hellbender.utils.mcmc.ParameterSampler;
import org.broadinstitute.hellbender.utils.mcmc.ParameterizedModel;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Given segments and counts of alt and ref reads over a list of het sites,
 * infers the minor-allele fraction of each segment.  For example, a segment
 * with (alt,ref) counts (10,90), (11,93), (88,12), (90,10) probably has a minor-allele fraction
 * somewhere around 0.1.  The model takes into account allelic reference bias due to mapping etc. by learning
 * a global gamma distribution on allelic bias ratios.
 *
 * <p>
 *     We define the bias ratio of each het locus to be the expected ratio of
 *     mapped ref reads to mapped alt reads given equal amounts of DNA (that is, given
 *     a germline het).  The model learns a common gamma distribution:
 *     bias ratio ~ Gamma(alpha = mu^2 / sigma^2, beta = mu / sigma^2)
 *     where mu and sigma^2 are the global mean and variance of bias ratios, and
 *     alpha, beta are the natural parameters of the gamma distribution.
 * </p>
 * <p>
 *     Each segment has a minor-allele fraction f, and for each het within the locus
 *     the number of alt reads is drawn from a binomial distribution with total count
 *     n = #alt reads + #ref reads and alt probability f / (f + (1 - f) * bias ratio) if the
 *     locus is alt minor and (1 - f) / (1 - f + f * bias ratio) if the locus is ref minor.
 *     We also allow a prior on minor-allele fraction to be specified by the alpha parameter of the 4-parameter
 *     beta distribution Beta(alpha, 1, 0, 1/2).
 * </p>
 * <p>
 *     Conceptually, the model contains latent variables corresponding to the bias ratio
 *     and indicators for alt minor/ref minor at each het locus.  However, we integrate them
 *     out and the MCMC model below only contains the minor-allele fractions and
 *     the three hyperparameters of the model: the two parameters of the gamma distribution
 *     along with the global outlier probability.
 * </p>
 *
 * See docs/CNVs/CNV-methods.pdf for a thorough description of the model.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AlleleFractionModeller {
    private static final double MAX_REASONABLE_MEAN_BIAS = AlleleFractionInitializer.MAX_REASONABLE_MEAN_BIAS;
    private static final double MAX_REASONABLE_BIAS_VARIANCE = AlleleFractionInitializer.MAX_REASONABLE_BIAS_VARIANCE;
    private static final double MAX_REASONABLE_OUTLIER_PROBABILITY = AlleleFractionInitializer.MAX_REASONABLE_OUTLIER_PROBABILITY;
    private static final double MIN_MINOR_FRACTION_SAMPLING_WIDTH = 1E-3;

    private final SampleLocatableMetadata metadata;
    private final ParameterizedModel<AlleleFractionParameter, AlleleFractionState, AlleleFractionSegmentedData> model;

    private final List<Double> meanBiasSamples = new ArrayList<>();
    private final List<Double> biasVarianceSamples = new ArrayList<>();
    private final List<Double> outlierProbabilitySamples = new ArrayList<>();
    private final List<AlleleFractionState.MinorFractions> minorFractionsSamples = new ArrayList<>();

    /**
     * Constructs an allele-fraction model given allelic counts and segments.
     * {@link AlleleFractionInitializer} is used for initialization and slice-sampling widths are estimated.
     */
    AlleleFractionModeller(final AllelicCountCollection allelicCounts,
                           final SimpleIntervalCollection segments,
                           final AlleleFractionPrior prior) {
        Utils.nonNull(allelicCounts);
        Utils.nonNull(segments);
        Utils.validateArg(allelicCounts.getMetadata().getSequenceDictionary().equals(segments.getMetadata().getSequenceDictionary()),
                "Metadata of the allelic counts and the segments do not match.");
        Utils.nonNull(prior);

        metadata = allelicCounts.getMetadata();
        final AlleleFractionSegmentedData data = new AlleleFractionSegmentedData(allelicCounts, segments);

        //initialization gets us to the mode of the likelihood
        final AlleleFractionState initialState = new AlleleFractionInitializer(data).getInitializedState();
        final AlleleFractionGlobalParameters initialParameters = initialState.globalParameters();
        final AlleleFractionState.MinorFractions initialMinorFractions = initialState.minorFractions();

        //if we approximate conditionals as normal, we can guess the width from the curvature at the mode and use as the slice-sampling widths
        final double meanBiasSamplingWidths = approximatePosteriorWidthAtMode(meanBias ->
                AlleleFractionLikelihoods.logLikelihood(initialParameters.copyWithNewMeanBias(meanBias), initialMinorFractions, data), initialParameters.getMeanBias());
        final double biasVarianceSamplingWidths = approximatePosteriorWidthAtMode(biasVariance ->
                AlleleFractionLikelihoods.logLikelihood(initialParameters.copyWithNewBiasVariance(biasVariance), initialMinorFractions, data), initialParameters.getBiasVariance());
        final double outlierProbabilitySamplingWidths = approximatePosteriorWidthAtMode(outlierProbability ->
                AlleleFractionLikelihoods.logLikelihood(initialParameters.copyWithNewOutlierProbability(outlierProbability), initialMinorFractions, data), initialParameters.getOutlierProbability());

        final List<Double> minorFractionsSliceSamplingWidths = IntStream.range(0, data.getNumSegments()).boxed()
                .map(segment -> approximatePosteriorWidthAtMode(
                        f -> AlleleFractionLikelihoods.segmentLogLikelihood(initialParameters, f, data.getIndexedAllelicCountsInSegment(segment)), initialMinorFractions.get(segment)))
                .map(w -> Math.max(w, MIN_MINOR_FRACTION_SAMPLING_WIDTH))
                .collect(Collectors.toList());

        final ParameterSampler<Double, AlleleFractionParameter, AlleleFractionState, AlleleFractionSegmentedData> meanBiasSampler =
                new AlleleFractionSamplers.MeanBiasSampler(MAX_REASONABLE_MEAN_BIAS, meanBiasSamplingWidths);
        final ParameterSampler<Double, AlleleFractionParameter, AlleleFractionState, AlleleFractionSegmentedData> biasVarianceSampler =
                new AlleleFractionSamplers.BiasVarianceSampler(MAX_REASONABLE_BIAS_VARIANCE, biasVarianceSamplingWidths);
        final ParameterSampler<Double, AlleleFractionParameter, AlleleFractionState, AlleleFractionSegmentedData> outlierProbabilitySampler =
                new AlleleFractionSamplers.OutlierProbabilitySampler(MAX_REASONABLE_OUTLIER_PROBABILITY, outlierProbabilitySamplingWidths);
        final ParameterSampler<AlleleFractionState.MinorFractions, AlleleFractionParameter, AlleleFractionState, AlleleFractionSegmentedData> minorFractionsSampler =
                new AlleleFractionSamplers.MinorFractionsSampler(prior, minorFractionsSliceSamplingWidths);

        model = new ParameterizedModel.GibbsBuilder<>(initialState, data)
                .addParameterSampler(AlleleFractionParameter.MEAN_BIAS, meanBiasSampler, Double.class)
                .addParameterSampler(AlleleFractionParameter.BIAS_VARIANCE, biasVarianceSampler, Double.class)
                .addParameterSampler(AlleleFractionParameter.OUTLIER_PROBABILITY, outlierProbabilitySampler, Double.class)
                .addParameterSampler(AlleleFractionParameter.MINOR_ALLELE_FRACTIONS, minorFractionsSampler, AlleleFractionState.MinorFractions.class)
                .build();
    }

    /**
     * Adds {@code numSamples - numBurnIn} Markov-Chain Monte-Carlo samples of the parameter posteriors (generated using
     * Gibbs sampling) to the collections held internally.  The current {@link AlleleFractionState} held internally is used
     * to initialize the Markov Chain.
     * @param numSamples    total number of samples per posterior
     * @param numBurnIn     number of burn-in samples to discard
     */
    void fitMCMC(final int numSamples, final int numBurnIn) {
        //run MCMC
        final GibbsSampler<AlleleFractionParameter, AlleleFractionState, AlleleFractionSegmentedData> gibbsSampler = new GibbsSampler<>(numSamples, model);
        gibbsSampler.runMCMC();

        //update posterior samples
        meanBiasSamples.addAll(gibbsSampler.getSamples(AlleleFractionParameter.MEAN_BIAS, Double.class, numBurnIn));
        biasVarianceSamples.addAll(gibbsSampler.getSamples(AlleleFractionParameter.BIAS_VARIANCE, Double.class, numBurnIn));
        outlierProbabilitySamples.addAll(gibbsSampler.getSamples(AlleleFractionParameter.OUTLIER_PROBABILITY, Double.class, numBurnIn));
        minorFractionsSamples.addAll(gibbsSampler.getSamples(AlleleFractionParameter.MINOR_ALLELE_FRACTIONS, AlleleFractionState.MinorFractions.class, numBurnIn));
    }

    List<Double> getMeanBiasSamples() {
        return Collections.unmodifiableList(meanBiasSamples);
    }

    List<Double> getBiasVarianceSamples() {
        return Collections.unmodifiableList(biasVarianceSamples);
    }

    List<Double> getOutlierProbabilitySamples() {
        return Collections.unmodifiableList(outlierProbabilitySamples);
    }

    List<AlleleFractionState.MinorFractions> getMinorFractionsSamples() {
        return Collections.unmodifiableList(minorFractionsSamples);
    }

    /**
     * Should only be called after {@link #fitMCMC} has been called.
     */
    List<ModeledSegment.SimplePosteriorSummary> getMinorAlleleFractionsPosteriorSummaries() {
        if (minorFractionsSamples.isEmpty()) {
            throw new IllegalStateException("Attempted to get posterior summaries for minor-allele fractions before MCMC was performed.");
        }
        final int numSegments = minorFractionsSamples.get(0).size();
        final List<ModeledSegment.SimplePosteriorSummary> posteriorSummaries = new ArrayList<>(numSegments);
        for (int segment = 0; segment < numSegments; segment++) {
            final int j = segment;
            final List<Double> minorFractionSamples =
                    minorFractionsSamples.stream().map(s -> s.get(j)).collect(Collectors.toList());
            posteriorSummaries.add(new ModeledSegment.SimplePosteriorSummary(minorFractionSamples));
        }
        return posteriorSummaries;
    }

    /**
     * Should only be called after {@link #fitMCMC} has been called.
     */
    ParameterDecileCollection<AlleleFractionParameter> getGlobalParameterDeciles() {
        if (meanBiasSamples.isEmpty()) {
            throw new IllegalStateException("Attempted to get posterior summaries for global parameters before MCMC was performed.");
        }
        final Map<AlleleFractionParameter, DecileCollection> parameterToDecilesMap = new LinkedHashMap<>();
        parameterToDecilesMap.put(AlleleFractionParameter.MEAN_BIAS, new DecileCollection(meanBiasSamples));
        parameterToDecilesMap.put(AlleleFractionParameter.BIAS_VARIANCE, new DecileCollection(biasVarianceSamples));
        parameterToDecilesMap.put(AlleleFractionParameter.OUTLIER_PROBABILITY, new DecileCollection(outlierProbabilitySamples));
        return new ParameterDecileCollection<>(new SimpleSampleMetadata(metadata.getSampleName()), parameterToDecilesMap, AlleleFractionParameter.class);
    }

    //use width of a probability distribution given the position of its mode (estimated from Gaussian approximation) as step size
    private static double approximatePosteriorWidthAtMode(final Function<Double, Double> logPDF,
                                                          final double mode) {
        final double absMode = Math.abs(mode);
        final double epsilon = Math.min(1E-6, absMode / 2);    //adjust scale if mode is very near zero
        final double defaultWidth = absMode / 10;              //if "mode" is not close to true mode of logPDF, approximation may not apply; just use 1 / 10 of absMode in this case
        final double secondDerivative = (logPDF.apply(mode + epsilon) - 2 * logPDF.apply(mode) + logPDF.apply(mode - epsilon)) / (epsilon * epsilon);
        return secondDerivative < 0 ? Math.sqrt(-1.0 / secondDerivative) : defaultWidth;
    }
}
