package org.broadinstitute.hellbender.tools.copynumber.models;

import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CopyRatioCollection;
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
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Represents a segmented model for copy ratio fit to denoised log2 copy-ratio data.
 * The log2 copy ratios in each segment are fit by a mixture model with a normal-distribution component
 * and a uniform outlier component.  The variance of the normal-distribution component and the relative
 * contribution of the uniform outlier component in all segments are both assumed to be global parameters.
 * The mean of the normal-distribution component in each segment is taken to be a segment-level parameter.
 * The component (i.e., normal or outlier) that each copy-ratio point is drawn from is determined by a latent
 * point-level indicator.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class CopyRatioModeller {
    private static final double EPSILON = 1E-6;
    static final double LOG2_COPY_RATIO_MIN = -50.;
    static final double LOG2_COPY_RATIO_MAX = 10.;
    private static final double LOG2_COPY_RATIO_RANGE = LOG2_COPY_RATIO_MAX - LOG2_COPY_RATIO_MIN;
    private static final double VARIANCE_MIN = EPSILON;

    private static final double OUTLIER_PROBABILITY_INITIAL = 0.05;
    private static final double OUTLIER_PROBABILITY_PRIOR_ALPHA = 5.;
    private static final double OUTLIER_PROBABILITY_PRIOR_BETA = 95.;

    private final SampleLocatableMetadata metadata;
    private final ParameterizedModel<CopyRatioParameter, CopyRatioState, CopyRatioSegmentedData> model;

    private final List<Double> varianceSamples = new ArrayList<>();
    private final List<Double> outlierProbabilitySamples = new ArrayList<>();
    private final List<CopyRatioState.SegmentMeans> segmentMeansSamples = new ArrayList<>();

    /**
     * Constructs a copy-ratio model given copy ratios and segments.
     * Initial point estimates of parameters are set to empirical estimates where available.
     */
    CopyRatioModeller(final CopyRatioCollection copyRatios,
                      final SimpleIntervalCollection segments) {
        Utils.nonNull(copyRatios);
        Utils.nonNull(segments);
        Utils.validateArg(copyRatios.getMetadata().getSequenceDictionary().equals(segments.getMetadata().getSequenceDictionary()),
                "Metadata of the copy ratios and the segments do not match.");
        Utils.nonEmpty(segments.getRecords());

        metadata = copyRatios.getMetadata();
        final CopyRatioSegmentedData data = new CopyRatioSegmentedData(copyRatios, segments);

        //set widths for slice sampling of variance and segment-mean posteriors using empirical variance estimate.
        //variance posterior is inverse chi-squared, segment-mean posteriors are Gaussian; the below expressions
        //approximate the standard deviations of these distributions.
        //we also make sure all initial values are within appropriate bounds
        final double dataRangeOrNaN = data.getMaxLog2CopyRatioValue() - data.getMinLog2CopyRatioValue();
        final double dataRange = Double.isNaN(dataRangeOrNaN) ? LOG2_COPY_RATIO_RANGE : dataRangeOrNaN;
        final double varianceEstimateOrNaN = data.estimateVariance();
        final double varianceEstimate = Double.isNaN(varianceEstimateOrNaN) ? VARIANCE_MIN : Math.max(varianceEstimateOrNaN, VARIANCE_MIN);
        final double varianceSliceSamplingWidth = 2. * varianceEstimate;
        final double varianceMax = Math.max(10. * varianceEstimate, dataRange * dataRange);
        final double meanSliceSamplingWidth = Math.sqrt(varianceEstimate * data.getNumSegments() / data.getNumPoints());
        final List<Double> segmentMeans = data.estimateSegmentMeans().stream()
                .map(m -> Math.max(LOG2_COPY_RATIO_MIN, Math.min(LOG2_COPY_RATIO_MAX, m)))
                .collect(Collectors.toList());

        //the uniform log-likelihood for outliers is determined by the minimum and maximum coverages in the dataset;
        //the outlier-probability parameter should be interpreted accordingly
        final double outlierUniformLogLikelihood = -Math.log(dataRange);

        //use empirical segment means and empirical average variance across segments to initialize CopyRatioState
        final CopyRatioState initialState = new CopyRatioState(varianceEstimate, CopyRatioModeller.OUTLIER_PROBABILITY_INITIAL,
                new CopyRatioState.SegmentMeans(segmentMeans), new CopyRatioState.OutlierIndicators(Collections.nCopies(data.getNumPoints(), false)));

        //define ParameterSamplers
        final ParameterSampler<Double, CopyRatioParameter, CopyRatioState, CopyRatioSegmentedData> varianceSampler =
                new CopyRatioSamplers.VarianceSampler(VARIANCE_MIN, varianceMax, varianceSliceSamplingWidth);
        final ParameterSampler<Double, CopyRatioParameter, CopyRatioState, CopyRatioSegmentedData> outlierProbabilitySampler =
                new CopyRatioSamplers.OutlierProbabilitySampler(OUTLIER_PROBABILITY_PRIOR_ALPHA, OUTLIER_PROBABILITY_PRIOR_BETA);
        final ParameterSampler<CopyRatioState.SegmentMeans, CopyRatioParameter, CopyRatioState, CopyRatioSegmentedData> segmentMeansSampler =
                new CopyRatioSamplers.SegmentMeansSampler(LOG2_COPY_RATIO_MIN, LOG2_COPY_RATIO_MAX, meanSliceSamplingWidth);
        final ParameterSampler<CopyRatioState.OutlierIndicators, CopyRatioParameter, CopyRatioState, CopyRatioSegmentedData> outlierIndicatorsSampler =
                new CopyRatioSamplers.OutlierIndicatorsSampler(outlierUniformLogLikelihood);

        model = new ParameterizedModel.GibbsBuilder<>(initialState, data)
                .addParameterSampler(CopyRatioParameter.VARIANCE, varianceSampler, Double.class)
                .addParameterSampler(CopyRatioParameter.OUTLIER_PROBABILITY, outlierProbabilitySampler, Double.class)
                .addParameterSampler(CopyRatioParameter.SEGMENT_MEANS, segmentMeansSampler, CopyRatioState.SegmentMeans.class)
                .addParameterSampler(CopyRatioParameter.OUTLIER_INDICATORS, outlierIndicatorsSampler, CopyRatioState.OutlierIndicators.class)
                .build();
    }

    /**
     * Adds {@code numSamples - numBurnIn} Markov-Chain Monte-Carlo samples of the parameter posteriors (generated using
     * Gibbs sampling) to the collections held internally.  The current {@link CopyRatioState} held internally is used
     * to initialize the Markov Chain.
     * @param numSamples    total number of samples per posterior
     * @param numBurnIn     number of burn-in samples to discard
     */
    void fitMCMC(final int numSamples,
                        final int numBurnIn) {
        ParamUtils.isPositiveOrZero(numBurnIn, "Number of burn-in samples must be non-negative.");
        Utils.validateArg(numBurnIn < numSamples, "Number of samples must be greater than number of burn-in samples.");

        //run MCMC
        final GibbsSampler<CopyRatioParameter, CopyRatioState, CopyRatioSegmentedData> gibbsSampler = new GibbsSampler<>(numSamples, model);
        gibbsSampler.runMCMC();

        //update posterior samples
        varianceSamples.addAll(gibbsSampler.getSamples(CopyRatioParameter.VARIANCE, Double.class, numBurnIn));
        outlierProbabilitySamples.addAll(gibbsSampler.getSamples(CopyRatioParameter.OUTLIER_PROBABILITY, Double.class, numBurnIn));
        segmentMeansSamples.addAll(gibbsSampler.getSamples(CopyRatioParameter.SEGMENT_MEANS, CopyRatioState.SegmentMeans.class, numBurnIn));
    }

    List<Double> getVarianceSamples() {
        return Collections.unmodifiableList(varianceSamples);
    }

    List<Double> getOutlierProbabilitySamples() {
        return Collections.unmodifiableList(outlierProbabilitySamples);
    }

    List<CopyRatioState.SegmentMeans> getSegmentMeansSamples() {
        return Collections.unmodifiableList(segmentMeansSamples);
    }

    /**
     * Should only be called after {@link #fitMCMC} has been called.
     */
    List<ModeledSegment.SimplePosteriorSummary> getSegmentMeansPosteriorSummaries() {
        if (segmentMeansSamples.isEmpty()) {
            throw new IllegalStateException("Attempted to get posterior summaries for segment means before MCMC was performed.");
        }
        final int numSegments = segmentMeansSamples.get(0).size();
        final List<ModeledSegment.SimplePosteriorSummary> posteriorSummaries = new ArrayList<>(numSegments);
        for (int segment = 0; segment < numSegments; segment++) {
            final int j = segment;
            final List<Double> meanSamples =
                    segmentMeansSamples.stream().map(s -> s.get(j)).collect(Collectors.toList());
            posteriorSummaries.add(new ModeledSegment.SimplePosteriorSummary(meanSamples));
        }
        return posteriorSummaries;
    }

    /**
     * Should only be called after {@link #fitMCMC} has been called.
     */
    ParameterDecileCollection<CopyRatioParameter> getGlobalParameterDeciles() {
        if (varianceSamples.isEmpty()) {
            throw new IllegalStateException("Attempted to get posterior summaries for global parameters before MCMC was performed.");
        }
        final Map<CopyRatioParameter, DecileCollection> parameterToDecilesMap = new LinkedHashMap<>();
        parameterToDecilesMap.put(CopyRatioParameter.VARIANCE, new DecileCollection(varianceSamples));
        parameterToDecilesMap.put(CopyRatioParameter.OUTLIER_PROBABILITY, new DecileCollection(outlierProbabilitySamples));
        return new ParameterDecileCollection<>(new SimpleSampleMetadata(metadata.getSampleName()), parameterToDecilesMap, CopyRatioParameter.class);
    }
}
