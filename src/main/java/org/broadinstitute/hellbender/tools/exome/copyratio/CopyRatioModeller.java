package org.broadinstitute.hellbender.tools.exome.copyratio;

import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.tools.exome.Genome;
import org.broadinstitute.hellbender.tools.exome.SegmentedGenome;
import org.broadinstitute.hellbender.utils.mcmc.*;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Represents an ACNV segmented model for copy ratio fit to tangent-normalized log_2 target-coverage data.
 * The log_2 coverages in each segment are fit by a mixture model with a normal-distribution component
 * and a uniform outlier component.  The variance of the normal-distribution component and the relative
 * contribution of the uniform outlier component in all segments are both assumed to be global parameters.
 * The mean of the normal-distribution component in each segment is taken to be a segment-level parameter.
 * The component (i.e., normal or outlier) that each target coverage is drawn from is determined by a latent
 * target-level indicator parameter.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class CopyRatioModeller {
    private static final double EPSILON = 1E-10;
    private static final double VARIANCE_MIN = EPSILON;
    private static final double VARIANCE_MAX = 10.;

    private static final double OUTLIER_PROBABILITY_INITIAL = 0.05;
    private static final double OUTLIER_PROBABILITY_PRIOR_ALPHA = 5.;
    private static final double OUTLIER_PROBABILITY_PRIOR_BETA = 95.;

    private final SegmentedGenome segmentedGenome;
    private final ParameterizedModel<CopyRatioParameter, CopyRatioState, CopyRatioData> model;

    private final List<Double> varianceSamples = new ArrayList<>();
    private final List<Double> outlierProbabilitySamples = new ArrayList<>();
    private final List<CopyRatioState.SegmentMeans> segmentMeansSamples = new ArrayList<>();
    private final List<CopyRatioState.OutlierIndicators> outlierIndicatorsSamples = new ArrayList<>();

    /**
     * Constructs a copy-ratio model given a {@link SegmentedGenome} with segments and a {@link Genome}.
     * Initial point estimates of parameters are set to empirical estimates where available.
     * @param segmentedGenome    SegmentedGenome with segments and a Genome
     */
    public CopyRatioModeller(final SegmentedGenome segmentedGenome) {
        this.segmentedGenome = segmentedGenome;

        //load segmented coverages from SegmentedGenome into CopyRatioData
        final CopyRatioData data = new CopyRatioData(segmentedGenome);

        //set widths for slice sampling of variance and segment-mean posteriors using empirical variance estimate.
        //variance posterior is inverse chi-squared, segment-mean posteriors are Gaussian; the below expressions
        //approximate the standard deviations of these distributions.
        final double varianceEstimate = data.estimateVariance();
        final double varianceSliceSamplingWidth = Math.sqrt(2. * varianceEstimate / data.getNumTargets());
        final double coverageMin = data.getCoverageMin();
        final double coverageMax = data.getCoverageMax();
        final double meanSliceSamplingWidth = Math.sqrt(varianceEstimate * data.getNumSegments() / data.getNumTargets());

        //the uniform log-likelihood for outliers is determined by the minimum and maximum coverages in the dataset;
        //the outlier-probability parameter should be interpreted accordingly
        final double outlierUniformLogLikelihood = -Math.log(coverageMax - coverageMin);

        //use empirical segment means and empirical average variance across segments to initialize CopyRatioState
        final CopyRatioState initialState = new CopyRatioState(varianceEstimate, CopyRatioModeller.OUTLIER_PROBABILITY_INITIAL,
                        data.estimateSegmentMeans(), new CopyRatioState.OutlierIndicators(Collections.nCopies(data.getNumTargets(), false)));

        //define ParameterSamplers
        final ParameterSampler<Double, CopyRatioParameter, CopyRatioState, CopyRatioData> varianceSampler =
                new CopyRatioSamplers.VarianceSampler(VARIANCE_MIN, VARIANCE_MAX, varianceSliceSamplingWidth);
        final ParameterSampler<Double, CopyRatioParameter, CopyRatioState, CopyRatioData> outlierProbabilitySampler =
                new CopyRatioSamplers.OutlierProbabilitySampler(OUTLIER_PROBABILITY_PRIOR_ALPHA, OUTLIER_PROBABILITY_PRIOR_BETA);
        final ParameterSampler<CopyRatioState.SegmentMeans, CopyRatioParameter, CopyRatioState, CopyRatioData> segmentMeansSampler =
                new CopyRatioSamplers.SegmentMeansSampler(coverageMin, coverageMax, meanSliceSamplingWidth);
        final ParameterSampler<CopyRatioState.OutlierIndicators, CopyRatioParameter, CopyRatioState, CopyRatioData> outlierIndicatorsSampler =
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
    public void fitMCMC(final int numSamples, final int numBurnIn) {
        //run MCMC
        final GibbsSampler<CopyRatioParameter, CopyRatioState, CopyRatioData> gibbsSampler
                = new GibbsSampler<>(numSamples, model);
        gibbsSampler.runMCMC();
        //update posterior samples
        varianceSamples.addAll(gibbsSampler.getSamples(CopyRatioParameter.VARIANCE,
                Double.class, numBurnIn));
        outlierProbabilitySamples.addAll(gibbsSampler.getSamples(CopyRatioParameter.OUTLIER_PROBABILITY,
                Double.class, numBurnIn));
        segmentMeansSamples.addAll(gibbsSampler.getSamples(CopyRatioParameter.SEGMENT_MEANS,
                CopyRatioState.SegmentMeans.class, numBurnIn));
        outlierIndicatorsSamples.addAll(gibbsSampler.getSamples(CopyRatioParameter.OUTLIER_INDICATORS,
                CopyRatioState.OutlierIndicators.class, numBurnIn));
    }

    /**
     * Returns the {@link SegmentedGenome} held internally.
     * @return the {@link SegmentedGenome} held internally
     */
    public SegmentedGenome getSegmentedGenome() {
        return segmentedGenome;
    }

    /**
     * Returns an unmodifiable view of the list of samples of the variance posterior.
     * @return  unmodifiable view of the list of samples of the variance posterior
     */
    public List<Double> getVarianceSamples() {
        return Collections.unmodifiableList(varianceSamples);
    }

    /**
     * Returns an unmodifiable view of the list of samples of the outlier-probability posterior.
     * @return  unmodifiable view of the list of samples of the outlier-probability posterior
     */
    public List<Double> getOutlierProbabilitySamples() {
        return Collections.unmodifiableList(outlierProbabilitySamples);
    }

    /**
     * Returns an unmodifiable view of the list of samples of the segment-means posterior, represented as a list of
     * {@link CopyRatioState.SegmentMeans} objects.
     * @return  unmodifiable view of the list of samples of the segment-means posterior
     */
    public List<CopyRatioState.SegmentMeans> getSegmentMeansSamples() {
        return Collections.unmodifiableList(segmentMeansSamples);
    }

    /**
     * Returns an unmodifiable view of the list of samples of the outlier-indicators posterior, represented as a list of
     * {@link CopyRatioState.OutlierIndicators} objects.
     * @return  unmodifiable view of the list of samples of the outlier-indicators posterior
     */
    public List<CopyRatioState.OutlierIndicators> getOutlierIndicatorsSamples() {
        return Collections.unmodifiableList(outlierIndicatorsSamples);
    }

    /**
     * Returns a list of {@link PosteriorSummary} elements summarizing the segment-mean posterior for each segment.
     * Should only be called after {@link CopyRatioModeller#fitMCMC(int, int)} has been called.
     * @param credibleIntervalAlpha credible-interval alpha, must be in (0, 1)
     * @param ctx                   {@link JavaSparkContext} used for mllib kernel density estimation
     * @return                      list of {@link PosteriorSummary} elements summarizing the
     *                              segment-mean posterior for each segment
     */
    public List<PosteriorSummary> getSegmentMeansPosteriorSummaries(final double credibleIntervalAlpha,
                                                                    final JavaSparkContext ctx) {
        final int numSegments = segmentedGenome.getSegments().size();
        final List<PosteriorSummary> posteriorSummaries = new ArrayList<>(numSegments);
        for (int segment = 0; segment < numSegments; segment++) {
            final int j = segment;
            final List<Double> meanSamples =
                    segmentMeansSamples.stream().map(s -> s.get(j)).collect(Collectors.toList());
            posteriorSummaries.add(PosteriorSummaryUtils.calculateHighestPosteriorDensityAndDecilesSummary(meanSamples, credibleIntervalAlpha, ctx));
        }
        return posteriorSummaries;
    }

    /**
     * Returns a Map of {@link PosteriorSummary} elements summarizing the global parameters.
     * Should only be called after {@link CopyRatioModeller#fitMCMC(int, int)} has been called.
     * @param credibleIntervalAlpha credible-interval alpha, must be in (0, 1)
     * @param ctx                   {@link JavaSparkContext} used for mllib kernel density estimation
     * @return                      list of {@link PosteriorSummary} elements summarizing the global parameters
     */
    public Map<CopyRatioParameter, PosteriorSummary> getGlobalParameterPosteriorSummaries(final double credibleIntervalAlpha, final JavaSparkContext ctx) {
        final Map<CopyRatioParameter, PosteriorSummary> posteriorSummaries = new LinkedHashMap<>();
        posteriorSummaries.put(CopyRatioParameter.VARIANCE, PosteriorSummaryUtils.calculateHighestPosteriorDensityAndDecilesSummary(varianceSamples, credibleIntervalAlpha, ctx));
        posteriorSummaries.put(CopyRatioParameter.OUTLIER_PROBABILITY, PosteriorSummaryUtils.calculateHighestPosteriorDensityAndDecilesSummary(outlierProbabilitySamples, credibleIntervalAlpha, ctx));
        return posteriorSummaries;
    }
}
