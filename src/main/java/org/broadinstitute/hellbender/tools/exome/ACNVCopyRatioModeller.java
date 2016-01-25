package org.broadinstitute.hellbender.tools.exome;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.primitives.Doubles;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.mcmc.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Represents an ACNV segmented model for copy ratio fit to target-coverage data.
 * The log_2 coverages in each segment are fit by a mixture model with a normal-distribution component
 * and a uniform outlier component.  The variance of the normal-distribution component and the relative
 * contribution of the uniform outlier component in all segments are both assumed to be global parameters.
 * The mean of the normal-distribution component in each segment is taken to be a segment-level parameter.
 * The component (i.e., normal or outlier) that each target coverage is drawn from is determined by a latent
 * target-level indicator parameter.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class ACNVCopyRatioModeller {
    private static final double EPSILON = 1E-10;
    private static final double VARIANCE_MIN = EPSILON;
    private static final double VARIANCE_MAX = 10.;

    private static final double OUTLIER_PROBABILITY_INITIAL = 0.05;
    private static final double OUTLIER_PROBABILITY_PRIOR_ALPHA = 5.;
    private static final double OUTLIER_PROBABILITY_PRIOR_BETA = 95.;

    private final SegmentedModel segmentedModel;

    private final double outlierUniformLogLikelihood;
    private final double varianceSliceSamplingWidth;
    private final double meanSliceSamplingWidth;

    private final ParameterizedModel<CopyRatioState, CopyRatioDataCollection> model;

    private final List<Double> varianceSamples = new ArrayList<>();
    private final List<Double> outlierProbabilitySamples = new ArrayList<>();
    private final List<SegmentMeans> segmentMeansSamples = new ArrayList<>();
    private final List<OutlierIndicators> outlierIndicatorsSamples = new ArrayList<>();

    //Calculates the exponent for a normal distribution; used in log-likelihood calculation below.
    private static double normalTerm(final double quantity, final double mean, final double variance) {
        return (quantity - mean) * (quantity - mean) / (2. * variance);
    }

    /**
     * Constructs a copy-ratio model given a SegmentedModel with segments and a Genome.  Initial point estimates of
     * parameters are set to empirical estimates where available.
     * @param segmentedModel    SegmentedModel with segments and a Genome
     */
    public ACNVCopyRatioModeller(final SegmentedModel segmentedModel) {
        this.segmentedModel = segmentedModel;

        //load segmented coverages from SegmentedModel into CopyRatioDataCollection
        final CopyRatioDataCollection data = new CopyRatioDataCollection(segmentedModel);

        //the uniform log-likelihood for outliers is determined by the minimum and maximum coverages in the dataset;
        //the outlier-probability parameter should be interpreted accordingly
        outlierUniformLogLikelihood = -Math.log(data.coverageRange);

        //set widths for slice sampling of variance and segment-mean posteriors using empirical variance estimate.
        //variance posterior is inverse chi-squared, segment-mean posteriors are Gaussian; the below expressions
        //approximate the standard deviations of these distributions.
        final double varianceEstimate = data.estimateVariance();
        varianceSliceSamplingWidth = Math.sqrt(2. * varianceEstimate / data.numTargets);
        meanSliceSamplingWidth = Math.sqrt(varianceEstimate * data.numSegments / data.numTargets);

        //use empirical segment means and empirical average variance across segments to initialize CopyRatioState
        final CopyRatioState initialState = new CopyRatioState(data);

        //samples log conditional posterior for the variance parameter, assuming uniform prior; this is given by
        //the product of Gaussian likelihoods for each non-outlier target t:
        //  log[product_{non-outlier t} variance^(-1/2) * exp(-(coverage_t - mean_t)^2 / (2 * variance))] + constant
        //where mean_t is identical for all targets in a segment
        final Sampler<Double, CopyRatioState, CopyRatioDataCollection> varianceSampler =
                (rng, state, dataCollection) -> {
            final Function<Double, Double> logConditionalPDF = newVariance -> {
                double ll = 0.;
                int numNotOutliers = 0;
                for (int segment = 0; segment < dataCollection.numSegments; segment++) {
                    final int numTargetsInSegment = dataCollection.getNumTargetsInSegment(segment);
                    final int targetStartInSegment = dataCollection.getStartTargetInSegment(segment);
                    for (int target = 0; target < numTargetsInSegment; target++) {
                        final double coverage = dataCollection.getCoveragesInSegment(segment).get(target);
                        if (!state.getOutlierIndicator(targetStartInSegment + target)) {
                            ll -= normalTerm(coverage, state.getMeanInSegment(segment), newVariance);
                            numNotOutliers++;
                        }
                    }
                }
                ll -= 0.5 * Math.log(newVariance) * numNotOutliers;
                return ll;
            };

            return new SliceSampler(rng, logConditionalPDF, state.variance(),
                    VARIANCE_MIN, VARIANCE_MAX, varianceSliceSamplingWidth).sample();
        };

        //samples log conditional posterior for the outlier-probability parameter, assuming Beta(alpha, beta) prior;
        //this is given by:
        //  log Beta(alpha + number of outlier targets, beta + number of non-outlier targets) + constant
        final Sampler<Double, CopyRatioState, CopyRatioDataCollection> outlierProbabilitySampler =
                (rng, state, dataCollection) -> {
            final int numOutliers = (int) IntStream.range(0, dataCollection.numTargets).filter(state::getOutlierIndicator).count();
            return new BetaDistribution(
                    rng,
                    OUTLIER_PROBABILITY_PRIOR_ALPHA + numOutliers,
                    OUTLIER_PROBABILITY_PRIOR_BETA + dataCollection.numTargets - numOutliers).sample();
        };

        //samples log conditional posteriors for the segment-mean parameters, assuming uniform priors;
        //for each segment s, this is given by the product of Gaussian likelihoods for each non-outlier target t:
        //  log[product_{non-outlier t in s} exp(-(coverage_t - mean_s)^2 / (2 * variance))] + constant
        final Sampler<SegmentMeans, CopyRatioState, CopyRatioDataCollection> segmentMeansSampler =
                (rng, state, dataCollection) -> {
            final List<Double> means = new ArrayList<>(dataCollection.numSegments);
            for (int i = 0; i < dataCollection.numSegments; i++) {
                //segment, targetStartInSegment, and numTargetsInSegment need to be effectively final in lambda
                final int segment = i;
                final int targetStartInSegment = dataCollection.getStartTargetInSegment(segment);
                final int numTargetsInSegment = dataCollection.getNumTargetsInSegment(segment);
                if (numTargetsInSegment == 0) {
                    means.add(Double.NaN);
                } else {
                    final Function<Double, Double> logConditionalPDF = newMean -> IntStream.range(0, numTargetsInSegment)
                            .filter(target -> !state.getOutlierIndicator(targetStartInSegment + target))
                            .mapToDouble(target -> dataCollection.getCoveragesInSegment(segment).get(target))
                            .map(coverage -> -normalTerm(coverage, newMean, state.variance()))
                            .sum();

                    final SliceSampler sampler =
                            new SliceSampler(rng, logConditionalPDF, state.getMeanInSegment(segment), meanSliceSamplingWidth);
                    means.add(sampler.sample());
                }
            }
            return new SegmentMeans(means);
        };

        //samples log conditional posteriors for the outlier-indicator parameters; for each target t, this is given by:
        //          z_t * [log outlier_prob + outlierUniformLogLikelihood]
        //  + (1 - z_t) * [log(1 - outlier_prob) - log(2 * pi * variance)/2 - (coverage_t - mean_t)^2 / (2 * variance)]
        //  + const
        //where z_t is the indicator for target t, and outlier_prob is the outlier probability.
        //note that we compute the normalizing constant, so that we can sample a new indicator value by simply sampling
        //uniformly in [0, 1] and checking whether the resulting value is less than the probability of being an outlier
        //(corresponding to the first line in the unnormalized expression above)
        final Sampler<OutlierIndicators, CopyRatioState, CopyRatioDataCollection> outlierIndicatorsSampler =
                (rng, state, dataCollection) -> {
            final List<Boolean> indicators = new ArrayList<>();
            final double outlierUnnormalizedLogProbability =
                    Math.log(state.outlierProbability()) + outlierUniformLogLikelihood;
            final double notOutlierUnnormalizedLogProbabilityPrefactor =
                    Math.log(1. - state.outlierProbability()) - 0.5 * Math.log(2 * Math.PI * state.variance());
            for (int segment = 0; segment < dataCollection.numSegments; segment++) {
                final int numTargetsInSegment = dataCollection.getNumTargetsInSegment(segment);
                for (int target = 0; target < numTargetsInSegment; target++) {
                    final double coverage = dataCollection.getCoveragesInSegment(segment).get(target);
                    final double notOutlierUnnormalizedLogProbability =
                            notOutlierUnnormalizedLogProbabilityPrefactor
                                    - normalTerm(coverage, state.getMeanInSegment(segment), state.variance());
                    //note: we are working in natural log space, but this differs from log10 by a multiplicative
                    //constant which is absorbed into the normalization.  Thus normalizeFromLog10 works here.
                    final double conditionalProbability =
                            MathUtils.normalizeFromLog10(new double[]{
                                    outlierUnnormalizedLogProbability,
                                    notOutlierUnnormalizedLogProbability})[0];
                    indicators.add(rng.nextDouble() < conditionalProbability);
                }
            }
            return new OutlierIndicators(indicators);
        };

        model = new ParameterizedModel.GibbsBuilder<>(initialState, data, CopyRatioState.class)
                .addParameterSampler(CopyRatioState.VARIANCE_NAME, varianceSampler, Double.class)
                .addParameterSampler(CopyRatioState.OUTLIER_PROBABILITY_NAME, outlierProbabilitySampler, Double.class)
                .addParameterSampler(CopyRatioState.SEGMENT_MEANS_NAME, segmentMeansSampler, SegmentMeans.class)
                .addParameterSampler(CopyRatioState.OUTLIER_INDICATORS_NAME, outlierIndicatorsSampler,
                        OutlierIndicators.class)
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
        final GibbsSampler<CopyRatioState, CopyRatioDataCollection> gibbsSampler
                = new GibbsSampler<>(numSamples, model);
        gibbsSampler.runMCMC();
        //update posterior samples
        varianceSamples.addAll(gibbsSampler.getSamples(CopyRatioState.VARIANCE_NAME,
                Double.class, numBurnIn));
        outlierProbabilitySamples.addAll(gibbsSampler.getSamples(CopyRatioState.OUTLIER_PROBABILITY_NAME,
                Double.class, numBurnIn));
        segmentMeansSamples.addAll(gibbsSampler.getSamples(CopyRatioState.SEGMENT_MEANS_NAME,
                SegmentMeans.class, numBurnIn));
        outlierIndicatorsSamples.addAll(gibbsSampler.getSamples(CopyRatioState.OUTLIER_INDICATORS_NAME,
                OutlierIndicators.class, numBurnIn));
    }

    /**
     * Returns the {@link SegmentedModel} held internally.
     * @return the {@link SegmentedModel} held internally
     */
    public SegmentedModel getSegmentedModel() {
        return segmentedModel;
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
     * {@link ACNVCopyRatioModeller.SegmentMeans} objects.
     * @return  unmodifiable view of the list of samples of the segment-means posterior
     */
    public List<SegmentMeans> getSegmentMeansSamples() {
        return Collections.unmodifiableList(segmentMeansSamples);
    }

    /**
     * Returns a list of {@link PosteriorSummary} elements summarizing the segment-mean posterior for each segment.
     * Should only be called after {@link ACNVCopyRatioModeller#fitMCMC(int, int)} has been called.
     * @param credibleIntervalAlpha credible-interval alpha, must be in (0, 1)
     * @param ctx                   {@link JavaSparkContext} used for mllib kernel density estimation
     * @return                      list of {@link PosteriorSummary} elements summarizing the
     *                              segment-mean posterior for each segment
     */
    public List<PosteriorSummary> getSegmentMeansPosteriorSummaries(final double credibleIntervalAlpha,
                                                                    final JavaSparkContext ctx) {
        final int numSegments = segmentedModel.getSegments().size();
        final List<PosteriorSummary> posteriorSummaries = new ArrayList<>(numSegments);
        for (int segment = 0; segment < numSegments; segment++) {
            final int j = segment;
            final List<Double> meanSamples =
                    segmentMeansSamples.stream().map(s -> s.getMeanInSegment(j)).collect(Collectors.toList());
            posteriorSummaries.add(PosteriorSummary.calculateHighestPosteriorDensitySummary(meanSamples, credibleIntervalAlpha, ctx));
        }
        return posteriorSummaries;
    }

    /**
     * Returns an unmodifiable view of the list of samples of the outlier-indicators posterior, represented as a list of
     * {@link ACNVCopyRatioModeller.OutlierIndicators} objects.
     * @return  unmodifiable view of the list of samples of the outlier-indicators posterior
     */
    public List<OutlierIndicators> getOutlierIndicatorsSamples() {
        return Collections.unmodifiableList(outlierIndicatorsSamples);
    }

    private final class SegmentMeans {
        private final List<Double> segmentMeans;

        private SegmentMeans(final List<Double> segmentMeans) {
            this.segmentMeans = new ArrayList<>(segmentMeans);
        }

        //returns mean of segment (indexed by segment = 0, ...., numSegments - 1)
        protected double getMeanInSegment(final int segment) {
            return segmentMeans.get(segment);
        }
    }

    @VisibleForTesting
    protected final class OutlierIndicators {
        private final List<Boolean> outlierIndicators;

        private OutlierIndicators(final List<Boolean> outlierIndicators) {
            this.outlierIndicators = new ArrayList<>(outlierIndicators);
        }

        //returns outlier indicator of target (indexed by target = 0, ..., numTargets - 1)
        protected boolean getOutlierIndicator(final int target) {
            return outlierIndicators.get(target);
        }
    }

    private final class CopyRatioState extends AbstractParameterizedState {
        private static final String VARIANCE_NAME = "var";
        private static final String OUTLIER_PROBABILITY_NAME = "pi";
        private static final String SEGMENT_MEANS_NAME = "mu";
        private static final String OUTLIER_INDICATORS_NAME = "z";

        @Override
        protected <S extends AbstractParameterizedState> S copy(final Class<S> stateClass) {
            return stateClass.cast(new CopyRatioState(this));
        }

        private CopyRatioState(final CopyRatioState state) {
            super(state);
        }

        private CopyRatioState(final double variance, final double outlierProbability,
                               final SegmentMeans segmentMeans, final OutlierIndicators outlierIndicators) {
            super(Arrays.asList(
                    new Parameter<>(VARIANCE_NAME, variance),
                    new Parameter<>(OUTLIER_PROBABILITY_NAME, outlierProbability),
                    new Parameter<>(SEGMENT_MEANS_NAME, segmentMeans),
                    new Parameter<>(OUTLIER_INDICATORS_NAME, outlierIndicators)));
        }

        private CopyRatioState(final CopyRatioDataCollection dataCollection) {
            this(dataCollection.estimateVariance(), OUTLIER_PROBABILITY_INITIAL, dataCollection.estimateSegmentMeans(),
                    new OutlierIndicators(Collections.nCopies(dataCollection.numTargets, false)));
        }

        private double variance() {
            return get(VARIANCE_NAME, Double.class);
        }

        private double outlierProbability()  { return get(OUTLIER_PROBABILITY_NAME, Double.class); }

        private double getMeanInSegment(final int segment) {
            return get(SEGMENT_MEANS_NAME, SegmentMeans.class).getMeanInSegment(segment);
        }

        private boolean getOutlierIndicator(final int target) {
            return get(OUTLIER_INDICATORS_NAME, OutlierIndicators.class).getOutlierIndicator(target);
        }
    }

    private final class CopyRatioDataCollection implements DataCollection {
        private final int numSegments;
        private final int numTargets;
        private final double coverageRange;

        private final List<List<Double>> coveragesPerSegment = new ArrayList<>();
        private final List<Integer> numTargetsPerSegment;
        private final List<Integer> startTargetsPerSegment = new ArrayList<>();

        private CopyRatioDataCollection(final SegmentedModel segmentedModel) {
            if (segmentedModel.getGenome().getTargets().targetCount() == 0) {
                throw new IllegalArgumentException("Cannot construct CopyRatioDataCollection with no target-coverage data.");
            }
            //construct coverages and number of targets per segment from the segments and Genome in SegmentedModel
            final TargetCollection<TargetCoverage> targetCoverages = segmentedModel.getGenome().getTargets();
            //construct list of coverages (in order corresponding to that of segments in SegmentedModel;
            //this may not be in genomic order, depending on how the segments are sorted in the segment file,
            //so we cannot simply take the list of coverages in the order from TargetCollection.targets()
            final List<Double> coverages = segmentedModel.getSegments().stream()
                    .flatMap(s -> targetCoverages.targets(s).stream())
                    .map(TargetCoverage::getCoverage)
                    .collect(Collectors.toList());
            numTargetsPerSegment = segmentedModel.getSegments().stream()
                    .map(s -> targetCoverages.targets(s).size())
                    .collect(Collectors.toList());
            numSegments = numTargetsPerSegment.size();
            numTargets = coverages.size();
            coverageRange = Collections.max(coverages) - Collections.min(coverages);
            //partition coverages by segment
            int startTarget = 0;
            for (int segment = 0; segment < numSegments; segment++) {
                startTargetsPerSegment.add(startTarget);
                final int numTargetsInSegment = numTargetsPerSegment.get(segment);
                coveragesPerSegment.add(coverages.subList(startTarget, startTarget + numTargetsInSegment));
                startTarget += numTargetsInSegment;
            }
        }

        private List<Double> getCoveragesInSegment(final int segment) {
            return coveragesPerSegment.get(segment);
        }

        private int getNumTargetsInSegment(final int segment) {
            return numTargetsPerSegment.get(segment);
        }

        private int getStartTargetInSegment(final int segment) {
            return startTargetsPerSegment.get(segment);
        }

        //estimate global variance empirically by taking average of all per-segment variances
        private double estimateVariance() {
            final double[] variancesPerSegment = new double[numSegments];
            for (int segment = 0; segment < numSegments; segment++) {
                variancesPerSegment[segment] = new Variance().evaluate(Doubles.toArray(getCoveragesInSegment(segment)));
            }
            return new Mean().evaluate(variancesPerSegment);
        }

        //estimate segment means empirically by taking averages of coverages in each segment
        private SegmentMeans estimateSegmentMeans() {
            final List<Double> means = new ArrayList<>(numSegments);
            for (int segment = 0; segment < numSegments; segment++) {
                means.add(new Mean().evaluate(Doubles.toArray(getCoveragesInSegment(segment))));
            }
            return new SegmentMeans(means);
        }
    }
}
