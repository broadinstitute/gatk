package org.broadinstitute.hellbender.tools.copynumber.models;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.OverlapDetector;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.*;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.AllelicCount;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyRatio;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.ModeledSegment;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Represents a segmented model for copy ratio and allele fraction.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class MultidimensionalModeller {
    private static final Logger logger = LogManager.getLogger(MultidimensionalModeller.class);

    private final SampleLocatableMetadata metadata;
    private final CopyRatioCollection denoisedCopyRatios;
    private final OverlapDetector<CopyRatio> copyRatioMidpointOverlapDetector;
    private final AllelicCountCollection allelicCounts;
    private final OverlapDetector<AllelicCount> allelicCountOverlapDetector;
    private final AlleleFractionPrior alleleFractionPrior;

    private CopyRatioModeller copyRatioModeller;
    private AlleleFractionModeller alleleFractionModeller;

    private SimpleIntervalCollection currentSegments;
    private final List<ModeledSegment> modeledSegments = new ArrayList<>();

    //similar-segment merging may leave model in a state where it is not properly fit (deciles may be estimated naively)
    private boolean isModelFit;

    private final int numSamplesCopyRatio;
    private final int numBurnInCopyRatio;
    private final int numSamplesAlleleFraction;
    private final int numBurnInAlleleFraction;

    /**
     * Constructs a copy-ratio and allele-fraction modeller, specifying number of total samples
     * and number of burn-in samples for Markov-Chain Monte Carlo model fitting.
     * An initial model fit is performed.
     */
    public MultidimensionalModeller(final MultidimensionalSegmentCollection multidimensionalSegments,
                                    final CopyRatioCollection denoisedCopyRatios,
                                    final AllelicCountCollection allelicCounts,
                                    final AlleleFractionPrior alleleFractionPrior,
                                    final int numSamplesCopyRatio,
                                    final int numBurnInCopyRatio,
                                    final int numSamplesAlleleFraction,
                                    final int numBurnInAlleleFraction) {
        Utils.validateArg(Stream.of(
                Utils.nonNull(multidimensionalSegments).getMetadata(),
                Utils.nonNull(denoisedCopyRatios).getMetadata(),
                Utils.nonNull(allelicCounts).getMetadata()).distinct().count() == 1,
                "Metadata from all inputs must match.");
        ParamUtils.isPositive(multidimensionalSegments.size(), "Number of segments must be positive.");
        metadata = multidimensionalSegments.getMetadata();
        currentSegments = new SimpleIntervalCollection(
                new SimpleLocatableMetadata(metadata.getSequenceDictionary()),
                multidimensionalSegments.getIntervals());
        this.denoisedCopyRatios = denoisedCopyRatios;
        copyRatioMidpointOverlapDetector = denoisedCopyRatios.getMidpointOverlapDetector();
        this.allelicCounts = allelicCounts;
        allelicCountOverlapDetector = allelicCounts.getOverlapDetector();
        this.alleleFractionPrior = Utils.nonNull(alleleFractionPrior);
        this.numSamplesCopyRatio = numSamplesCopyRatio;
        this.numBurnInCopyRatio = numBurnInCopyRatio;
        this.numSamplesAlleleFraction = numSamplesAlleleFraction;
        this.numBurnInAlleleFraction = numBurnInAlleleFraction;
        logger.info("Fitting initial model...");
        fitModel();
    }

    public ModeledSegmentCollection getModeledSegments() {
        return new ModeledSegmentCollection(metadata, modeledSegments);
    }

    /**
     * Performs Markov-Chain Monte Carlo model fitting using the
     * number of total samples and number of burn-in samples specified at construction.
     */
    private void fitModel() {
        //perform MCMC to generate posterior samples
        logger.info("Fitting copy-ratio model...");
        copyRatioModeller = new CopyRatioModeller(denoisedCopyRatios, currentSegments);
        copyRatioModeller.fitMCMC(numSamplesCopyRatio, numBurnInCopyRatio);
        logger.info("Fitting allele-fraction model...");
        alleleFractionModeller = new AlleleFractionModeller(allelicCounts, currentSegments, alleleFractionPrior);
        alleleFractionModeller.fitMCMC(numSamplesAlleleFraction, numBurnInAlleleFraction);

        //update list of ModeledSegment with new PosteriorSummaries
        modeledSegments.clear();
        final List<ModeledSegment.SimplePosteriorSummary> segmentMeansPosteriorSummaries =
                copyRatioModeller.getSegmentMeansPosteriorSummaries();
        final List<ModeledSegment.SimplePosteriorSummary> minorAlleleFractionsPosteriorSummaries =
                alleleFractionModeller.getMinorAlleleFractionsPosteriorSummaries();
        for (int segmentIndex = 0; segmentIndex < currentSegments.size(); segmentIndex++) {
            final SimpleInterval segment = currentSegments.getRecords().get(segmentIndex);
            final int numPointsCopyRatio = copyRatioMidpointOverlapDetector.getOverlaps(segment).size();
            final int numPointsAlleleFraction = allelicCountOverlapDetector.getOverlaps(segment).size();
            final ModeledSegment.SimplePosteriorSummary segmentMeansPosteriorSummary = segmentMeansPosteriorSummaries.get(segmentIndex);
            final ModeledSegment.SimplePosteriorSummary minorAlleleFractionPosteriorSummary = minorAlleleFractionsPosteriorSummaries.get(segmentIndex);
            modeledSegments.add(new ModeledSegment(
                    segment, numPointsCopyRatio, numPointsAlleleFraction, segmentMeansPosteriorSummary, minorAlleleFractionPosteriorSummary));
        }
        isModelFit = true;
    }

    /**
     * @param numSmoothingIterationsPerFit  if this is zero, no refitting will be performed between smoothing iterations
     */
    public void smoothSegments(final int maxNumSmoothingIterations,
                               final int numSmoothingIterationsPerFit,
                               final double smoothingCredibleIntervalThresholdCopyRatio,
                               final double smoothingCredibleIntervalThresholdAlleleFraction) {
        ParamUtils.isPositiveOrZero(maxNumSmoothingIterations,
                "The maximum number of smoothing iterations must be non-negative.");
        ParamUtils.isPositiveOrZero(smoothingCredibleIntervalThresholdCopyRatio,
                "The number of smoothing iterations per fit must be non-negative.");
        ParamUtils.isPositiveOrZero(smoothingCredibleIntervalThresholdAlleleFraction,
                "The allele-fraction credible-interval threshold for segmentation smoothing must be non-negative.");
        logger.info(String.format("Initial number of segments before smoothing: %d", modeledSegments.size()));
        //perform iterations of similar-segment merging until all similar segments are merged
        for (int numIterations = 1; numIterations <= maxNumSmoothingIterations; numIterations++) {
            logger.info(String.format("Smoothing iteration: %d", numIterations));
            final int prevNumSegments = modeledSegments.size();
            if (numSmoothingIterationsPerFit > 0 && numIterations % numSmoothingIterationsPerFit == 0) {
                //refit model after this merge iteration
                performSmoothingIteration(smoothingCredibleIntervalThresholdCopyRatio, smoothingCredibleIntervalThresholdAlleleFraction, true);
            } else {
                //do not refit model after this merge iteration (posterior modes will be identical to posterior medians)
                performSmoothingIteration(smoothingCredibleIntervalThresholdCopyRatio, smoothingCredibleIntervalThresholdAlleleFraction, false);
            }
            if (modeledSegments.size() == prevNumSegments) {
                break;
            }
        }
        if (!isModelFit) {
            //make sure final model is completely fit (i.e., posterior modes are specified)
            fitModel();
        }
        logger.info(String.format("Final number of segments after smoothing: %d", modeledSegments.size()));
    }

    /**
     * Performs one iteration of similar-segment merging on the list of {@link ModeledSegment} held internally.
     * Markov-Chain Monte Carlo model fitting is optionally performed after each iteration using the
     * number of total samples and number of burn-in samples specified at construction.
     * @param intervalThresholdSegmentMean         threshold number of credible intervals for segment-mean similarity
     * @param intervalThresholdMinorAlleleFraction threshold number of credible intervals for minor-allele-fraction similarity
     * @param doModelFit                           if true, refit MCMC model after merging
     */
    private void performSmoothingIteration(final double intervalThresholdSegmentMean,
                                           final double intervalThresholdMinorAlleleFraction,
                                           final boolean doModelFit) {
        logger.info("Number of segments before smoothing iteration: " + modeledSegments.size());
        final List<ModeledSegment> mergedSegments = SimilarSegmentUtils.mergeSimilarSegments(
                modeledSegments, intervalThresholdSegmentMean, intervalThresholdMinorAlleleFraction);
        logger.info("Number of segments after smoothing iteration: " + mergedSegments.size());
        currentSegments = new SimpleIntervalCollection(
                new SimpleLocatableMetadata(metadata.getSequenceDictionary()),
                mergedSegments.stream().map(ModeledSegment::getInterval).collect(Collectors.toList()));
        if (doModelFit) {
            fitModel();
        } else {
            modeledSegments.clear();
            modeledSegments.addAll(mergedSegments);
            isModelFit = false;
        }
    }

    /**
     * Writes posterior summaries for the global model parameters to a file.
     */
    public void writeModelParameterFiles(final File copyRatioParameterFile,
                                         final File alleleFractionParameterFile) {
        Utils.nonNull(copyRatioParameterFile);
        Utils.nonNull(alleleFractionParameterFile);
        ensureModelIsFit();
        logger.info("Writing posterior summaries for copy-ratio global parameters to " + copyRatioParameterFile);
        copyRatioModeller.getGlobalParameterDeciles().write(copyRatioParameterFile);
        logger.info("Writing posterior summaries for allele-fraction global parameters to " + alleleFractionParameterFile);
        alleleFractionModeller.getGlobalParameterDeciles().write(alleleFractionParameterFile);
    }

    @VisibleForTesting
    CopyRatioModeller getCopyRatioModeller() {
        return copyRatioModeller;
    }

    @VisibleForTesting
    AlleleFractionModeller getAlleleFractionModeller() {
        return alleleFractionModeller;
    }

    private void ensureModelIsFit() {
        if (!isModelFit) {
            logger.warn("Attempted to write results to file when model was not completely fit. Performing model fit now.");
            fitModel();
        }
    }

    /**
     * Contains private methods for similar-segment merging.
     */
    private static final class SimilarSegmentUtils {
        /**
         * Returns a new, modifiable list of segments with similar segments (i.e., adjacent segments with both
         * segment-mean and minor-allele-fractions posteriors similar; posteriors are similar if the difference between
         * posterior central tendencies is less than intervalThreshold times the posterior credible interval of either summary)
         * merged.  The list of segments is traversed once from beginning to end, and each segment is checked for similarity
         * with the segment to the right and merged until it is no longer similar.
         * @param intervalThresholdSegmentMean         threshold number of credible intervals for segment-mean similarity
         * @param intervalThresholdMinorAlleleFraction threshold number of credible intervals for minor-allele-fraction similarity
         */
        private static List<ModeledSegment> mergeSimilarSegments(final List<ModeledSegment> segments,
                                                                 final double intervalThresholdSegmentMean,
                                                                 final double intervalThresholdMinorAlleleFraction) {
            final List<ModeledSegment> mergedSegments = new ArrayList<>(segments);
            int index = 0;
            while (index < mergedSegments.size() - 1) {
                final ModeledSegment segment1 = mergedSegments.get(index);
                final ModeledSegment segment2 = mergedSegments.get(index + 1);
                if (segment1.getContig().equals(segment2.getContig()) &&
                        areSimilar(segment1, segment2,
                                intervalThresholdSegmentMean, intervalThresholdMinorAlleleFraction)) {
                    mergedSegments.set(index, merge(segment1, segment2));
                    mergedSegments.remove(index + 1);
                    index--; //if merge performed, stay on current segment during next iteration
                }
                index++; //if no merge performed, go to next segment during next iteration
            }
            return mergedSegments;
        }

        //checks similarity of posterior summaries to within a credible-interval threshold;
        //posterior summaries are similar if the difference between posterior central tendencies is less than
        //intervalThreshold times the credible-interval width for both summaries
        private static boolean areSimilar(final ModeledSegment.SimplePosteriorSummary summary1,
                                          final ModeledSegment.SimplePosteriorSummary summary2,
                                          final double intervalThreshold) {
            if (Double.isNaN(summary1.getDecile50()) || Double.isNaN(summary2.getDecile50())) {
                return true;
            }
            final double absoluteDifference = Math.abs(summary1.getDecile50() - summary2.getDecile50());
            return absoluteDifference < intervalThreshold * (summary1.getDecile90() - summary1.getDecile10()) ||
                    absoluteDifference < intervalThreshold * (summary2.getDecile90() - summary2.getDecile10());
        }

        //checks similarity of modeled segments to within credible-interval thresholds for segment mean and minor allele fraction
        private static boolean areSimilar(final ModeledSegment segment1,
                                          final ModeledSegment segment2,
                                          final double intervalThresholdSegmentMean,
                                          final double intervalThresholdMinorAlleleFraction) {
            return areSimilar(segment1.getLog2CopyRatioSimplePosteriorSummary(), segment2.getLog2CopyRatioSimplePosteriorSummary(), intervalThresholdSegmentMean) &&
                    areSimilar(segment1.getMinorAlleleFractionSimplePosteriorSummary(), segment2.getMinorAlleleFractionSimplePosteriorSummary(), intervalThresholdMinorAlleleFraction);
        }

        //merges posterior summaries naively by approximating posteriors as normal
        private static ModeledSegment.SimplePosteriorSummary merge(final ModeledSegment.SimplePosteriorSummary summary1,
                                                                   final ModeledSegment.SimplePosteriorSummary summary2) {
            if (Double.isNaN(summary1.getDecile50()) && !Double.isNaN(summary2.getDecile50())) {
                return summary2;
            }
            if ((!Double.isNaN(summary1.getDecile50()) && Double.isNaN(summary2.getDecile50())) ||
                    (Double.isNaN(summary1.getDecile50()) && Double.isNaN(summary2.getDecile50()))) {
                return summary1;
            }
            //use credible half-interval as standard deviation
            final double standardDeviation1 = (summary1.getDecile90() - summary1.getDecile10()) / 2.;
            final double standardDeviation2 = (summary2.getDecile90() - summary2.getDecile10()) / 2.;
            final double variance = 1. / (1. / Math.pow(standardDeviation1, 2.) + 1. / Math.pow(standardDeviation2, 2.));
            final double mean =
                    (summary1.getDecile50() / Math.pow(standardDeviation1, 2.) + summary2.getDecile50() / Math.pow(standardDeviation2, 2.))
                            * variance;
            final double standardDeviation = Math.sqrt(variance);
            return new ModeledSegment.SimplePosteriorSummary(mean, mean - standardDeviation, mean + standardDeviation);
        }

        private static ModeledSegment merge(final ModeledSegment segment1,
                                            final ModeledSegment segment2) {
            return new ModeledSegment(mergeSegments(segment1.getInterval(), segment2.getInterval()),
                    segment1.getNumPointsCopyRatio() + segment2.getNumPointsCopyRatio(),
                    segment1.getNumPointsAlleleFraction() + segment2.getNumPointsAlleleFraction(),
                    merge(segment1.getLog2CopyRatioSimplePosteriorSummary(), segment2.getLog2CopyRatioSimplePosteriorSummary()),
                    merge(segment1.getMinorAlleleFractionSimplePosteriorSummary(), segment2.getMinorAlleleFractionSimplePosteriorSummary()));
        }

        private static SimpleInterval mergeSegments(final SimpleInterval segment1,
                                                    final SimpleInterval segment2) {
            Utils.validateArg(segment1.getContig().equals(segment2.getContig()),
                    String.format("Cannot join segments %s and %s on different chromosomes.", segment1.toString(), segment2.toString()));
            final int start = Math.min(segment1.getStart(), segment2.getStart());
            final int end = Math.max(segment1.getEnd(), segment2.getEnd());
            return new SimpleInterval(segment1.getContig(), start, end);
        }
    }
}
