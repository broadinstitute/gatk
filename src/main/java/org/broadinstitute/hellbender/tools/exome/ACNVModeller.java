package org.broadinstitute.hellbender.tools.exome;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionModeller;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.mcmc.PosteriorSummary;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Represents an ACNV segmented model for copy ratio and allele fraction.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class ACNVModeller {
    protected static final String INITIAL_SEG_FILE_TAG = "sim-0";
    protected static final String INTERMEDIATE_SEG_FILE_TAG = "sim";
    protected static final String FINAL_SEG_FILE_TAG = "sim-final";

    private static final int MAX_SIMILAR_SEGMENT_MERGE_ITERATIONS = 25;

    //use 95% HPD interval to construct {@link PosteriorSummary} for segment means and minor allele fractions
    private static final double CREDIBLE_INTERVAL_ALPHA = 0.05;

    public static final Logger logger = LogManager.getLogger(ACNVModeller.class);

    private SegmentedModel segmentedModel;
    private final List<ACNVModeledSegment> segments = new ArrayList<>();

    private final String outputPrefix;

    private final int numSamplesCopyRatio;
    private final int numBurnInCopyRatio;
    private final int numSamplesAlleleFraction;
    private final int numBurnInAlleleFraction;
    private final JavaSparkContext ctx;

    /**
     * Constructs a copy-ratio and allele-fraction modeller for a {@link SegmentedModel},
     * specifying number of total samples and number of burn-in samples for Markov-Chain Monte Carlo model fitting.
     * An initial model fit is performed.
     *
     * @param segmentedModel            contains segments, target coverages, and SNP counts to model
     * @param outputPrefix              output prefix string (file directory + name prefix) for similar-segment merging
     * @param numSamplesCopyRatio       number of total samples for copy-ratio model MCMC
     * @param numBurnInCopyRatio        number of burn-in samples to discard for copy-ratio model MCMC
     * @param numSamplesAlleleFraction  number of total samples for allele-fraction model MCMC
     * @param numBurnInAlleleFraction   number of burn-in samples to discard for allele-fraction model MCMC
     * @param ctx                       JavaSparkContext, used for kernel density estimation in {@link PosteriorSummary}
     */
    public ACNVModeller(final SegmentedModel segmentedModel, final String outputPrefix,
                        final int numSamplesCopyRatio, final int numBurnInCopyRatio,
                        final int numSamplesAlleleFraction, final int numBurnInAlleleFraction,
                        final JavaSparkContext ctx) {
        this.segmentedModel = segmentedModel;
        this.outputPrefix = outputPrefix;
        this.numSamplesCopyRatio = numSamplesCopyRatio;
        this.numBurnInCopyRatio = numBurnInCopyRatio;
        this.numSamplesAlleleFraction = numSamplesAlleleFraction;
        this.numBurnInAlleleFraction = numBurnInAlleleFraction;
        this.ctx = ctx;
        fitModel(numSamplesCopyRatio, numBurnInCopyRatio, numSamplesAlleleFraction, numBurnInAlleleFraction);
    }

    /**
     * Performs similar-segment merging on the list of {@link ACNVModeledSegment} held internally until convergence,
     * specifying number of total samples and number of burn-in samples for Markov-Chain Monte Carlo model fitting
     * performed after each iteration.
     * <p>
     *     First, the list of {@link ACNVModeledSegment} (i.e., the list of {@link PosteriorSummary} for each
     *     segment copy-ratio mean and minor allele fraction) is taken from the initial model fit.
     *     Then, {@link SegmentMergeUtils#mergeSimilarSegments} is used to merge segments in this
     *     list given specified copy-ratio and minor-allele-fraction merging thresholds, after which another MCMC
     *     fit is performed; this is repeated until convergence (i.e., no more similar segments are found).
     * </p>
     * Segment files (with filenames specified by {@link ACNVModeller#outputPrefix} are output for the initial
     * model fit and the model fits after each merge iteration as a side effect.
     *
     * @param intervalThresholdSegmentMean         threshold number of credible intervals for segment-mean similarity
     * @param intervalThresholdMinorAlleleFraction threshold number of credible intervals for minor-allele-fraction similarity
     * @param numSamplesCopyRatio       number of total samples for copy-ratio model MCMC (for each iteration)
     * @param numBurnInCopyRatio        number of burn-in samples to discard for copy-ratio model MCMC (for each iteration)
     * @param numSamplesAlleleFraction  number of total samples for allele-fraction model MCMC (for each iteration)
     * @param numBurnInAlleleFraction   number of burn-in samples to discard for allele-fraction model MCMC (for each iteration)
     */
    public void mergeSimilarSegments(final double intervalThresholdSegmentMean,
                                     final double intervalThresholdMinorAlleleFraction,
                                     final int numSamplesCopyRatio, final int numBurnInCopyRatio,
                                     final int numSamplesAlleleFraction, final int numBurnInAlleleFraction) {
        logger.info("Starting similar-segment merging...");
        logger.info("Initial number of segments: " + segments.size());
        List<ACNVModeledSegment> mergedSegments = new ArrayList<>(segments);

        //write initial model fit to file
        final File initialModeledSegmentsFile = new File(outputPrefix + "-" + INITIAL_SEG_FILE_TAG + ".seg");
        SegmentUtils.writeACNVModeledSegmentFile(initialModeledSegmentsFile, mergedSegments, segmentedModel.getGenome());

        //perform iterations of similar-segment merging until all similar segments are merged
        int prevNumSegments;
        for (int numIterations = 1; numIterations <= MAX_SIMILAR_SEGMENT_MERGE_ITERATIONS; numIterations++) {
            logger.info("Similar-segment merging iteration: " + numIterations);
            logger.info("Number of segments before merging: " + mergedSegments.size());
            prevNumSegments = mergedSegments.size();
            mergedSegments =
                    SegmentMergeUtils.mergeSimilarSegments(segments,
                            intervalThresholdSegmentMean, intervalThresholdMinorAlleleFraction);
            logger.info("Number of segments after merging: " + mergedSegments.size());
            segmentedModel = new SegmentedModel(toUnmodeledSegments(mergedSegments), segmentedModel.getGenome());
            if (mergedSegments.size() == prevNumSegments) {
                break;
            }
            //refit model and write to file after each iteration
            fitModel(numSamplesCopyRatio, numBurnInCopyRatio, numSamplesAlleleFraction, numBurnInAlleleFraction);
            final File modeledSegmentsFile = new File(outputPrefix + "-" + INTERMEDIATE_SEG_FILE_TAG + "-" + numIterations + ".seg");
            SegmentUtils.writeACNVModeledSegmentFile(modeledSegmentsFile, segments, segmentedModel.getGenome());
        }

        //perform final model fit and write to file
        fitModel(this.numSamplesCopyRatio, this.numBurnInCopyRatio, this.numSamplesAlleleFraction, this.numBurnInAlleleFraction);
        final File finalModeledSegmentsFile = new File(outputPrefix + "-" + FINAL_SEG_FILE_TAG + ".seg");
        SegmentUtils.writeACNVModeledSegmentFile(finalModeledSegmentsFile, segments, segmentedModel.getGenome());
    }

    //fits copy-ratio and allele-fraction models using MCMC and reinitializes internal list of ACNVModeledSegments
    private void fitModel(final int numSamplesCopyRatio, final int numBurnInCopyRatio,
                          final int numSamplesAlleleFraction, final int numBurnInAlleleFraction) {
        logger.info("Fitting copy-ratio model...");
        final ACNVCopyRatioModeller copyRatioModeller = new ACNVCopyRatioModeller(segmentedModel);
        copyRatioModeller.fitMCMC(numSamplesCopyRatio, numBurnInCopyRatio);
        logger.info("Fitting allele-fraction model...");
        final AlleleFractionModeller alleleFractionModeller = new AlleleFractionModeller(segmentedModel);
        alleleFractionModeller.fitMCMC(numSamplesAlleleFraction, numBurnInAlleleFraction);

        segments.clear();
        final List<SimpleInterval> unmodeledSegments = copyRatioModeller.getSegmentedModel().getSegments();
        final List<PosteriorSummary> segmentMeansPosteriorSummaries =
                copyRatioModeller.getSegmentMeansPosteriorSummaries(CREDIBLE_INTERVAL_ALPHA, ctx);
        final List<PosteriorSummary> minorAlleleFractionsPosteriorSummaries =
                alleleFractionModeller.getMinorAlleleFractionsPosteriorSummaries(CREDIBLE_INTERVAL_ALPHA, ctx);
        for (int segment = 0; segment < unmodeledSegments.size(); segment++) {
            segments.add(new ACNVModeledSegment(unmodeledSegments.get(segment),
                    segmentMeansPosteriorSummaries.get(segment), minorAlleleFractionsPosteriorSummaries.get(segment)));
        }
    }

    //converts list of ACNVModeledSegments to list of SimpleIntervals
    private List<SimpleInterval> toUnmodeledSegments(final List<ACNVModeledSegment> segments) {
        final List<SimpleInterval> unmodeledSegments = new ArrayList<>();
        for (final ACNVModeledSegment segment : segments) {
            unmodeledSegments.add(segment.getInterval());
        }
        return unmodeledSegments;
    }
}
