package org.broadinstitute.hellbender.tools.archive;


import org.apache.commons.math3.distribution.MixtureMultivariateNormalDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.fitting.MultivariateNormalMixtureExpectationMaximization;
import org.apache.commons.math3.exception.ConvergenceException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.linear.SingularMatrixException;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.exome.ModeledSegment;
import org.broadinstitute.hellbender.tools.exome.ReadCountRecord;
import org.broadinstitute.hellbender.tools.exome.SegmentUtils;
import org.broadinstitute.hellbender.tools.exome.TargetCollection;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;


/**
 * Responsible for detecting coverage dropout on select targets.
 * <p>
 * Summary: uses long segments (ones with many targets) to model the distribution of target copy ratios.  If the
 *  distribution looks like two Gaussians with a distance (of means) higher than threshold, assume that there is an
 *  error.
 * </p>
 * <p>Through out this code, GMM = "Gaussian Mixture Model"</p>
 *
 * @author lichtens &lt;lichtens@broadinstitute.org&gt;
 */
public final class CoverageDropoutDetector {

    private static final int NUM_COMPONENTS_GMM = 2;
    private static final String RESULT_COMMENT = "These results were calibrated against ICE with coverage collapsing turned on.  Failures should be manually inspected.";
    public static final int RANDOM_SEED = 337;
    public static final int ABSOLUTE_MIN_NUM_TARGETS = 100;
    private static final Logger logger = LogManager.getLogger(CoverageDropoutDetector.class);

    /**
     * Private helper class to support complex return information regarding the GMM
     *
     * <p>Dev note: Setters have been removed, but should be added if needed later.</p>
     */
    private static final class MixtureMultivariateNormalFitResult {
        private boolean isConverged;
        private MixtureMultivariateNormalDistribution mixtureMultivariateNormalDistribution;
        private boolean wasFitAttempted;

        public MixtureMultivariateNormalFitResult(final MixtureMultivariateNormalDistribution mixtureMultivariateNormalDistribution, final boolean isConverged, final boolean wasFitAttempted) {
            this.isConverged = isConverged;
            this.mixtureMultivariateNormalDistribution = mixtureMultivariateNormalDistribution;
            this.wasFitAttempted = wasFitAttempted;
        }

        public boolean isConverged() {
            return isConverged;
        }

        public MixtureMultivariateNormalDistribution getMixtureMultivariateNormalDistribution() {
            return mixtureMultivariateNormalDistribution;
        }

        public boolean wasFitAttempted() {
            return wasFitAttempted;
        }
    }


    /** Same as {@link CoverageDropoutDetector::retrieveGaussianMixtureModelForFilteredTargets}, but for a single
     * segment.
     *
     * @return never {@code null}.
     */
    private MixtureMultivariateNormalFitResult retrieveGaussianMixtureModelForFilteredTargets(final ModeledSegment segment,
                                                                                              final TargetCollection<ReadCountRecord.SingleSampleRecord> targets, final double minProportion, final int numComponents){
        final List<ModeledSegment> segmentSingleton = Collections.singletonList(segment);
        return retrieveGaussianMixtureModelForFilteredTargets(segmentSingleton, targets, minProportion, numComponents);
    }

    /** <p>Produces a Gaussian mixture model based on the difference between targets and segment means.</p>
     * <p>Filters targets to populations where more than the minProportion lie in a single segment.</p>
     * <p>Returns null if no pass filtering.  Please note that in these cases,
     * in the rest of this class, we use this to assume that a GMM is not a good model.</p>
     *
     * @param segments  -- segments with segment mean in log2 copy ratio space
     * @param targets -- targets with a log2 copy ratio estimate
     * @param minProportion -- minimum proportion of all targets that a given segment must have in order to be used
     *                      in the evaluation
     * @param numComponents -- number of components to use in the GMM.  Usually, this is 2.
     * @return  never {@code null}.  Fitting result with indications whether it converged or was even attempted.
     */
    private MixtureMultivariateNormalFitResult retrieveGaussianMixtureModelForFilteredTargets(final List<ModeledSegment> segments,
                                                                                              final TargetCollection<ReadCountRecord.SingleSampleRecord> targets, double minProportion, int numComponents){

        // For each target in a segment that contains enough targets, normalize the difference against the segment mean
        //  and collapse the filtered targets into the copy ratio estimates.
        final List<Double> filteredTargetsSegDiff = getNumProbeFilteredTargetList(segments, targets, minProportion);

        if (filteredTargetsSegDiff.size() < numComponents) {
            return new MixtureMultivariateNormalFitResult(null, false, false);
        }

        // Assume that Apache Commons wants data points in the first dimension.
        // Note that second dimension of length 2 (instead of 1) is to wrok around funny Apache commons API.
        final double[][] filteredTargetsSegDiff2d = new double[filteredTargetsSegDiff.size()][2];

        // Convert the filtered targets into 2d array (even if second dimension is length 1).  The second dimension is
        //  uncorrelated Gaussian.  This is only to get around funny API in Apache Commons, which will throw an
        //  exception if the length of the second dimension is < 2
        final RandomGenerator rng = RandomGeneratorFactory.createRandomGenerator(new Random(RANDOM_SEED));
        final NormalDistribution nd = new NormalDistribution(rng, 0, .1);
        for (int i = 0; i < filteredTargetsSegDiff.size(); i++) {
            filteredTargetsSegDiff2d[i][0] = filteredTargetsSegDiff.get(i);
            filteredTargetsSegDiff2d[i][1] = nd.sample();
        }

        final MixtureMultivariateNormalDistribution estimateEM0 = MultivariateNormalMixtureExpectationMaximization.estimate(filteredTargetsSegDiff2d, numComponents);
        final MultivariateNormalMixtureExpectationMaximization multivariateNormalMixtureExpectationMaximization = new MultivariateNormalMixtureExpectationMaximization(filteredTargetsSegDiff2d);

        try {
            multivariateNormalMixtureExpectationMaximization.fit(estimateEM0);
        } catch (final MaxCountExceededException | ConvergenceException | SingularMatrixException e) {
            // We are done, we cannot make a fitting.  We should return a result that we attempted a fitting, but it
            //  did not converge.  Include the model as it was when the exception was thrown.
            return new MixtureMultivariateNormalFitResult(multivariateNormalMixtureExpectationMaximization.getFittedModel(), false, true);
        }
        return new MixtureMultivariateNormalFitResult(multivariateNormalMixtureExpectationMaximization.getFittedModel(), true, true);
    }

    /**
     * Retrieve the difference between targets and corresponding segment mean, in CR space.
     *
     * @param segments -- all regions to consider
     * @param targets -- all targets
     * @param minProportion -- minimum proportion (0, 1) of targets in the segment in order for the targets to be considered.
     * @return never {@code null}.  If no targets remain after filtering, an empty list will be returned.
     */
    private List<Double> getNumProbeFilteredTargetList(final List<ModeledSegment> segments,
                                                       final TargetCollection<ReadCountRecord.SingleSampleRecord> targets, final double minProportion) {
        final int allTargetCount = targets.targetCount();
        final double minProportionCount = Math.max(Math.ceil(allTargetCount * minProportion), ABSOLUTE_MIN_NUM_TARGETS);

        // Return a single list of Doubles that contains all targets from the filtered segments.  The collect statement
        //  collapses a list of lists into one list.
        return segments.stream()
                .filter(s -> targets.targetCount(s.getSimpleInterval()) >= minProportionCount)
                .map(s -> SegmentUtils.segmentMeanTargetDifference(s, targets))
                .flatMap(diff -> diff.stream())
                .collect(Collectors.toList());
    }

    /**
     * Using given modeled segments and associated targets, determine whether the given case sample + PoN
     *  combination is appropriate.
     *
     * @param segments -- segments that were produced by a model, Please note that the model must include a segment mean.
     * @param targets -- targets used to create the segments with a copy ratio estimate.
     * @param minTargetProportion -- the minimum proportion (0.0, 1.0) of targets that needs to be included in a segment
     *                      ("num probes") in order to use it for the determination.
     * @param thresholdDistancePerSegment -- distance (of the mean between each component)that must be exceeded to be marked as
     *                          "drop out"
     * @param minProportionGoodSegments -- the proportion of good segments to classify as a "good sample".
     * @param minWeightRestriction -- the minimum weight in the two component gaussian mixture model to be considered
     *                             an actual two component GMM.
     * @return never {@code null}.  Whether this sample should be treated as a failure and additional data for analysis.
     */
    public CoverageDropoutResult determineCoverageDropoutDetected(final List<ModeledSegment> segments, final TargetCollection<ReadCountRecord.SingleSampleRecord> targets,
        final double minTargetProportion, final double thresholdDistancePerSegment, final double minProportionGoodSegments, final double minWeightRestriction) {

        Utils.nonNull(segments, "Segments cannot be null");
        ParamUtils.inRange(minTargetProportion, 0, 1, "minimum target proportion must be in the range (0,1).");
        ParamUtils.inRange(minProportionGoodSegments, 0, 1, "minimum proportion of good segments must be in the range (0,1).");
        ParamUtils.inRange(minWeightRestriction, 0, 1, "minimum component weight must be in the range (0,1).");

        // Prune out the nulls since this indicates that there were not enough targets to make a fitting.
        // Get an individual fitting for each segment remaining.
        final List<MixtureMultivariateNormalFitResult> distros = segments.stream()
                .map(s -> retrieveGaussianMixtureModelForFilteredTargets(s, targets, minTargetProportion, NUM_COMPONENTS_GMM))
                .filter(mm -> mm.wasFitAttempted()).collect(Collectors.toList());

        List<Double> segmentDistances = new ArrayList<>();
        for (MixtureMultivariateNormalFitResult distro: distros) {

            // If we do not meet a minimum weight restriction or we did not converge, then make the distance = 0, since that indicates
            //  that a single Gaussian is the right choice.
            final boolean isWeightLimitMet = distro.getMixtureMultivariateNormalDistribution().getComponents().stream()
                    .allMatch(c -> c.getFirst() > minWeightRestriction);

            if (isWeightLimitMet && distro.isConverged()) {
                final double d = Math.abs(distro.getMixtureMultivariateNormalDistribution().getComponents().get(0).getSecond().getMeans()[0]
                        - distro.getMixtureMultivariateNormalDistribution().getComponents().get(1).getSecond().getMeans()[0]);
                segmentDistances.add(d);
            } else {

                // If we do not converge or weight limit is not met, we know that the 2-GMM is a bad fit, so treat distance
                //  between the two components as zero
                segmentDistances.add(0.0);
            }
        }

        final int numGoodSegments = (int) segmentDistances.stream()
                .filter(d -> d < thresholdDistancePerSegment)
                .count();

        final double proportionGoodSegments = (double) numGoodSegments / (double) segmentDistances.size();
        logger.info("Good segments: " + numGoodSegments + " / " + segmentDistances.size() + "  proportion: " + String.format("%2.2f", proportionGoodSegments));

        return createCoverageDropoutResult(segments, minTargetProportion, thresholdDistancePerSegment, minProportionGoodSegments, minWeightRestriction, numGoodSegments, proportionGoodSegments, segmentDistances.size());
    }

    private CoverageDropoutResult createCoverageDropoutResult(final List<ModeledSegment> segments, final double minTargetProportion, final double thresholdDistancePerSegment, final double minProportionGoodSegments, final double minWeightRestriction, final int numGoodSegments, final double proportionGoodSegments, final long numSegsAfterFiltering) {
        return new CoverageDropoutResult(proportionGoodSegments < minProportionGoodSegments,
                minProportionGoodSegments, minWeightRestriction, minTargetProportion, numGoodSegments, segments.size(),
                numSegsAfterFiltering, thresholdDistancePerSegment, RESULT_COMMENT);
    }
}
