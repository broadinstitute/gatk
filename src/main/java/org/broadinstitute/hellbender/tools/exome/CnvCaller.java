package org.broadinstitute.hellbender.tools.exome;


import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.Arrays;
import java.util.List;

/**
 * This caller is considered experimental.
 *
 * Simple caller that
 * 1) estimates the value of normalized coverage corresponding to copy neutral segments
 * This value is not necessarily equal to zero.  For example, a large amplification on one segment
 * depresses the proportional coverage of all other segments.
 * 2) estimates the variance of targets' normalized coverage within a segment, assuming that
 * this variance is globally shared by all segments
 * 3) calls each segment based on a p-value for the null hypothesis that its targets'
 * normalized coverages are drawn from a normal distribution defined by the mean from 1)
 * and the variance from 2)
 *
 * Note:  Assumes that the incoming target file is in log_2(CR)
 * Note:  Assumes that the incoming segment means are CR.
 *
 * @author David Benjamin
 */
public final class CnvCaller {
    private static final Logger logger = LogManager.getLogger(CnvCaller.class);
    public static String AMPLIFICATION_CALL = "+";
    public static String DELETION_CALL = "-";
    public static String NEUTRAL_CALL = "0";

    // default number of standard deviations from mean required to call an amplification or deletion
    public static final double DEFAULT_Z_SCORE_THRESHOLD = 2.0;

    //bounds on log_2 coverage we take as non-neutral a priori
    public static final double NON_NEUTRAL_THRESHOLD = 0.3;

    //bounds on log_2 coverage for high-confidence neutral segments
    public static final double COPY_NEUTRAL_CUTOFF = 0.1;

    private CnvCaller() {} // prevent instantiation

    /**
     * Make calls for a list of segments based on the coverage data in a set of targets.
     *
     * @param targets the collection representing all targets
     * @param segments segments, each of which holds a reference to these same targets
     * @param zThreshold number of standard deviations from mean required to call an amplification or deletion
     */
    public static List<ModeledSegment> makeCalls(final TargetCollection<TargetCoverage> targets, final List<ModeledSegment> segments, final double zThreshold) {
        Utils.nonNull(segments, "Can't make calls on a null list of segments.");
        Utils.nonNull(targets, "Can't make calls on a null list of targets.");
        ParamUtils.isPositiveOrZero(zThreshold, "zThreshold must be positive or zero.");

        logger.warn("This caller is still considered experimental.");

        /**
         * estimate the copy-number neutral log_2 coverage as the median over all targets.
         * As a precaution we take the median after filtering extreme coverages that are definitely not neutral
         * i.e. those below -NON_NEUTRAL_COVERAGE or above +NON_NEUTRAL_COVERAGE
         */
        final double []  nonExtremeCoverage = targets.targets()
                .stream()
                .mapToDouble(TargetCoverage::getCoverage)
                .filter(c -> Math.abs(c) < NON_NEUTRAL_THRESHOLD)
                .toArray();
        final double neutralCoverage = new Median().evaluate(nonExtremeCoverage);

        logger.info(String.format("True copy neutral estimate (median): %.3f", neutralCoverage));
        logger.info(String.format("True copy neutral estimate (median) in CR: %.3f", Math.pow(2, neutralCoverage)));

        // Get the standard deviation of targets belonging to high-confidence neutral segments
        final double[] nearNeutralCoverage = segments.stream()
                .filter(s -> Math.abs(s.getSegmentMean() - neutralCoverage) < COPY_NEUTRAL_CUTOFF)
                .flatMap(s -> targets.targets(s).stream())
                .mapToDouble(TargetCoverage::getCoverage)
                .toArray();
        final double neutralSigma = new StandardDeviation().evaluate(nearNeutralCoverage);
        logger.info(String.format("Neutral Sigma (unfiltered): %.4f", neutralSigma));

        final double filterThreshold = new Percentile(2).evaluate(nearNeutralCoverage);
        final double[] nearNeutralCoverageFiltered = Arrays.stream(nearNeutralCoverage)
                .filter(c -> c > filterThreshold)
                .toArray();
        final double neutralSigmaFiltered = new StandardDeviation().evaluate(nearNeutralCoverageFiltered);
        logger.info(String.format("Neutral Sigma (filtered -- keep > %.4f): %.4f", filterThreshold, neutralSigmaFiltered));


        for (final ModeledSegment segment : segments) {
            //Get the number of standard deviations the mean falls from the copy neutral value
            //we should really use the standard deviation of the mean of segment.numTargets() targets
            //which is sigma/ sqrt(numTargets).  However, the underlying segment means vary due to noise not
            //removed by tangent normalization, and such a caller would be too strict.
            //Thus we just use the population standard deviation of targets, and defer a mathematically
            //principled treatment for a later caller
            final double z = (segment.getSegmentMean() - neutralCoverage)/neutralSigmaFiltered;

            final String call = z < -zThreshold ? DELETION_CALL: (z > zThreshold ? AMPLIFICATION_CALL : NEUTRAL_CALL);

            // Attach the call to this modeled segment
            segment.setCall(call);
        }

        // Log some information about thresholds chosen for the segments.
        final double zThresholdInLog2SpaceHigh = zThreshold * neutralSigmaFiltered + neutralCoverage;
        final double zThresholdInLog2SpaceLow = -zThreshold * neutralSigmaFiltered + neutralCoverage;
        logger.info(String.format("Thresholds Hi, Low: %.4f, %.4f", zThresholdInLog2SpaceHigh, zThresholdInLog2SpaceLow));
        logger.info(String.format("Thresholds (CR space) Hi, Low: %.4f, %.4f", Math.pow(2, zThresholdInLog2SpaceHigh), Math.pow(2, zThresholdInLog2SpaceLow)));

        return segments;
    }
}
