package org.broadinstitute.hellbender.tools.exome;


import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

/**
 * <p>This caller mimics the legacy ReCapSeg Caller that was originally implemented in ReCapSeg v1.4.5.0.</p>
 *
 * <p>There is a small difference.  The python code was using the same algorithm as intersectBed, which was causing it to drop
 *  the first target of each segment in calculations of the copy neutral targets.  The code here (targets.targets(s)) does
 *  not do this.  This difference in the two codebases can cause a slight difference in the T calculation.  Hence, the
 *  results of this code and the python code will not be exactly the same, but will be
 *  very close.  A fix (to make this code match the python) has been deemed unworthy of our time.</p>
 *
 *  Note: Assumes that the input @{link ReadCountCollection} contains only one sample -- calls will be made for only one
 *        counts column and other columns ignored.
 *
 */
public final class ReCapSegCaller {
    private static final Logger logger = LogManager.getLogger(ReCapSegCaller.class);
    public static String AMPLIFICATION_CALL = "+";
    public static String DELETION_CALL = "-";
    public static String NEUTRAL_CALL = "0";

    //bounds on log_2 coverage for high-confidence neutral segments
    public static final double COPY_NEUTRAL_CUTOFF = 0.1;

    // Number of standard deviations before assuming that a target was an outlier in a segment
    public static final double Z_THRESHOLD = 2;

    private ReCapSegCaller() {} // prevent instantiation

    private static double calculateT(final ReadCountCollection tangentNormalizedCoverage, final List<ModeledSegment> segments) {

        //Get the segments that are likely copy neutral.
        // Math.abs removed to mimic python...
        final List<ModeledSegment> copyNeutralSegments = segments.stream().filter(s -> s.getSegmentMean() < COPY_NEUTRAL_CUTOFF).collect(Collectors.toList());

        // Get the targets that correspond to the copyNeutralSegments... note that individual targets, due to noise,
        //  can be far away from copy neutral
        final TargetCollection<ReadCountRecord.SingleSampleRecord> targetsWithCoverage =
                new HashedListTargetCollection<>(tangentNormalizedCoverage.records().stream().map(ReadCountRecord::asSingleSampleRecord).collect(Collectors.toList()));
        final double[] copyNeutralTargetsCopyRatio = copyNeutralSegments.stream()
                .flatMap(s -> targetsWithCoverage.targets(s).stream())
                .mapToDouble(ReadCountRecord.SingleSampleRecord::getCount).toArray();

        final double meanCopyNeutralTargets = new Mean().evaluate(copyNeutralTargetsCopyRatio);
        final double sigmaCopyNeutralTargets = new StandardDeviation().evaluate(copyNeutralTargetsCopyRatio);

        // Now we filter outliers by only including those w/in 2 standard deviations.
        final double [] filteredCopyNeutralTargetsCopyRatio = Arrays.stream(copyNeutralTargetsCopyRatio)
                .filter(c -> Math.abs(c - meanCopyNeutralTargets) < sigmaCopyNeutralTargets * Z_THRESHOLD).toArray();

        return new StandardDeviation().evaluate(filteredCopyNeutralTargetsCopyRatio);
    }

    /**
     * Make calls for a list of segments based on the coverage data in a set of tangentNormalizedCoverage.
     *  @param tangentNormalizedCoverage the collection representing all targets with their tangent-normalized coverage
     * @param segments segments, each of which holds a reference to these same tangentNormalizedCoverage
     */
    public static List<ModeledSegment> makeCalls(final ReadCountCollection tangentNormalizedCoverage, final List<ModeledSegment> segments) {
        Utils.nonNull(segments, "Can't make calls on a null list of segments.");
        Utils.nonNull(tangentNormalizedCoverage, "Can't make calls on a null list of tangentNormalizedCoverage.");

        final double t = calculateT(tangentNormalizedCoverage, segments);

        logger.info("Running caller that mimics the ReCapSeg 1.4.5.0 (python) caller.");
        logger.info(String.format("T in log2CR: %.4f", t));

        for (final ModeledSegment segment : segments) {

            String call = NEUTRAL_CALL;
            if (segment.getSegmentMean() < -t){
                call = DELETION_CALL;
            }
            if (segment.getSegmentMean() > t){
                call = AMPLIFICATION_CALL;
            }

            // Attach the call to this modeled segment
            segment.setCall(call);
        }

        // Log some information about thresholds chosen for the segments.
        logger.info(String.format("Thresholds Hi, Low: %.4f, %.4f", t, -t));
        logger.info(String.format("Thresholds (CR space) Hi, Low: %.4f, %.4f", Math.pow(2, t), Math.pow(2, -t)));

        return segments;
    }
}
