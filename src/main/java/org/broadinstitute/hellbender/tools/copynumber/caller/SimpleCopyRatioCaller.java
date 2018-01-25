package org.broadinstitute.hellbender.tools.copynumber.caller;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.copynumber.formats.CopyNumberFormatsUtils;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CalledCopyRatioSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CopyRatioSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CalledCopyRatioSegment;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyRatioSegment;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * This caller is loosely based on the legacy ReCapSeg caller that was originally implemented in ReCapSeg v1.4.5.0,
 * but introduces major changes.  The method is as follows:
 * 1) use the non-log2 mean copy ratio to determine copy-neutral segments,
 * 2) weight segments by length for determining the mean and standard deviation of the non-log2 copy ratio in copy-neutral segments,
 * 3) filter outlier copy-neutral segments by non-log2 copy ratio z-score,
 * 4) use the filtered copy-neutral segments to determine a length-weighted mean and standard deviation,
 * 5) call segments using z-score based on this mean and standard deviation.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class SimpleCopyRatioCaller {
    private static final Logger logger = LogManager.getLogger(SimpleCopyRatioCaller.class);

    private final double neutralSegmentCopyRatioLowerBound;
    private final double neutralSegmentCopyRatioUpperBound;
    private final double outlierNeutralSegmentCopyRatioZScoreThreshold;
    private final double callingCopyRatioZScoreThreshold;
    private final Statistics callingStatistics;

    private final CopyRatioSegmentCollection copyRatioSegments;

    /**
     * @param neutralSegmentCopyRatioLowerBound             non-log2 copy ratio must be within [lower bound, upper bound] for a segment to be copy neutral
     * @param neutralSegmentCopyRatioUpperBound             non-log2 copy ratio must be within [lower bound, upper bound] for a segment to be copy neutral
     * @param outlierNeutralSegmentCopyRatioZScoreThreshold z-score on non-log2 copy ratio above which a copy-neutral segment is assumed to be an outlier
     *                                                      and not included in the calculation of the length-weighted standard deviation of
     *                                                      non-log2 copy ratio in copy-neutral segments
     * @param callingCopyRatioZScoreThreshold               z-score with respect to length-weighted standard deviation of non-log2 copy ratio
     *                                                      in non-outlier copy-neutral segments used for calling segments
     */
    public SimpleCopyRatioCaller(final CopyRatioSegmentCollection copyRatioSegments,
                                 final double neutralSegmentCopyRatioLowerBound,
                                 final double neutralSegmentCopyRatioUpperBound,
                                 final double outlierNeutralSegmentCopyRatioZScoreThreshold,
                                 final double callingCopyRatioZScoreThreshold) {
        ParamUtils.isPositiveOrZero(neutralSegmentCopyRatioLowerBound, "Copy-neutral lower bound must be non-negative.");
        Utils.validateArg(neutralSegmentCopyRatioLowerBound < neutralSegmentCopyRatioUpperBound, "Copy-neutral lower bound must be less than upper bound.");
        ParamUtils.isPositive(outlierNeutralSegmentCopyRatioZScoreThreshold, "Outlier z-score threshold must be positive.");
        ParamUtils.isPositive(callingCopyRatioZScoreThreshold, "Calling z-score threshold must be positive.");
        this.copyRatioSegments = Utils.nonNull(copyRatioSegments);
        this.neutralSegmentCopyRatioLowerBound = neutralSegmentCopyRatioLowerBound;
        this.neutralSegmentCopyRatioUpperBound = neutralSegmentCopyRatioUpperBound;
        this.outlierNeutralSegmentCopyRatioZScoreThreshold = outlierNeutralSegmentCopyRatioZScoreThreshold;
        this.callingCopyRatioZScoreThreshold = callingCopyRatioZScoreThreshold;
        callingStatistics = calculateCallingStatistics();
    }

    public CalledCopyRatioSegmentCollection makeCalls() {
        final List<CopyRatioSegment> segments = copyRatioSegments.getRecords();
        final List<CalledCopyRatioSegment> calledSegments = new ArrayList<>(segments.size());
        for (final CopyRatioSegment segment : segments) {
            final double copyRatioMean = Math.pow(2., segment.getMeanLog2CopyRatio());
            if (neutralSegmentCopyRatioLowerBound <= copyRatioMean && copyRatioMean <= neutralSegmentCopyRatioUpperBound) {
                calledSegments.add(new CalledCopyRatioSegment(segment, CalledCopyRatioSegment.Call.NEUTRAL));
            } else {
                final double copyRatioDeviation = copyRatioMean - callingStatistics.mean;
                if (copyRatioDeviation < -callingStatistics.standardDeviation * callingCopyRatioZScoreThreshold) {
                    calledSegments.add(new CalledCopyRatioSegment(segment, CalledCopyRatioSegment.Call.DELETION));
                } else if (copyRatioDeviation > callingStatistics.standardDeviation * callingCopyRatioZScoreThreshold) {
                    calledSegments.add(new CalledCopyRatioSegment(segment, CalledCopyRatioSegment.Call.AMPLIFICATION));
                } else {
                    calledSegments.add(new CalledCopyRatioSegment(segment, CalledCopyRatioSegment.Call.NEUTRAL));
                }
            }
        }
        return new CalledCopyRatioSegmentCollection(copyRatioSegments.getMetadata(), calledSegments);
    }

    private Statistics calculateCallingStatistics() {
        //get the segments that fall within the copy-neutral region
        final List<CopyRatioSegment> copyNeutralSegments = copyRatioSegments.getRecords().stream()
                .filter(s -> {
                    final double copyRatioMean = Math.pow(2., s.getMeanLog2CopyRatio());
                    return neutralSegmentCopyRatioLowerBound <= copyRatioMean && copyRatioMean <= neutralSegmentCopyRatioUpperBound;})
                .collect(Collectors.toList());
        logger.info(String.format("%d segments in copy-neutral region [%s, %s]...", copyNeutralSegments.size(),
                CopyNumberFormatsUtils.formatDouble(neutralSegmentCopyRatioLowerBound),
                CopyNumberFormatsUtils.formatDouble(neutralSegmentCopyRatioUpperBound)));

        //calculate length-weighted statistics of unfiltered copy-neutral segments
        final Statistics unfilteredStatistics = calculateLengthWeightedStatistics(copyNeutralSegments);
        logger.info(String.format("Length-weighted mean of segments in copy-neutral region (CR space): %s",
                CopyNumberFormatsUtils.formatDouble(unfilteredStatistics.mean)));
        logger.info(String.format("Length-weighted standard deviation for segments in copy-neutral region : %s",
                CopyNumberFormatsUtils.formatDouble(unfilteredStatistics.standardDeviation)));

        //filter outlier segments by only including those within 2 standard deviations
        final List<CopyRatioSegment> filteredCopyNeutralSegments = copyNeutralSegments.stream()
                .filter(s -> Math.abs(Math.pow(2., s.getMeanLog2CopyRatio()) - unfilteredStatistics.mean)
                        <= unfilteredStatistics.standardDeviation * outlierNeutralSegmentCopyRatioZScoreThreshold)
                .collect(Collectors.toList());
        logger.info(String.format("%d / %d segments in copy-neutral region remain after outliers filtered using z-score threshold (%s)...",
                filteredCopyNeutralSegments.size(), copyNeutralSegments.size(),
                CopyNumberFormatsUtils.formatDouble(outlierNeutralSegmentCopyRatioZScoreThreshold)));

        final Statistics statistics = calculateLengthWeightedStatistics(filteredCopyNeutralSegments);
        logger.info(String.format("Length-weighted mean for z-score calling (CR space): %s",
                CopyNumberFormatsUtils.formatDouble(statistics.mean)));
        logger.info(String.format("Length-weighted standard deviation for z-score calling (CR space): %s",
                CopyNumberFormatsUtils.formatDouble(statistics.standardDeviation)));

        return statistics;
    }

    private static Statistics calculateLengthWeightedStatistics(final List<CopyRatioSegment> copyRatioSegments) {
        final List<Integer> segmentLengths = copyRatioSegments.stream()
                .map(c -> c.getInterval().getLengthOnReference())
                .collect(Collectors.toList());
        final double totalLength = segmentLengths.stream().mapToDouble(Integer::doubleValue).sum();
        final int numSegments = segmentLengths.size();
        final double lengthWeightedCopyRatioMean = IntStream.range(0, numSegments)
                .mapToDouble(i -> segmentLengths.get(i) * Math.pow(2., copyRatioSegments.get(i).getMeanLog2CopyRatio()))
                .sum() / totalLength;
        final double lengthWeightedCopyRatioStandardDeviation = Math.sqrt(IntStream.range(0, numSegments)
                .mapToDouble(i -> segmentLengths.get(i) * Math.pow(Math.pow(2., copyRatioSegments.get(i).getMeanLog2CopyRatio()) - lengthWeightedCopyRatioMean, 2))
                .sum() / (((double) (numSegments - 1) / numSegments) * totalLength));
        return new Statistics(lengthWeightedCopyRatioMean, lengthWeightedCopyRatioStandardDeviation);
    }

    private static final class Statistics {
        private final double mean;
        private final double standardDeviation;

        private Statistics(final double mean,
                           final double standardDeviation) {
            this.mean = mean;
            this.standardDeviation = standardDeviation;
        }
    }
}
