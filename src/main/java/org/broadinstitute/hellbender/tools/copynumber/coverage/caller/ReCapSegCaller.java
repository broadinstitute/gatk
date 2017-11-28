package org.broadinstitute.hellbender.tools.copynumber.coverage.caller;

import htsjdk.samtools.util.OverlapDetector;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.copynumber.coverage.copyratio.CopyRatio;
import org.broadinstitute.hellbender.tools.copynumber.coverage.copyratio.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.coverage.segmentation.CopyRatioSegment;
import org.broadinstitute.hellbender.tools.copynumber.coverage.segmentation.CopyRatioSegmentCollection;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;

/**
 * <p>This caller mimics the legacy ReCapSeg Caller that was originally implemented in ReCapSeg v1.4.5.0.</p>
 *
 * <p>There is a small difference.  The python code was using the same algorithm as intersectBed, which was causing it to drop
 *  the first interval of each segment in calculations of the copy-neutral intervals.  The code here does
 *  not do this.  This difference in the two codebases can cause a slight difference in the T calculation.  Hence, the
 *  results of this code and the python code will not be exactly the same, but will be
 *  very close.  A fix (to make this code match the python) has been deemed unworthy of our time.</p>
 */
public final class ReCapSegCaller {
    private static final Logger logger = LogManager.getLogger(ReCapSegCaller.class);

    //bounds on log_2 coverage for high-confidence neutral segments
    private static final double COPY_NEUTRAL_CUTOFF = 0.1;
    // Number of standard deviations before assuming that an interval was an outlier in a segment
    private static final double Z_THRESHOLD = 2;

    private final CopyRatioSegmentCollection copyRatioSegments;
    private final LinkedHashMap<CopyRatioSegment, Set<CopyRatio>> segmentToCopyRatiosMap;

    /**
     * @param denoisedCopyRatios in log2 space
     */
    public ReCapSegCaller(final CopyRatioCollection denoisedCopyRatios,
                          final CopyRatioSegmentCollection copyRatioSegments) {
        this.copyRatioSegments = Utils.nonNull(copyRatioSegments);
        Utils.validateArg(denoisedCopyRatios.getSampleName().equals(copyRatioSegments.getSampleName()),
                "Denoised copy ratios and copy-ratio segments do not have the same sample name.");
        segmentToCopyRatiosMap = constructSegmentToCopyRatiosMap(denoisedCopyRatios, copyRatioSegments);
    }

    private static LinkedHashMap<CopyRatioSegment,Set<CopyRatio>> constructSegmentToCopyRatiosMap(final CopyRatioCollection denoisedCopyRatios,
                                                                                                  final CopyRatioSegmentCollection copyRatioSegments) {
        final LinkedHashMap<CopyRatioSegment, Set<CopyRatio>> segmentToCopyRatiosMap = new LinkedHashMap<>();
        final OverlapDetector<CopyRatio> copyRatioMidpointOverlapDetector = denoisedCopyRatios.getMidpointOverlapDetector();
        for (final CopyRatioSegment segment : copyRatioSegments.getRecords()) {
            final int numPointsExpected = segment.getNumPoints();
            final Set<CopyRatio> copyRatiosInSegment = copyRatioMidpointOverlapDetector.getOverlaps(segment);
            if (copyRatiosInSegment.size() != numPointsExpected) {
                throw new IllegalArgumentException("Denoised copy ratios and copy-ratio segments are not consistent.");
            }
            segmentToCopyRatiosMap.put(segment, copyRatiosInSegment);
        }
        return segmentToCopyRatiosMap;
    }

    private double calculateT() {
        //Get the segments that are likely copy neutral.
        //Math.abs removed to mimic python...
        final List<CopyRatioSegment> copyNeutralSegments = segmentToCopyRatiosMap.keySet().stream()
                .filter(s -> s.getMeanLog2CopyRatio() < COPY_NEUTRAL_CUTOFF).collect(Collectors.toList());

        //Get the intervals that correspond to the copyNeutralSegments... note that individual intervals, due to noise,
        //can be far away from copy neutral
        final double[] copyNeutralIntervals = copyNeutralSegments.stream()
                .flatMap(s -> segmentToCopyRatiosMap.get(s).stream())
                .mapToDouble(CopyRatio::getLog2CopyRatioValue).toArray();

        final double meanCopyNeutralIntervals = new Mean().evaluate(copyNeutralIntervals);
        final double sigmaCopyNeutralIntervals = new StandardDeviation().evaluate(copyNeutralIntervals);

        // Now we filter outliers by only including those w/in 2 standard deviations.
        final double [] filteredCopyNeutralIntervals = Arrays.stream(copyNeutralIntervals)
                .filter(c -> Math.abs(c - meanCopyNeutralIntervals) < sigmaCopyNeutralIntervals * Z_THRESHOLD).toArray();

        return new StandardDeviation().evaluate(filteredCopyNeutralIntervals);
    }

    public CalledCopyRatioSegmentCollection makeCalls() {
        final double t = calculateT();

        logger.info("Running caller that mimics the ReCapSeg 1.4.5.0 (python) caller.");
        // Log some information about thresholds chosen for the segments.
        logger.info(String.format("Copy neutral (log2CR space) [%.4f, %.4f]", -t, t));
        logger.info(String.format("Copy neutral (CR space) [%.4f, %.4f]", Math.pow(2, -t), Math.pow(2, t)));

        final Set<CopyRatioSegment> segments = segmentToCopyRatiosMap.keySet();
        final List<CalledCopyRatioSegment> calledSegments = new ArrayList<>(segments.size());
        for (final CopyRatioSegment segment : segments) {
            if (segment.getMeanLog2CopyRatio() < -t) {
                calledSegments.add(new CalledCopyRatioSegment(segment, CalledCopyRatioSegment.Call.DELETION));
            } else if (segment.getMeanLog2CopyRatio() > t) {
                calledSegments.add(new CalledCopyRatioSegment(segment, CalledCopyRatioSegment.Call.AMPLIFICATION));
            } else {
                calledSegments.add(new CalledCopyRatioSegment(segment, CalledCopyRatioSegment.Call.NEUTRAL));
            }
        }

        return new CalledCopyRatioSegmentCollection(copyRatioSegments.getSampleMetadata(), calledSegments);
    }
}
