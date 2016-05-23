package org.broadinstitute.hellbender.tools.exome.acnvconversion;

import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.hellbender.tools.exome.*;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.List;
import java.util.stream.Collectors;

/**
 * Utilities for converting ACNV modeled segments.
 */
public final class ACNVModeledSegmentConversionUtils {
    private ACNVModeledSegmentConversionUtils() {}

    /**
     *  Convert our ACNV modeled segments to a generic modeled segments.
     *
     *  <p>Output modeled segments will never have a segment mean of zero, in CR space.  This will be reflected in the
     *  {@code ModeledSegment :: getSegmentMean()} value.</p>
     *
     * @param acnvModeledSegment Never {@code null}
     * @param genome Never {@code null}
     * @return Never {@code null}, but empty list if input list is empty.  Note that the segment mean (in CR space) will
     *  never be zero.
     */
    public static List<ModeledSegment> convertACNVModeledSegmentsToModeledSegments(final List<ACNVModeledSegment> acnvModeledSegment, final Genome genome) {
        Utils.nonNull(acnvModeledSegment);
        Utils.nonNull(genome);
        return acnvModeledSegment.stream().map(s -> convertACNVModeledSegmentToModeledSegment(s, genome)).collect(Collectors.toList());
    }


    private static ModeledSegment convertACNVModeledSegmentToModeledSegment(final ACNVModeledSegment acnvModeledSegment, final Genome genome) {

        final TargetCollection<ReadCountRecord.SingleSampleRecord> targets = genome.getTargets();
        return convertACNVModeledSegmentToModeledSegment(acnvModeledSegment, targets);
    }

    @VisibleForTesting
    static ModeledSegment convertACNVModeledSegmentToModeledSegment(ACNVModeledSegment acnvModeledSegment, TargetCollection<ReadCountRecord.SingleSampleRecord> targets) {

        // Make sure that we do not let segment mean become zero
        final double updatedCenter = Math.max(acnvModeledSegment.getSegmentMeanPosteriorSummary().getCenter(), ParamUtils.log2(TangentNormalizer.EPSILON));

        return new ModeledSegment(acnvModeledSegment.getInterval(), ModeledSegment.NO_CALL,
                targets.targetCount(acnvModeledSegment.getInterval()), updatedCenter);
    }
}
