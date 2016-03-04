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
     * @param acnvModeledSegment Never {@code null}
     * @param genome Never {@code null}
     * @return Never {@code null}, but empty list if input list is empty
     */
    public static List<ModeledSegment> convertACNVModeledSegmentsToModeledSegments(final List<ACNVModeledSegment> acnvModeledSegment, final Genome genome) {
        Utils.nonNull(acnvModeledSegment);
        Utils.nonNull(genome);
        return acnvModeledSegment.stream().map(s -> convertACNVModeledSegmentToModeledSegment(s, genome)).collect(Collectors.toList());
    }


    private static ModeledSegment convertACNVModeledSegmentToModeledSegment(final ACNVModeledSegment acnvModeledSegment, final Genome genome) {

        final TargetCollection<TargetCoverage> targets = genome.getTargets();
        return convertACNVModeledSegmentToModeledSegment(acnvModeledSegment, targets);
    }

    @VisibleForTesting
    static ModeledSegment convertACNVModeledSegmentToModeledSegment(ACNVModeledSegment acnvModeledSegment, TargetCollection<TargetCoverage> targets) {

        // Make sure that we do not let segment mean become zero
        double updatedCenter = acnvModeledSegment.getSegmentMeanPosteriorSummary().center();
        if (Math.pow(2, updatedCenter) <= 0) {
            updatedCenter = ParamUtils.log2(TangentNormalizer.EPSILON);
        }

        return new ModeledSegment(acnvModeledSegment.getInterval(), ModeledSegment.NO_CALL,
                targets.targetCount(acnvModeledSegment.getInterval()), updatedCenter);
    }
}
