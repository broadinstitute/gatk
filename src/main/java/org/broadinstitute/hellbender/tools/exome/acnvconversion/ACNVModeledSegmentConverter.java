package org.broadinstitute.hellbender.tools.exome.acnvconversion;

import org.broadinstitute.hellbender.tools.exome.ACNVModeledSegment;
import org.broadinstitute.hellbender.tools.exome.Genome;
import org.broadinstitute.hellbender.tools.exome.ModeledSegment;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.List;
import java.util.stream.Collectors;

/**
 * Utilities for converting ACNV modeled segments.
 */
public class ACNVModeledSegmentConverter {
    private ACNVModeledSegmentConverter() {}

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
        return new ModeledSegment(acnvModeledSegment.getInterval(), ModeledSegment.NO_CALL,
                genome.getTargets().targetCount(acnvModeledSegment.getInterval()),
                acnvModeledSegment.getSegmentMeanPosteriorSummary().center());
    }
}
