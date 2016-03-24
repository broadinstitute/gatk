package org.broadinstitute.hellbender.tools.exome.acsconversion;


import org.broadinstitute.hellbender.tools.exome.ACNVModeledSegment;

public class SimpleBalancedSegmentCaller implements BalancedSegmentCaller{

    public SimpleBalancedSegmentCaller() {}

    /**
     * Extremely simple method to classify an ACNVModeledSegment as balanced (i.e. MAF is exactly 0.5) or not.
     *
     * @param segment Never {@code null}
     * @return whether this segment should be considered balanced
     */
    public boolean isSegmentBalanced(final ACNVModeledSegment segment) {
        return segment.getMinorAlleleFractionPosteriorSummary().getUpper() > 0.49;
    }
}
