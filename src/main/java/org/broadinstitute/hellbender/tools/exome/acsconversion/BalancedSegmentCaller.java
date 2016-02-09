package org.broadinstitute.hellbender.tools.exome.acsconversion;

import org.broadinstitute.hellbender.tools.exome.ACNVModeledSegment;


public interface BalancedSegmentCaller {
    boolean isSegmentBalanced(final ACNVModeledSegment segment);
}
