package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.read.markduplicates.MarkDuplicatesScoringStrategy;

import java.io.Serializable;


/**
 * An argument collection for use with tools that mark optical
 * duplicates.
 */
public final class MarkDuplicatesSparkArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    public static final String DO_NOT_MARK_UNMAPPED_MATES_LONG_NAME = "do-not-mark-unmapped-mates";

    @Argument(shortName = StandardArgumentDefinitions.DUPLICATE_SCORING_STRATEGY_SHORT_NAME, fullName = StandardArgumentDefinitions.DUPLICATE_SCORING_STRATEGY_LONG_NAME, doc = "The scoring strategy for choosing the non-duplicate among candidates.")
    public MarkDuplicatesScoringStrategy duplicatesScoringStrategy = MarkDuplicatesScoringStrategy.SUM_OF_BASE_QUALITIES;

    @Argument(fullName = MarkDuplicatesSparkArgumentCollection.DO_NOT_MARK_UNMAPPED_MATES_LONG_NAME, doc = "Enabling this option will mean unmapped mates of duplicate marked reads will not be marked as duplicates.")
    public boolean dontMarkUnmappedMates = false;
}