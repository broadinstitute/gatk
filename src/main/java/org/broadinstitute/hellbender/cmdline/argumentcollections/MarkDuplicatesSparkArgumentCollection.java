package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.tools.spark.transforms.markduplicates.MarkDuplicatesSpark;
import org.broadinstitute.hellbender.utils.read.markduplicates.MarkDuplicatesScoringStrategy;
import org.broadinstitute.hellbender.utils.read.markduplicates.OpticalDuplicateFinder;

import java.io.Serializable;


/**
 * An argument collection for use with tools that mark optical
 * duplicates.
 */
public final class MarkDuplicatesSparkArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    @Argument(shortName = "DS", fullName = "DUPLICATE_SCORING_STRATEGY", doc = "The scoring strategy for choosing the non-duplicate among candidates.")
    public MarkDuplicatesScoringStrategy duplicatesScoringStrategy = MarkDuplicatesScoringStrategy.SUM_OF_BASE_QUALITIES;

    @Argument(fullName = MarkDuplicatesSpark.DO_NOT_MARK_UNMAPPED_MATES, doc = "Enabling this option will mean unmapped mates of duplicate marked reads will not be marked as duplicates.")
    public boolean dontMarkUnmappedMates = false;
}