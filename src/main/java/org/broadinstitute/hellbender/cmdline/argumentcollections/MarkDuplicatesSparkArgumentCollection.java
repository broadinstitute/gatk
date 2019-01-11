package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.read.markduplicates.MarkDuplicatesScoringStrategy;
import picard.sam.markduplicates.MarkDuplicates;

import java.io.Serializable;


/**
 * An argument collection for use with tools that mark optical
 * duplicates.
 */
public final class MarkDuplicatesSparkArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    public static final String DO_NOT_MARK_UNMAPPED_MATES_LONG_NAME = "do-not-mark-unmapped-mates";
    public static final String DUPLICATE_TAGGING_POLICY_LONG_NAME = "duplicate-tagging-policy";
    public static final String REMOVE_ALL_DUPLICATE_READS = "remove-all-duplicates";
    public static final String REMOVE_SEQUENCING_DUPLICATE_READS = "remove-sequencing-duplicates";

    @Argument(shortName = StandardArgumentDefinitions.DUPLICATE_SCORING_STRATEGY_SHORT_NAME, fullName = StandardArgumentDefinitions.DUPLICATE_SCORING_STRATEGY_LONG_NAME, doc = "The scoring strategy for choosing the non-duplicate among candidates.")
    public MarkDuplicatesScoringStrategy duplicatesScoringStrategy = MarkDuplicatesScoringStrategy.SUM_OF_BASE_QUALITIES;

    @Argument(fullName = MarkDuplicatesSparkArgumentCollection.DO_NOT_MARK_UNMAPPED_MATES_LONG_NAME, doc = "Enabling this option will mean unmapped mates of duplicate marked reads will not be marked as duplicates.")
    public boolean dontMarkUnmappedMates = false;

    @Argument(fullName = MarkDuplicatesSparkArgumentCollection.DUPLICATE_TAGGING_POLICY_LONG_NAME, doc = "Determines how duplicate types are recorded in the DT optional attribute.", optional = true,
              mutex = {REMOVE_ALL_DUPLICATE_READS, REMOVE_SEQUENCING_DUPLICATE_READS})
    public MarkDuplicates.DuplicateTaggingPolicy taggingPolicy = MarkDuplicates.DuplicateTaggingPolicy.DontTag;

    @Argument(fullName = MarkDuplicatesSparkArgumentCollection.REMOVE_ALL_DUPLICATE_READS, doc = "If true do not write duplicates to the output file instead of writing them with appropriate flags set.",
            mutex = {MarkDuplicatesSparkArgumentCollection.DUPLICATE_TAGGING_POLICY_LONG_NAME, MarkDuplicatesSparkArgumentCollection.REMOVE_SEQUENCING_DUPLICATE_READS}, optional = true)
    public boolean removeAllDuplicates = false;

    @Argument(fullName = MarkDuplicatesSparkArgumentCollection.REMOVE_SEQUENCING_DUPLICATE_READS, doc = "If true do not write optical/sequencing duplicates to the output file instead of writing them with appropriate flags set.",
            mutex = {MarkDuplicatesSparkArgumentCollection.DUPLICATE_TAGGING_POLICY_LONG_NAME, MarkDuplicatesSparkArgumentCollection.REMOVE_ALL_DUPLICATE_READS}, optional = true)
    public boolean removeSequencingDuplicates = false;
}