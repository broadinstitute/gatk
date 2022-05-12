package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.barclay.argparser.Advanced;
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

    public static final String FLOW_MD_MODE_LONG_NAME = "flowbased";

    public static final String FLOW_QUALITY_SUM_STRATEGY_LONG_NAME = "flow-quality-sum-strategy";
    public static final String SINGLE_END_READS_END_POSITION_SIGNIFICANT = "single-end-reads-end-position-significant";
    public static final String FLOW_END_POS_UNCERTAINTY_LONG_NAME = "flow-end-pos-uncertainty";
    public static final String SINGLE_END_READS_CLIPPING_IS_END_LONG_NAME = "single-end-reads-clipping-is-end";
    public static final String FLOW_SKIP_START_HOMOPOLYMERS_LONG_NAME = "flow-skip-start-homopolymers";
    public static final String FLOW_Q_IS_KNOWN_END_LONG_NAME = "flow-q-is-known-end";

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

    @Advanced
    @Argument(fullName = FLOW_QUALITY_SUM_STRATEGY_LONG_NAME, doc = "Use specific quality summing strategy for flow based reads. The strategy ensures that the same " +
            "(and correct) quality value is used for all bases of the same homopolymer. Default false.", optional = true)
    public boolean FLOW_QUALITY_SUM_STRATEGY = false;

    @Advanced
    @Argument(fullName = SINGLE_END_READS_END_POSITION_SIGNIFICANT, doc = "Make end location of read (fragment) be significant when considering duplicates, " +
            "in addition to the start location, which is always significant (should only be applied to flow based reads). Default false.", optional = true)
    public boolean FLOW_END_LOCATION_SIGNIFICANT = false;

    @Advanced
    @Argument(fullName = FLOW_END_POS_UNCERTAINTY_LONG_NAME, doc = "Maximal number of bases of reads (fragment) ends difference that is marked as match (should only be applied to flow based reads). Default 0.", optional = true)
    public int ENDS_READ_UNCERTAINTY = 0;

    @Advanced
    @Argument(fullName = SINGLE_END_READS_CLIPPING_IS_END_LONG_NAME, doc = "Use clipped, rather than unclipped, when considering duplicates (should only be applied to flow based reads). Default false.", optional = true)
    public boolean FLOW_USE_CLIPPED_LOCATIONS = false;

    @Advanced
    @Argument(fullName = FLOW_SKIP_START_HOMOPOLYMERS_LONG_NAME, doc = "Skip first N flows, when considering duplicates (should only be applied to flow based reads). Default 0.", optional = true)
    public int FLOW_SKIP_START_HOMOPOLYMERS = 0;

    @Advanced
    @Argument(fullName = FLOW_Q_IS_KNOWN_END_LONG_NAME, doc = "Treat reads (fragment) clipped on tm:Q as known end position (should only be applied to flow based reads) (default: false)", optional = true)
    public boolean FLOW_Q_IS_KNOWN_END = false;

    @Advanced
    @Argument(fullName = FLOW_MD_MODE_LONG_NAME, optional = true, doc="Single argument for enabling the bulk of flow based features (should only be applied to flow based reads).")
    public Boolean useFlowFragments = false;

    public boolean isFlowEnabled() {
        return FLOW_QUALITY_SUM_STRATEGY || FLOW_END_LOCATION_SIGNIFICANT || FLOW_USE_CLIPPED_LOCATIONS || FLOW_SKIP_START_HOMOPOLYMERS != 0;
    }

    public String[] getFlowModeArgValues() {
        return new String[] {
                MarkDuplicatesSparkArgumentCollection.SINGLE_END_READS_END_POSITION_SIGNIFICANT, "true",
                MarkDuplicatesSparkArgumentCollection.SINGLE_END_READS_CLIPPING_IS_END_LONG_NAME, "true",
                MarkDuplicatesSparkArgumentCollection.FLOW_END_POS_UNCERTAINTY_LONG_NAME, "1",
                MarkDuplicatesSparkArgumentCollection.FLOW_SKIP_START_HOMOPOLYMERS_LONG_NAME, "0"
        };

    }
}