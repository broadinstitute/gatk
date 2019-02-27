package org.broadinstitute.hellbender.cmdline.argumentcollections;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.engine.TraversalParameters;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.*;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

/**
 * Intended to be used as an @ArgumentCollection for specifying intervals at the command line.
 * Subclasses must override getIntervalStrings and addToIntervalStrings().
 */
public abstract class IntervalArgumentCollection implements Serializable {
    private static final Logger logger = LogManager.getLogger(IntervalArgumentCollection.class);
    private static final long serialVersionUID = 1L;

    public static final String EXCLUDE_INTERVALS_LONG_NAME = "exclude-intervals";
    public static final String EXCLUDE_INTERVALS_SHORT_NAME = "XL";
    public static final String INTERVAL_SET_RULE_LONG_NAME = "interval-set-rule";
    public static final String INTERVAL_PADDING_LONG_NAME = "interval-padding";
    public static final String INTERVAL_EXCLUSION_PADDING_LONG_NAME = "interval-exclusion-padding";
    public static final String INTERVAL_MERGING_RULE_LONG_NAME = "interval-merging-rule";

    /**
     * Subclasses must provide a -L argument and override this to return the results of that argument.
     *
     * The -L argument specifies intervals to include in analysis and has the following semantics
     * It can use samtools-style intervals either explicitly on the command line (e.g. -L 1 or -L 1:100-200) or
     * by loading in a file containing a list of intervals (e.g. -L myFile.intervals).
     * It can be specified multiple times.
     *
     * @return string gathered from the command line -L argument to be parsed into intervals to include
     */
    protected abstract List<String> getIntervalStrings();

    /**
     * Add an extra interval string to the intervals to include.
     * ONLY for testing -- will throw if called after interval parsing has been performed.
     */
    @VisibleForTesting
    protected abstract void addToIntervalStrings(String newInterval);

    /**
     * Use this argument to exclude certain parts of the genome from the analysis (like -L, but the opposite).
     * This argument can be specified multiple times. You can use samtools-style intervals either explicitly on the
     * command line (e.g. -XL 1 or -XL 1:100-200) or by loading in a file containing a list of intervals
     * (e.g. -XL myFile.intervals).
     * @return strings gathered from the command line -XL argument to be parsed into intervals to exclude
     */
    // Interval list files such as Picard interval lists are structured and require specialized parsing that
    // is handled by IntervalUtils, so use suppressFileExpansion to bypass command line parser auto-expansion.
    @Argument(fullName = EXCLUDE_INTERVALS_LONG_NAME, shortName = EXCLUDE_INTERVALS_SHORT_NAME, doc = "One or more genomic intervals to exclude from processing",
            suppressFileExpansion = true, optional = true, common = true)
    protected final List<String> excludeIntervalStrings = new ArrayList<>();

    /**
     * By default, the program will take the UNION of all intervals specified using -L and/or -XL. However, you can
     * change this setting for -L, for example if you want to take the INTERSECTION of the sets instead. E.g. to
     * perform the analysis only on chromosome 1 exomes, you could specify -L exomes.intervals -L 1 --interval-set-rule
     * INTERSECTION. However, it is not possible to modify the merging approach for intervals passed using -XL (they will
     * always be merged using UNION).
     *
     * Note that if you specify both -L and -XL, the -XL interval set will be subtracted from the -L interval set.
     */
    @Argument(fullName = INTERVAL_SET_RULE_LONG_NAME, shortName = "isr", doc = "Set merging approach to use for combining interval inputs", common = true)
    protected IntervalSetRule intervalSetRule = IntervalSetRule.UNION;
    /**
     * Use this to add padding to the intervals specified using -L. For example, '-L 1:100' with a
     * padding value of 20 would turn into '-L 1:80-120'. This is typically used to add padding around targets when
     * analyzing exomes.
     */
    @Argument(fullName = INTERVAL_PADDING_LONG_NAME, shortName = "ip", doc = "Amount of padding (in bp) to add to each interval you are including.", common = true)
    protected int intervalPadding = 0;
    /**
     * Use this to add padding to the intervals specified using -XL. For example, '-XL 1:100' with a
     * padding value of 20 would turn into '-XL 1:80-120'. This is typically used to add padding around targets when
     * analyzing exomes.
     */
    @Argument(fullName = INTERVAL_EXCLUSION_PADDING_LONG_NAME, shortName= "ixp", doc = "Amount of padding (in bp) to add to each interval you are excluding.", common = true)
    protected int intervalExclusionPadding = 0;

    /**
     * By default, the program merges abutting intervals (i.e. intervals that are directly side-by-side but do not
     * actually overlap) into a single continuous interval. However you can change this behavior if you want them to be
     * treated as separate intervals instead.
     */
    @Argument(fullName = INTERVAL_MERGING_RULE_LONG_NAME, shortName = "imr", doc = "Interval merging rule for abutting intervals", optional = true)
    protected IntervalMergingRule intervalMergingRule = IntervalMergingRule.ALL;

    /**
     * Full parameters for traversal, including our parsed intervals and a flag indicating whether unmapped records
     * should be returned. Lazily initialized.
     */
    protected TraversalParameters traversalParameters = null;

    /**
     * Get the intervals specified on the command line.
     * @param sequenceDict used to validate intervals
     * @return a list of the given intervals after processing and validation
     */
    public List<SimpleInterval> getIntervals( final SAMSequenceDictionary sequenceDict ){
        return getTraversalParameters(sequenceDict).getIntervalsForTraversal();
    }

    /**
     * Get the largest interval per contig that contains the intervals specified on the command line.
     * @param sequenceDict used to validate intervals
     * @return a list of one interval per contig spanning the input intervals after processing and validation
     */
    public List<SimpleInterval> getSpanningIntervals( final SAMSequenceDictionary sequenceDict ){
        return IntervalUtils.getSpanningIntervals(getIntervals(sequenceDict), sequenceDict);
    }

    /**
     * Get the interval set rule specified on the command line.
     */
    public IntervalSetRule getIntervalSetRule() {
        return intervalSetRule;
    }

    /**
     * Get the interval padding specified on the command line.
     */
    public int getIntervalPadding() {
        return intervalPadding;
    }

    /**
     * Get the interval exclusion padding specified on the command line.
     */
    public int getIntervalExclusionPadding() {
        return intervalExclusionPadding;
    }

    /**
     * Get the interval merging rule specified on the command line.
     */
    public IntervalMergingRule getIntervalMergingRule() {
        return intervalMergingRule;
    }

    /**
     * Returns the full set of traversal parameters specified on the command line, including the parsed intervals
     * and a flag indicating whether unmapped records were requested.
     *
     * @param sequenceDict used to validate intervals
     * @return the full set of traversal parameters specified on the command line
     */
    public TraversalParameters getTraversalParameters( final SAMSequenceDictionary sequenceDict ) {
        if ( ! intervalsSpecified() ) {
            throw new GATKException("Cannot call getTraversalParameters() without specifying either intervals to include or exclude.");
        }

        if ( traversalParameters == null ) {
            parseIntervals(new GenomeLocParser(sequenceDict));
        }

        return traversalParameters;
    }

    private void parseIntervals(final GenomeLocParser genomeLocParser) {
        // return if no interval arguments at all
        if (!intervalsSpecified()) {
            throw new GATKException("Cannot call parseIntervals() without specifying either intervals to include or exclude.");
        }

        GenomeLocSortedSet includeSortedSet;
        if (getIntervalStrings().isEmpty()){
            // the -L argument isn't specified, which means that -XL was, since we checked intervalsSpecified()
            // therefore we set the include set to be the entire reference territory
            includeSortedSet = GenomeLocSortedSet.createSetFromSequenceDictionary(genomeLocParser.getSequenceDictionary());
        } else {
            try {
                includeSortedSet = IntervalUtils.loadIntervals(getIntervalStrings(), intervalSetRule, intervalMergingRule, intervalPadding, genomeLocParser);
            } catch( UserException.EmptyIntersection e) {
                throw new CommandLineException.BadArgumentValue("-L, --" + IntervalArgumentCollection.INTERVAL_SET_RULE_LONG_NAME, getIntervalStrings()+","+intervalSetRule, "The specified intervals had an empty intersection");
            }
        }

        final GenomeLocSortedSet excludeSortedSet = IntervalUtils.loadIntervals(excludeIntervalStrings, IntervalSetRule.UNION, intervalMergingRule, intervalExclusionPadding, genomeLocParser);
        if ( excludeSortedSet.contains(GenomeLoc.UNMAPPED) ) {
            throw new UserException("-XL unmapped is not currently supported");
        }

        GenomeLocSortedSet intervals;
        // if no exclude arguments, can return the included set directly
        if ( excludeSortedSet.isEmpty() ) {
            intervals = includeSortedSet;
        }// otherwise there are exclude arguments => must merge include and exclude GenomeLocSortedSets
        else {
            intervals = includeSortedSet.subtractRegions(excludeSortedSet);

            if( intervals.isEmpty()){
                throw new CommandLineException.BadArgumentValue("-L,-XL",getIntervalStrings().toString() + ", "+excludeIntervalStrings.toString(),"The intervals specified for exclusion with -XL removed all territory specified by -L.");
            }
            // logging messages only printed when exclude (-XL) arguments are given
            final long toPruneSize = includeSortedSet.coveredSize();
            final long toExcludeSize = excludeSortedSet.coveredSize();
            final long intervalSize = intervals.coveredSize();
            logger.info(String.format("Initial include intervals span %d loci; exclude intervals span %d loci", toPruneSize, toExcludeSize));
            logger.info(String.format("Excluding %d loci from original intervals (%.2f%% reduction)",
                    toPruneSize - intervalSize, (toPruneSize - intervalSize) / (0.01 * toPruneSize)));
        }

        logger.info(String.format("Processing %d bp from intervals", intervals.coveredSize()));

        // Separate out requests for unmapped records from the rest of the intervals.
        boolean traverseUnmapped = false;
        if ( intervals.contains(GenomeLoc.UNMAPPED) ) {
            traverseUnmapped = true;
            intervals.remove(GenomeLoc.UNMAPPED);
        }

        traversalParameters = new TraversalParameters(IntervalUtils.convertGenomeLocsToSimpleIntervals(intervals.toList()), traverseUnmapped);
    }


    /**
     * Have any intervals been specified for inclusion or exclusion
     */
    public boolean intervalsSpecified() {
        return !( getIntervalStrings().isEmpty() && excludeIntervalStrings.isEmpty());
    }
}
