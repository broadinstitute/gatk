
package org.broadinstitute.hellbender.cmdline.argumentcollections;

import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollectionDefinition;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.*;

import java.util.ArrayList;
import java.util.List;

public final class IntervalArgumentCollection implements ArgumentCollectionDefinition {
    /**
     * Use this argument to perform the analysis over only part of the genome. This argument can be specified multiple times.
     * You can use samtools-style intervals either explicitly on the command line (e.g. -L 1 or -L 1:100-200) or
     * by loading in a file containing a list of intervals (e.g. -L myFile.intervals).
     */
    @Argument(fullName = "intervals", shortName = "L", doc = "One or more genomic intervals over which to operate", optional = true)
    final protected List<String> intervalStrings = new ArrayList<>();

    /**
     * Use this argument to exclude certain parts of the genome from the analysis (like -L, but the opposite).
     * This argument can be specified multiple times. You can use samtools-style intervals either explicitly on the
     * command line (e.g. -XL 1 or -XL 1:100-200) or by loading in a file containing a list of intervals
     * (e.g. -XL myFile.intervals).
     *
     * */
    @Argument(fullName = "excludeIntervals", shortName = "XL", doc = "One or more genomic intervals to exclude from processing", optional = true)
    final protected List<String> excludeIntervalStrings = new ArrayList<>();

    /**
     * By default, the program will take the UNION of all intervals specified using -L and/or -XL. However, you can
     * change this setting for -L, for example if you want to take the INTERSECTION of the sets instead. E.g. to
     * perform the analysis only on chromosome 1 exomes, you could specify -L exomes.intervals -L 1 --interval_set_rule
     * INTERSECTION. However, it is not possible to modify the merging approach for intervals passed using -XL (they will
     * always be merged using UNION).
     *
     * Note that if you specify both -L and -XL, the -XL interval set will be subtracted from the -L interval set.
     */
    @Argument(fullName = "interval_set_rule", shortName = "isr", doc = "Set merging approach to use for combining interval inputs")
    protected IntervalSetRule intervalSetRule = IntervalSetRule.UNION;

    /**
     * Use this to add padding to the intervals specified using -L and/or -XL. For example, '-L 1:100' with a
     * padding value of 20 would turn into '-L 1:80-120'. This is typically used to add padding around exons when
     * analyzing exomes.
     */
    @Argument(fullName = "interval_padding", shortName = "ip", doc = "Amount of padding (in bp) to add to each interval")
    protected int intervalPadding = 0;

    final protected IntervalMergingRule intervalMerging = IntervalMergingRule.ALL;
    static private Logger logger = LogManager.getLogger(IntervalArgumentCollection.class);

    /**
     * Get the intervals specified on the command line.
     * @param sequenceDict used to validate intervals
     * @return a list of the given intervals after processing and validation
     */
    public List<SimpleInterval> getIntervals(SAMSequenceDictionary sequenceDict){
        return getIntervals(new GenomeLocParser(sequenceDict));
    }

    /**
     * Get the intervals specified on the command line.
     * @param genomeLocParser used to validate the intervals
     * @return list of the given intervals after processing and validation
     */
    public List<SimpleInterval> getIntervals(final GenomeLocParser genomeLocParser) {
        // return if no interval arguments at all
        if (!intervalsSpecified()) {
            throw new GATKException("Cannot call getIntervals() without specifying either intervals to include or exclude.");
        }


        GenomeLocSortedSet includeSortedSet;
        if (intervalStrings.isEmpty()){
            // if -L argument isn't given but -XL is, set include to all of the territory
            includeSortedSet = GenomeLocSortedSet.createSetFromSequenceDictionary(genomeLocParser.getSequenceDictionary());
        } else {
            try {
                includeSortedSet = IntervalUtils.loadIntervals(intervalStrings, intervalSetRule, intervalMerging, intervalPadding, genomeLocParser);
            } catch( UserException.EmptyIntersection e) {
                throw new UserException.BadArgumentValue("-L, --interval_set_rule", intervalStrings+","+intervalSetRule, "The specified intervals had an empty intersection");
            }
        }
        final GenomeLocSortedSet excludeSortedSet = IntervalUtils.loadIntervals(excludeIntervalStrings, IntervalSetRule.UNION, intervalMerging, 0, genomeLocParser);


        GenomeLocSortedSet intervals;
        // if no exclude arguments, can return the included set directly
        if ( excludeSortedSet.isEmpty() ) {
            intervals = includeSortedSet;
        }// otherwise there are exclude arguments => must merge include and exclude GenomeLocSortedSets
        else {
            intervals = includeSortedSet.subtractRegions(excludeSortedSet);

            if( intervals.isEmpty()){
                throw new UserException.BadArgumentValue("-L,-XL",intervalStrings.toString() + ", "+excludeIntervalStrings.toString(),"The intervals specified for exclusion with -XL removed all territory specified by -L.");
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
        return convertGenomeLocsToSimpleIntervals(intervals.toList());
    }

    /**
     * Convert a List of intervals in GenomeLoc format into a List of intervals in SimpleInterval format.
     *
     * @param genomeLocIntervals list of GenomeLoc intervals to convert
     * @return equivalent List of SimpleIntervals
     */
    private List<SimpleInterval> convertGenomeLocsToSimpleIntervals( final List<GenomeLoc> genomeLocIntervals ) {
        List<SimpleInterval> convertedIntervals = new ArrayList<>(genomeLocIntervals.size());
        for ( GenomeLoc genomeLoc : genomeLocIntervals ) {
            if ( genomeLoc.isUnmapped() ) {
                throw new UserException("Unmapped intervals are not currently supported");
            }

            convertedIntervals.add(new SimpleInterval(genomeLoc.getContig(), genomeLoc.getStart(), genomeLoc.getStop()));
        }
        return convertedIntervals;
    }

    /**
     * Have any intervals been specified for inclusion or exclusion
     */
    public boolean intervalsSpecified() {
        return !( intervalStrings.isEmpty() && excludeIntervalStrings.isEmpty());
    }
}

