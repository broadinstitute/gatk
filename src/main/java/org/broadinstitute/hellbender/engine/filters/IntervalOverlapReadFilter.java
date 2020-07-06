package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.util.OverlapDetector;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.ReadFilterArgumentDefinitions;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.GenomeLocSortedSet;
import org.broadinstitute.hellbender.utils.IntervalMergingRule;
import org.broadinstitute.hellbender.utils.IntervalSetRule;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.logging.OneShotLogger;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.lang.annotation.Documented;
import java.util.ArrayList;
import java.util.List;

/**
 * A simple read filter that allows for the user to specify intervals at the filtering stage.
 *
 * NOTE: This class is intended to be a convenience method for the very specific case where a user might want to subset
 *       a bam file by intervals but they cannot for whatever reason index/resort the file to be sorted by position. As
 *       a consequence a file subsetted using this filter will involve reading over the entire bam input and will consequently
 *       be very slow. The preferred method for subsetting a bam file in this case is to use the -L command to subset using
 *       an index to avoid reading parts of the bam outside of the specified subset.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
@DocumentedFeature(groupName= HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Filters out reads that don't overlap the specified region. NOTE: This approach to extracting overlapping reads is very slow compared to using PrintReads and -L on an indexed bam file.")
public final class IntervalOverlapReadFilter extends ReadFilter {
    private static final long serialVersionUID = 1L;

    // Interval list files such as Picard interval lists are structured and require specialized parsing that
    // is handled by IntervalUtils, so use suppressFileExpansion to bypass command line parser auto-expansion.
    @Argument(fullName = ReadFilterArgumentDefinitions.KEEP_INTERVAL_NAME, suppressFileExpansion = true, doc = "One or more genomic intervals to keep", optional = false)
    protected final List<String> intervalStrings;

    private OverlapDetector<?> detector;
    
    protected final OneShotLogger warning = new OneShotLogger(this.getClass());

    /**
     * Default constructor.
     */
    public IntervalOverlapReadFilter() {
        this.intervalStrings = new ArrayList<>();
    }

    /**
     * Constructor with provided intervals.
     *
     * @param intervalStrings intervals to use.
     */
    public IntervalOverlapReadFilter(final List<String> intervalStrings) {
        this.intervalStrings = new ArrayList<>(intervalStrings);
    }

    @Override
    public boolean test(final GATKRead read) {
        return getDetector().overlapsAny(read);
    }

    private synchronized OverlapDetector<?>  getDetector() {
        if (detector == null) {
            warning.warn("You are using the IntervalOverlapReadFilter to subset your input, this is a very slow operation and it will usually be preferable to use the '-L' command or provide an interval list file for a sorted and indexed input file");
            // first load the intervals
            final GenomeLocParser genomeLocParser = new GenomeLocParser(this.samHeader.getSequenceDictionary());
            final GenomeLocSortedSet intervals = IntervalUtils.loadIntervals(intervalStrings, IntervalSetRule.UNION, IntervalMergingRule.ALL, 0, genomeLocParser);
            // then initializes the overlap detector with the loaded intervals
            this.detector = OverlapDetector.create(intervals.toList());
        }
        return this.detector;
    }
}
