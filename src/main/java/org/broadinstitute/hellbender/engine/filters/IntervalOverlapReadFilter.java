package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.util.OverlapDetector;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.ReadFilterArgumentDefinitions;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.GenomeLocSortedSet;
import org.broadinstitute.hellbender.utils.IntervalMergingRule;
import org.broadinstitute.hellbender.utils.IntervalSetRule;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public final class IntervalOverlapReadFilter extends ReadFilter {
    private static final long serialVersionUID = 1L;

    // Interval list files such as Picard interval lists are structured and require specialized parsing that
    // is handled by IntervalUtils, so use suppressFileExpansion to bypass command line parser auto-expansion.
    @Argument(fullName = ReadFilterArgumentDefinitions.KEEP_INTERVAL_NAME, suppressFileExpansion = true, doc = "One or more genomic intervals to keep", optional = false)
    protected final List<String> intervalStrings;

    private OverlapDetector<?> detector;

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
        createDetector();
        return detector.overlapsAny(read);
    }

    private synchronized void createDetector() {
        if (detector == null) {
            // first load the intervals
            final GenomeLocParser genomeLocParser = new GenomeLocParser(this.samHeader.getSequenceDictionary());
            final GenomeLocSortedSet intervals = IntervalUtils.loadIntervals(intervalStrings, IntervalSetRule.UNION, IntervalMergingRule.ALL, 0, genomeLocParser);
            // then initializes the overlap detector with the loaded intervals
            this.detector = OverlapDetector.create(intervals.toList());
        }
    }
}
