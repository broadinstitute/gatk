package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.reference.ReferenceSequence;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.GenomeLocParser;

import java.util.Iterator;

/**
 * Wrapper around ReferenceDataSource that presents data from a specific interval/window to a client,
 * without improperly exposing the entire ReferenceDataSource interface.
 *
 * Reference data in the interval is lazily queried, so there's no overhead if the client chooses
 * not to examine contextual information from the reference.
 *
 * The reference interval can be optionally expanded by a configurable number of bases in each direction.
 * A windowStartOffset and windowStopOffset of [-3, 5] means 3 bases of extra context before
 * the start of the interval and 5 bases of extra context after the end of the interval.
 */
public class ReferenceContext {

    /**
     * Backing data source
     */
    private ReferenceDataSource dataSource;

    /**
     * Reference bases spanning this interval/window if a query has been performed. Null if we haven't been queried yet.
     */
    private ReferenceSequence cachedSequence;

    /**
     * Interval representing our location on the reference
     */
    private GenomeLoc interval;

    /**
     * Reference interval optionally expanded by requested window to produce the true query interval
     */
    private GenomeLoc window;

    /**
     * Offset from the interval start representing extra bases to include in the query. Must be negative.
     */
    private int windowStartOffset;

    /**
     * Offset from the interval stop representing extra bases to include in the query. Must be positive.
     */
    private int windowStopOffset;

    /**
     * Create a windowless ReferenceContext set up to lazily query the provided interval
     *
     * @param dataSource backing reference data source
     * @param interval interval to query, if we are accessed by a client
     */
    public ReferenceContext( final ReferenceDataSource dataSource, final GenomeLoc interval ) {
        this(dataSource, interval, 0, 0);
    }

    /**
     * Create a windowed ReferenceContext set up to lazily query the provided interval,
     * expanded according to the window start/stop offsets.
     *
     * Window boundaries are cropped at contig boundaries, if necessary.
     *
     * @param dataSource backing reference data source
     * @param interval our location on the reference
     * @param windowStartOffset Offset from the interval start representing extra bases to include in the context. Must be negative.
     * @param windowStopOffset Offset from the interval stop representing extra bases to include in the context. Must be positive.
     */
    public ReferenceContext( final ReferenceDataSource dataSource, final GenomeLoc interval, final int windowStartOffset, final int windowStopOffset ) {
        if ( dataSource == null || interval == null ) {
            throw new IllegalArgumentException("dataSource/interval must be non-null");
        }

        this.dataSource = dataSource;
        cachedSequence = null;
        this.interval = interval;
        setWindow(windowStartOffset, windowStopOffset);
    }

    /**
     * Get an iterator over the reference bases in this context
     *
     * @return iterator over the reference bases in this context
     */
    public Iterator<Byte> getBasesIterator() {
        return dataSource.query(window);
    }

    /**
     * Get all reference bases in this context
     *
     * @return reference bases in this context, as a byte array
     */
    public byte[] getBases() {
        // Only perform a query if we haven't fetched the bases in this context previously
        if ( cachedSequence == null ) {
            cachedSequence = dataSource.queryAndPrefetch(window);
        }
        return cachedSequence.getBases();
    }

    /**
     * Get the location on the reference represented by this context
     *
     * @return location on the reference represented by this context as a GenomeLoc
     */
    public GenomeLoc getInterval() {
        return interval;
    }

    /**
     * Get the full expanded window of bases spanned by this context
     *
     * @return full expanded window of bases spanned by this context as a GenomeLoc
     */
    public GenomeLoc getWindow() {
        return window;
    }

    /**
     * Get the number of extra bases of context before the start of our interval
     *
     * @return number of extra bases of context before the start of our interval
     */
    public int getLeadingWindowSize() {
        return interval.getStart() - window.getStart();
    }

    /**
     * Get the number of extra bases of context after the end of our interval
     *
     * @return number of extra bases of context after the end of our interval
     */
    public int getTrailingWindowSize() {
        return window.getStop() - interval.getStop();
    }

    /**
     * Set expanded window boundaries, subject to validation
     *
     * @param windowStartOffset Offset from interval start at which to start the window. Must be negative.
     * @param windowStopOffset Offset from interval stop at which to stop the window. Must be positive.
     */
    private void setWindow( final int windowStartOffset, final int windowStopOffset ) {
        if( windowStartOffset > 0 ) throw new GATKException("Reference window starts after the current locus");
        if( windowStopOffset < 0 ) throw new GATKException("Reference window ends before the current locus");

        this.windowStartOffset = windowStartOffset;
        this.windowStopOffset = windowStopOffset;

        if ( windowStartOffset == 0 && windowStopOffset == 0 ) {
            // the "windowless" case
            window = interval;
        }
        else {
            window = new GenomeLocParser(dataSource.getSequenceDictionary()).createGenomeLoc(interval.getContig(), getWindowStart(interval), getWindowStop(interval), true);
        }
    }

    /**
     * Gets the start of the expanded window, bounded if necessary by the contig.
     *
     * @param locus The locus to expand.
     * @return The start of the expanded window.
     */
    private int getWindowStart( final GenomeLoc locus ) {
        return Math.max(locus.getStart() + windowStartOffset, 1);
    }

    /**
     * Gets the stop of the expanded window, bounded if necessary by the contig.
     *
     * @param locus The locus to expand.
     * @return The end of the expanded window.
     */
    private int getWindowStop( final GenomeLoc locus ) {
        final int sequenceLength = dataSource.getSequenceDictionary().getSequence(locus.getContigIndex()).getSequenceLength();
        return Math.min(locus.getStop() + windowStopOffset, sequenceLength);
    }
}
