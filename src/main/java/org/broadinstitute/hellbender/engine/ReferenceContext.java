package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.reference.ReferenceSequence;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.iterators.ByteArrayIterator;

import java.util.Iterator;

/**
 * Wrapper around ReferenceDataSource that presents data from a specific interval/window to a client,
 * without improperly exposing the entire ReferenceDataSource interface.
 *
 * Reference data in the interval is lazily queried, so there's no overhead if the client chooses
 * not to examine contextual information from the reference.
 *
 * The reference interval can be optionally expanded by a configurable number of bases in each direction.
 * windowLeadingBases = 3 and windowTrailingBases = 5 means 3 bases of extra reference context before
 * the start of the interval and 5 bases of extra reference context after the end of the interval.
 *
 * Window boundaries can be set either at construction time or afterwards via {@link #setWindow}.
 *
 * A ReferenceContext may have no backing data source. In this case, queries on it will always
 * return empty arrays / iterators. You can determine whether there is a backing source of reference
 * data via {@link #hasBackingDataSource()}
 */
public final class ReferenceContext {

    /**
     * Backing data source. Null if there is no reference data.
     */
    private ReferenceDataSource dataSource;

    /**
     * Reference bases spanning this interval/window if a query has been performed. Null if we haven't been queried yet.
     * Cache is cleared if the window size changes between queries.
     */
    private ReferenceSequence cachedSequence;

    /**
     * Interval representing our location on the reference
     */
    private GenomeLoc interval;

    /**
     * Reference interval optionally expanded by configurable amount to produce the true query interval
     */
    private GenomeLoc window;


    /**
     * Create a ReferenceContext with no backing data source. This context will always return
     * empty arrays/iterators in response to queries.
     */
    public ReferenceContext() {
        this(null, null);
    }

    /**
     * Create a windowless ReferenceContext set up to lazily query the bases spanning just
     * the provided interval (with no extra bases of context)
     *
     * @param dataSource backing reference data source
     * @param interval interval to query, if we are accessed by a client
     */
    public ReferenceContext( final ReferenceDataSource dataSource, final GenomeLoc interval ) {
        this(dataSource, interval, 0, 0);
    }

    /**
     * Create a windowed ReferenceContext set up to lazily query the provided interval,
     * expanded by the specified number of bases in each direction.
     *
     * Window boundaries are cropped at contig boundaries, if necessary.
     *
     * @param dataSource backing reference data source (must be null if there is no reference)
     * @param interval our location on the reference (may be null if there is no reference)
     * @param windowLeadingBases Number of extra reference bases to include before the start of our interval. Must be >= 0.
     * @param windowTrailingBases Number of extra reference bases to include after the end of our interval. Must be >= 0.
     */
    public ReferenceContext( final ReferenceDataSource dataSource, final GenomeLoc interval, final int windowLeadingBases, final int windowTrailingBases ) {
        // If we have a backing data source, we must also have a query interval. If there's no data source,
        // we don't care about the interval (may be null or non-null).
        if ( dataSource != null && interval == null ) {
            throw new IllegalArgumentException("Must provide a non-null query interval for a ReferenceContext that has a backing ReferenceDataSource");
        }

        this.dataSource = dataSource;
        cachedSequence = null;
        this.interval = interval;
        setWindow(windowLeadingBases, windowTrailingBases);
    }

    /**
     * Determines whether this ReferenceContext has a backing reference data source. A ReferenceContext with
     * no backing data source will always return an empty bases array from {@link #getBases()} and an
     * empty iterator from {@link #getBasesIterator()}
     *
     * @return true if this ReferenceContext has a backing reference data source, otherwise false
     */
    public boolean hasBackingDataSource() {
        return dataSource != null;
    }

    /**
     * Get an iterator over the reference bases in this context
     *
     * Call {@link #setWindow} before calling this method if you want to configure the amount of extra reference context
     * to include around the current interval
     *
     * @return iterator over the reference bases in this context
     */
    public Iterator<Byte> getBasesIterator() {
        return dataSource != null ? dataSource.query(window) : new ByteArrayIterator(new byte[0]);
    }

    /**
     * Get all reference bases in this context. The results are cached in this object for future queries.
     * Will always return an empty array if there is no backing data source to query.
     *
     * Call {@link #setWindow} before calling this method if you want to configure the amount of extra reference context
     * to include around the current interval
     *
     * @return reference bases in this context, as a byte array
     */
    public byte[] getBases() {
        if ( dataSource == null ) {
            return new byte[0];
        }

        // Only perform a query if we haven't fetched the bases in this context previously
        if ( cachedSequence == null ) {
            cachedSequence = dataSource.queryAndPrefetch(window);
        }
        return cachedSequence.getBases();
    }

    /**
     * Get the location on the reference represented by this context, without including
     * any extra bases of requested context around this interval.
     *
     * @return location on the reference represented by this context as a GenomeLoc
     */
    public GenomeLoc getInterval() {
        return interval;
    }

    /**
     * Get the full expanded window of bases spanned by this context, including any extra
     * bases of requested context around the current interval.
     *
     * Note that the true window size may be smaller than originally requested due to cropping
     * at contig boundaries.
     *
     * @return full expanded window of bases spanned by this context as a GenomeLoc
     */
    public GenomeLoc getWindow() {
        return window;
    }

    /**
     * Set expanded window boundaries, subject to cropping at contig boundaries
     *
     * Allows the client to request a specific number of extra reference bases to include before
     * and after the bases within our interval. These extra bases will be returned by calls to
     * {@link #getBases} and {@link #getBasesIterator} in addition to the bases spanning our
     * actual interval.
     *
     * Note that the true window size may be smaller than requested due to cropping at contig boundaries.
     * Call {@link @numWindowLeadingBases} and {@link @numWindowTrailingBases} to get the actual
     * window dimensions.
     *
     * @param windowLeadingBases Number of extra reference bases to include before the start of our interval. Must be >= 0.
     * @param windowTrailingBases Number of extra reference bases to include after the end of our interval. Must be >= 0.
     */
    public void setWindow( final int windowLeadingBases, final int windowTrailingBases ) {
        if( windowLeadingBases < 0 ) throw new GATKException("Reference window starts after the current interval");
        if( windowTrailingBases < 0 ) throw new GATKException("Reference window ends before the current interval");

        if ( windowLeadingBases == 0 && windowTrailingBases == 0 ) {
            // the "windowless" case
            window = interval;
        }
        else {
            window = new GenomeLocParser(dataSource.getSequenceDictionary()).createGenomeLoc(interval.getContig(),
                                         calculateWindowStart(interval, windowLeadingBases), calculateWindowStop(interval, windowTrailingBases), true);
        }

        // Changing the window size invalidates our cached query result
        cachedSequence = null;
    }

    /**
     * Get the number of extra bases of context before the start of our interval, as configured
     * by a call to {@link #setWindow} or at construction time.
     *
     * Actual number of bases may be less than originally requested if the interval is near a contig boundary.
     *
     * @return number of extra bases of context before the start of our interval
     */
    public int numWindowLeadingBases() {
        return interval.getStart() - window.getStart();
    }

    /**
     * Get the number of extra bases of context after the end of our interval, as configured
     * by a call to {@link #setWindow} or at construction time.
     *
     * Actual number of bases may be less than originally requested if the interval is near a contig boundary.
     *
     * @return number of extra bases of context after the end of our interval
     */
    public int numWindowTrailingBases() {
        return window.getStop() - interval.getStop();
    }

    /**
     * Determines the start of the expanded reference window, bounded if necessary by the contig.
     *
     * @param locus The locus to expand.
     * @param windowLeadingBases number of bases to attempt to expand relative to the locus start (>= 0)
     * @return The start of the expanded window.
     */
    private int calculateWindowStart( final GenomeLoc locus, final int windowLeadingBases ) {
        return Math.max(locus.getStart() - windowLeadingBases, 1);
    }

    /**
     * Determines the stop of the expanded reference window, bounded if necessary by the contig.
     *
     * @param locus The locus to expand.
     * @param windowTrailingBases number of bases to attempt to expand relative to the locus end (>= 0)
     * @return The end of the expanded window.
     */
    private int calculateWindowStop( final GenomeLoc locus, final int windowTrailingBases ) {
        final int sequenceLength = dataSource.getSequenceDictionary().getSequence(locus.getContigIndex()).getSequenceLength();
        return Math.min(locus.getStop() + windowTrailingBases, sequenceLength);
    }
}
