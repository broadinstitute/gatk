package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.iterators.ByteArrayIterator;

import java.util.Arrays;
import java.util.Iterator;

/**
 * Wrapper around ReferenceDataSource that presents data from a specific interval/window to a client,
 * without improperly exposing the entire ReferenceDataSource interface.
 *
 * Reference data in the interval is lazily queried, so there's no overhead if the client chooses
 * not to examine contextual information from the reference.
 *
 * Reference bases are returned as bytes for integration with the existing {@link org.broadinstitute.hellbender.utils.BaseUtils}
 * implementation.
 *
 * The reference interval can be optionally expanded by a configurable number of bases in each direction.
 * windowLeadingBases = 3 and windowTrailingBases = 5 means 3 bases of extra reference context before
 * the start of the interval and 5 bases of extra reference context after the end of the interval.
 *
 * Window boundaries can be set either at construction time or afterwards via {@link #setWindow}.
 *
 * A ReferenceContext may have no backing data source and/or interval. In these cases, queries on it will always
 * return empty arrays / iterators. You can determine whether there is a backing source of reference
 * data via {@link #hasBackingDataSource()}, and whether there is an interval via {@link #getInterval}.
 */
public final class ReferenceContext implements Iterable<Byte> {

    /**
     * Backing data source. Null if there is no reference data.
     */
    private final ReferenceDataSource dataSource;

    /**
     * Interval representing our location on the reference. May be null if, eg., we're dealing with unmapped data.
     */
    private final SimpleInterval interval;

    /**
     * Reference interval optionally expanded by a configurable amount to produce the true query interval.
     * Will be null if this context lacks an interval.
     */
    private SimpleInterval window;

    /**
     * Reference bases spanning this interval/window if a query has been performed. Null if we haven't been queried yet.
     * Cache is cleared if the window size changes between queries.
     */
    private ReferenceSequence cachedSequence;


    /**
     * Create a ReferenceContext with no backing data source. This context will always return
     * empty arrays/iterators in response to queries.
     */
    public ReferenceContext() {
        this(null, null, 0, 0);
    }

    /**
     * Create a windowless ReferenceContext set up to lazily query the bases spanning just
     * the provided interval (with no extra bases of context)
     *
     * @param dataSource backing reference data source (may be null if there is no reference)
     * @param interval interval to query, if we are accessed by a client (may be null if our location is unknown)
     */
    public ReferenceContext( final ReferenceDataSource dataSource, final SimpleInterval interval ) {
        this(dataSource, interval, 0, 0);
    }

    /**
     * Create a windowed ReferenceContext set up to lazily query the provided interval,
     * expanded by the specified number of bases in each direction.
     *
     * Window boundaries are cropped at contig boundaries, if necessary.
     *
     * @param dataSource backing reference data source (may be null if there is no reference)
     * @param interval our location on the reference (may be null if our location is unknown)
     * @param windowLeadingBases Number of extra reference bases to include before the start of our interval. Must be >= 0.
     * @param windowTrailingBases Number of extra reference bases to include after the end of our interval. Must be >= 0.
     */
    public ReferenceContext( final ReferenceDataSource dataSource, final SimpleInterval interval, final int windowLeadingBases, final int windowTrailingBases ) {
        this.dataSource = dataSource;
        this.cachedSequence = null;
        this.interval = interval;
        setWindow(windowLeadingBases, windowTrailingBases);
    }

    /**
     * Create a windowed ReferenceContext set up to lazily query the provided interval.
     *
     * Window is preserved from {@code thatReferenceContext}.
     *
     * @param thatContext An existing {@link ReferenceContext} on which to base this new one.
     * @param interval our location on the reference (may be null if our location is unknown)
     */
    public ReferenceContext( final ReferenceContext thatContext, final SimpleInterval interval ) {
        this.dataSource = thatContext.dataSource;
        this.cachedSequence = null;
        this.interval = interval;

        // Determine the window:
        final int windowLeadingBases = thatContext.numWindowLeadingBases();
        final int windowTrailingBases = thatContext.numWindowTrailingBases();

        setWindow(windowLeadingBases, windowTrailingBases);
    }

    /**
     * Create a windowed ReferenceContext set up to lazily query the provided interval,
     * expanded by the specified number of bases in each direction.
     *
     * Window boundaries are cropped at contig boundaries, if necessary.
     *
     * @param dataSource backing reference data source (may be null if there is no reference)
     * @param interval our location on the reference (may be null if our location is unknown)
     * @param window the expanded location on the reference. May be null if our location is unknown or there is no expanded window
     *               (ie., the interval == window case). Must be null if interval is null. Must contain interval
     *               if both are non-null.
     */
    public ReferenceContext( final ReferenceDataSource dataSource, final SimpleInterval interval, final SimpleInterval window ) {
        this.dataSource = dataSource;
        this.cachedSequence = null;
        this.interval = interval;
        Utils.validateArg(interval != null || window == null, () -> "if interval is null then window must be null too but was " + window);
        Utils.validateArg( interval == null || window == null || window.contains(interval), () ->
                "window " + window + " does not contain the interval " + interval);

        // The "windowless" case
        if ( window == null ) {
            this.window = interval;
        } else {
            this.window = new SimpleInterval(interval.getContig(),
                    trimToContigStart(window.getStart()),
                    trimToContigLength(interval.getContig(), window.getEnd()));
        }
    }

    /**
     * Determines whether this ReferenceContext has a backing reference data source. A ReferenceContext with
     * no backing data source will always return an empty bases array from {@link #getBases()} and an
     * empty iterator from {@link #iterator()}
     *
     * @return true if this ReferenceContext has a backing reference data source, otherwise false
     */
    public boolean hasBackingDataSource() {
        return dataSource != null;
    }

    /**
     * Get an iterator over the reference bases in this context. Will return an empty iterator if this
     * context has no backing data source and/or interval.
     *
     * Call {@link #setWindow} before calling this method if you want to configure the amount of extra reference context
     * to include around the current interval
     *
     * @return iterator over the reference bases in this context
     */
    @Override
    public Iterator<Byte> iterator() {
        return dataSource != null && window != null ? dataSource.query(window) : new ByteArrayIterator(new byte[0]);
    }

    /**
     * Get all reference bases in this context. The results are cached in this object for future queries.
     * Will always return an empty array if there is no backing data source and/or interval to query.
     *
     * Call {@link #setWindow} before calling this method if you want to configure the amount of extra reference context
     * to include around the current interval
     *
     * @return reference bases in this context, as a byte array
     */
    public byte[] getBases() {
        if ( dataSource == null || window == null ) {
            return new byte[0];
        }

        // Only perform a query if we haven't fetched the bases in this context previously
        if ( cachedSequence == null ) {
            cachedSequence = dataSource.queryAndPrefetch(window);
        }
        return cachedSequence.getBases();
    }

    /**
     * Get all reference bases in this context with the given window.
     * Does not cache results or modify this {@link ReferenceContext} at all.
     * Will always return an empty array if there is no backing data source and/or interval to query.
     *
     * @return reference bases in this context, as a byte array
     */
    public byte[] getBases(final SimpleInterval window) {
        if ( dataSource == null || window == null ) {
            return new byte[0];
        }

        // Trim to the contig start/end:
        final SimpleInterval trimmedWindow = new SimpleInterval(
                window.getContig(),
                trimToContigStart(window.getStart()),
                trimToContigLength(window.getContig(), window.getEnd())
        );

        return dataSource.queryAndPrefetch(trimmedWindow).getBases();
    }

    /**
     * Get all reference bases in this context with the given leading / trailing bases as the window.
     * Uses the current {@link ReferenceContext#window} as a basis for the position.
     * Does not cache results or modify this {@link ReferenceContext} at all.
     * Will always return an empty array if there is no backing data source and/or interval to query.
     *
     * @return reference bases in this context, as a byte array
     */
    public byte[] getBases(final int windowLeadingBases, final int windowTrailingBases) {
        if ( dataSource == null || window == null ) {
            return new byte[0];
        }

        // Trim to the contig start/end:
        final SimpleInterval trimmedWindow = new SimpleInterval(
                window.getContig(),
                trimToContigStart(window.getStart() - windowLeadingBases),
                trimToContigLength(window.getContig(), window.getEnd() + windowTrailingBases)
        );

        return dataSource.queryAndPrefetch(trimmedWindow).getBases();
    }

    /**
     * Get the bases in this context, from the beginning of the interval to the end of the window.
     */
    public byte[] getForwardBases() {
        final byte[] bases = getBases();
        final int mid = interval.getStart() - window.getStart();
        return new String(bases).substring(mid).getBytes();
    }


    /**
     * Get the location on the reference represented by this context, without including
     * any extra bases of requested context around this interval.
     *
     * @return location on the reference represented by this context as a SimpleInterval
     *         (may be null if we have no known location)
     */
    public SimpleInterval getInterval() {
        return interval;
    }

    /**
     * Get the full expanded window of bases spanned by this context, including any extra
     * bases of requested context around the current interval.
     *
     * Note that the true window size may be smaller than originally requested due to cropping
     * at contig boundaries.
     *
     * @return full expanded window of bases spanned by this context as a SimpleInterval
     *         (will be null if this context has no interval)
     */
    public SimpleInterval getWindow() {
        return window;
    }

    /**
     * Set expanded window boundaries, subject to cropping at contig boundaries
     *
     * Allows the client to request a specific number of extra reference bases to include before
     * and after the bases within our interval. These extra bases will be returned by calls to
     * {@link #getBases} and {@link #iterator} in addition to the bases spanning our
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
        if( windowLeadingBases < 0 ) {
            throw new GATKException("Reference window starts after the current interval");
        }
        if( windowTrailingBases < 0 ) {
            throw new GATKException("Reference window ends before the current interval");
        }

        if ( interval == null || (windowLeadingBases == 0 && windowTrailingBases == 0) ) {
            // the "windowless" case
            window = interval;
        }
        else {
            window = new SimpleInterval(interval.getContig(),
                    calculateWindowStart(interval, windowLeadingBases),
                    calculateWindowStop(interval, windowTrailingBases));
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
        return window == null ? 0 : interval.getStart() - window.getStart();
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
        return window == null ? 0 : window.getEnd() - interval.getEnd();
    }

    /**
     * Determines the start of the expanded reference window, bounded by 1.
     *
     * @param locus The locus to expand.
     * @param windowLeadingBases number of bases to attempt to expand relative to the locus start (>= 0)
     * @return The start of the expanded window.
     */
    private int calculateWindowStart( final SimpleInterval locus, final int windowLeadingBases ) {
        return trimToContigStart(locus.getStart() - windowLeadingBases);
    }

    /**
     * Determines the start of the expanded reference window, bounded if necessary by the start of the contig.
     *
     * @param start the start that is to be trimmed to the contig's start
     * @return The start, potentially trimmed to the contig's start
     */
    private int trimToContigStart(final int start) {
        return Math.max(start, 1);
    }

    /**
     * Determines the stop of the expanded reference window, bounded if necessary by the contig.
     *
     * @param locus The locus to expand.
     * @param windowTrailingBases number of bases to attempt to expand relative to the locus end (>= 0)
     * @return The end of the expanded window.
     */
    private int calculateWindowStop( final SimpleInterval locus, final int windowTrailingBases ) {
        return trimToContigLength(locus.getContig(), locus.getEnd() + windowTrailingBases);
    }

    /**
     * Determines the end of the expanded reference window, bounded if necessary by the contig.
     *
     * @param contig contig on which the location is.
     * @param end the end that is to be trimmed to the contig's length
     * @return The end, potentially trimmed to the contig's length
     */
    private int trimToContigLength(final String contig, final int end){

        final SAMSequenceRecord sequence = dataSource.getSequenceDictionary().getSequence(contig);
        if ( sequence == null ) {
            throw new UserException("Given reference file does not have data at the requested contig(" + contig + ")!");
        }

        final int sequenceLength = dataSource.getSequenceDictionary().getSequence(contig).getSequenceLength();
        return Math.min(end, sequenceLength);
    }

    /**
     * @param contig
     * @return the length/end position of the contig
     */
    private int getContigLength(final String contig){
        return dataSource.getSequenceDictionary().getSequence(contig).getSequenceLength();
    }

    /**
     * Get the base at the given locus.
     * @return The base at the given locus from the reference.
     */
    public byte getBase() {
        return getBases()[interval.getStart() - window.getStart()];
    }

    /**
     * Get a kmer around a position in reference without altering the internal state of the object
     * The position must lie within the window
     *
     * Returns null when, at the ends of a contig, we cannot expand the window to the requested size
     */
    public String getKmerAround(final int center, final int numBasesOnEachSide){
        Utils.validateArg(center >= 1, () -> "start position must be positive");
        Utils.validateArg(window.getStart() <= center && center <= window.getEnd(), "position must be smaller than end position");

        final SimpleInterval newWindow = new SimpleInterval(window.getContig(), center, center)
                .expandWithinContig(numBasesOnEachSide, getContigLength(window.getContig()));

        if (newWindow.getEnd() - newWindow.getStart() < 2*numBasesOnEachSide){
            return null;
        }

        return new String(getBases(newWindow));
    }
}
