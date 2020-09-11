package org.broadinstitute.hellbender.utils.iterators;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloseableIterator;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.TraversalParameters;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

/**
 * An iterator that allows for traversals over a SamReader restricted to a set of intervals, unmapped reads,
 * or both. Reads overlapping the provided set of intervals will be returned first (if requested), followed
 * by the unmapped reads with no position assigned (if requested).
 *
 * Note that we can't use something like com.google.common.collect.Iterators.concat() for this purpose
 * because of the restriction in the SamReader API that "Only a single open iterator on a given SAMFileReader
 * may be extant at any one time.  If you want to start a second iteration, the first one must be closed first."
 * This custom iterator ensures that this invariant is maintained, and also handles conversion of the intervals
 * to htsjdk query format.
 */
public class SamReaderQueryingIterator implements CloseableIterator<SAMRecord>, Iterable<SAMRecord> {
    protected static final Logger logger = LogManager.getLogger(SamReaderQueryingIterator.class);

    private final SamReader reader;
    private final QueryInterval[] queryIntervals;
    private final boolean queryUnmapped;
    private CloseableIterator<SAMRecord> currentIterator;
    private SAMRecord nextRecord;
    private boolean intervalQueryPerformed;
    private boolean unmappedQueryPerformed;

    /**
     * Create a SamReaderQueryingIterator given a SamReader, query intervals, and a flag indicating whether unmapped
     * reads with no position assigned should be returned. The query intervals may be null or empty, indicating that no
     * mapped reads are to be returned. Unmapped reads assigned a position will be returned as part of queries
     * for intervals that overlap the assigned position. Unmapped reads with no position are returned last, and
     * only if the queryUnmapped flag is set to true.
     *
     * The queryIntervals will be optimized internally using {@link QueryInterval#optimizeIntervals}, as required by htsjdk,
     * but this will not alter the List passed in.
     *
     * @param reader Source of reads
     * @param queryIntervals Return reads overlapping these intervals first. May be null or empty, indicating that no
     *                       mapped reads are to be returned. Will return unmapped reads that are assigned a position
     *                       if the assigned position overlaps these intervals.
     * @param queryUnmapped If true, return reads that are unmapped and don't have a position assigned. These will be returned
     *                      after any reads that overlap the queryIntervals, if provided.
     */
    public SamReaderQueryingIterator( final SamReader reader, final List<SimpleInterval> queryIntervals, final boolean queryUnmapped ) {
        Utils.nonNull(reader);

        this.reader = reader;
        this.queryIntervals = prepareQueryIntervals(queryIntervals);
        this.queryUnmapped = queryUnmapped;
        this.intervalQueryPerformed = false;
        this.unmappedQueryPerformed = false;

        this.currentIterator = loadNextIterator();
        this.nextRecord = loadNextRecord();
    }

    /**
     * Create a SamReaderQueryingIterator given a SamReader and a set of traversal parameters governing
     * which reads are to be returned during iteration. See SamReaderQueryingIterator(SamReader, List<SimpleInterval>, boolean)
     * for details on how the settings affect traversal.
     *
     * @param reader Source of reads
     * @param traversalParameters Traversal parameters governing which reads are to be returned during iteration
     */
    public SamReaderQueryingIterator( final SamReader reader, final TraversalParameters traversalParameters ) {
        this(reader, traversalParameters.getIntervalsForTraversal(), traversalParameters.traverseUnmappedReads());
    }

    /**
     * Converts a List of SimpleIntervals into the format required by the SamReader query API
     * @param rawIntervals SimpleIntervals to be converted
     * @return A sorted, merged list of QueryIntervals suitable for passing to the SamReader query API
     */
    private QueryInterval[] prepareQueryIntervals( final List<SimpleInterval> rawIntervals ) {
        if ( rawIntervals == null || rawIntervals.isEmpty() ) {
            return null;
        }

        // This might take a while with large interval lists, so log a status message
        logger.debug("Preparing intervals for traversal");

        // Convert each SimpleInterval to a QueryInterval
        final QueryInterval[] convertedIntervals =
                rawIntervals.stream()
                .map(rawInterval -> IntervalUtils.convertSimpleIntervalToQueryInterval(rawInterval, reader.getFileHeader().getSequenceDictionary()))
                .toArray(QueryInterval[]::new);

        // Intervals must be optimized (sorted and merged) in order to use the htsjdk query API
        return QueryInterval.optimizeIntervals(convertedIntervals);
    }

    private SAMRecord loadNextRecord() {
        if ( currentIterator == null ) {
            return null;
        }

        while ( currentIterator != null && ! currentIterator.hasNext() ) {
            currentIterator = loadNextIterator();
        }

        return currentIterator != null ? currentIterator.next() : null;
    }

    private CloseableIterator<SAMRecord> loadNextIterator() {
        // The SamReader API requires us to close out the previous iterator over our reader before opening a new one.
        if ( currentIterator != null ) {
            currentIterator.close();
        }

        if ( hasQueryIntervals() && ! intervalQueryPerformed ) {
            intervalQueryPerformed = true;
            return reader.queryOverlapping(queryIntervals);
        }
        else if ( queryUnmapped && ! unmappedQueryPerformed ) {
            unmappedQueryPerformed = true;
            return reader.queryUnmapped();
        }

        return null;
    }

    private boolean hasQueryIntervals() {
        return queryIntervals != null && queryIntervals.length > 0;
    }

    @Override
    public boolean hasNext() {
        return nextRecord != null;
    }

    @Override
    public SAMRecord next() {
        if ( nextRecord == null ) {
            throw new NoSuchElementException("Iterator is exhausted");
        }

        final SAMRecord toReturn = nextRecord;
        nextRecord = loadNextRecord();
        return toReturn;
    }

    @Override
    public void close() {
        if ( currentIterator != null ) {
            currentIterator.close();
        }
    }

    @Override
    public Iterator<SAMRecord> iterator() {
        return this;
    }
}
