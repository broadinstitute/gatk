package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.GenomeLoc;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Manages traversals and queries over sources of reads (for now, SAM/BAM files only).
 *
 * Two basic operations are available:
 *
 * -Iteration over all reads, optionally restricted to reads that overlap a set of intervals
 * -Targeted queries by one interval at a time
 */
public class ReadsDataSource implements GATKDataSource<SAMRecord>, AutoCloseable {

    /**
     * Mapping from SamReaders to iterators over the reads from each reader. Only one
     * iterator can be open from a given reader at a time (this is a restriction
     * in htsjdk). Iterator is set to null for a reader if no iteration is currently
     * active on that reader.
     */
    private Map<SamReader, CloseableIterator<SAMRecord>> readers;

    /**
     * Interval set to bound iteration over this data source. Null if iteration is unbounded.
     * Only reads that overlap these intervals will be returned during a full iteration.
     * Individual queries are unaffected by these intervals, however.
     */
    private List<GenomeLoc> intervals;

    /**
     * Bounding intervals in htsjdk-compatible form after processing to sort and merge adjacent/overlapping
     * intervals (as required by htsjdk). Null if iteration is unbounded. Only reads that overlap these
     * intervals will be returned during a full iteration. Individual queries are unaffected by these
     * intervals, however.
     */
    private QueryInterval[] preparedIntervals;

    /**
     * Used to create a merged Sam header when we're dealing with multiple readers. Null if we only have a single reader.
     */
    private SamFileHeaderMerger headerMerger;

    /**
     * Merges iterators from multiple readers into a single sorted stream of reads. Null if we only have a single reader.
     */
    private MergingSamRecordIterator mergingIterator;

    /**
     * Initialize this data source with a single SAM/BAM file, and no intervals to restrict iteration
     *
     * @param samFile SAM/BAM file, not null. Must be indexed.
     */
    public ReadsDataSource( final File samFile ) {
        this(samFile != null ? Arrays.asList(samFile) : null, null);
    }

    /**
     * Initialize this data source with multiple SAM/BAM files, and no intervals to restrict iteration
     *
     * @param samFiles SAM/BAM files, not null. Each file must be indexed.
     */
    public ReadsDataSource( final List<File> samFiles ) {
        this(samFiles, null);
    }

    /**
     * Initialize this data source with a single SAM/BAM file, and a set of intervals to bound iteration.
     * Only reads that overlap these intervals will be returned during a full iteration (but individual
     * queries will be unaffected by these intervals).
     *
     * @param samFile SAM/BAM file, not null. Must be indexed.
     * @param intervals Intervals for iteration over this data source. Do not affect individual queries.
     */
    public ReadsDataSource( final File samFile, final List<GenomeLoc> intervals ) {
        this(samFile != null ? Arrays.asList(samFile) : null, intervals);
    }

    /**
     * Initialize this data source with multiple SAM/BAM files, and a set of intervals to bound iteration.
     * Only reads that overlap these intervals will be returned during a full iteration (but individual
     * queries will be unaffected by these intervals).
     *
     * @param samFiles SAM/BAM files, not null. Each file must be indexed.
     * @param intervals Intervals for iteration over this data source. Do not affect individual queries.
     */
    public ReadsDataSource( final List<File> samFiles, final List<GenomeLoc> intervals ) {
        if ( samFiles == null || samFiles.size() == 0 ) {
            throw new IllegalArgumentException("ReadsDataSource cannot be created from empty file list");
        }

        readers = new LinkedHashMap<>(samFiles.size() * 2);

        for ( File samFile : samFiles ) {
            // Ensure each file can be read
            try {
                IOUtil.assertFileIsReadable(samFile);
            }
            catch ( SAMException|IllegalArgumentException e ) {
                throw new UserException.CouldNotReadInputFile(samFile, e);
            }

            // TODO: allow SamReader settings to be customized by the client
            SamReader reader = SamReaderFactory.makeDefault().open(samFile);

            // Ensure that each file has an index
            if ( ! reader.hasIndex() ) {
                throw new UserException("File " + samFile.getAbsolutePath() + " has no index. Please index this file (can be done using \"samtools index\")");
            }

            readers.put(reader, null);
        }

        // Treat null and empty interval lists the same
        this.intervals = (intervals != null && ! intervals.isEmpty()) ? intervals : null;
        preparedIntervals = this.intervals != null ? prepareIntervalsForTraversal() : null;

        // Prepare a header merger only if we have multiple readers
        headerMerger = samFiles.size() > 1 ? createHeaderMerger() : null;
        mergingIterator = null;
    }

    /**
     * Iterate over all reads in this data source. If intervals were provided at construction time,
     * iteration is limited to reads that overlap that set of intervals.
     *
     * @return An iterator over the reads in this data source, limited to reads that overlap the intervals supplied at construction (if intervals were provided)
     */
    @Override
    public Iterator<SAMRecord> iterator() {
        return prepareIteratorsForTraversal(preparedIntervals);
    }

    /**
     * Query reads over a specific interval. This operation is not affected by the intervals supplied at construction.
     *
     * @param interval The interval over which to query
     * @return Iterator over reads overlapping the query interval
     */
    @Override
    public Iterator<SAMRecord> query( final GenomeLoc interval ) {
        final QueryInterval[] queryInterval = { new QueryInterval(interval.getContigIndex(), interval.getStart(), interval.getStop()) };
        return prepareIteratorsForTraversal(queryInterval);
    }

    /**
     * Prepare iterators over all readers in response to a request for a complete iteration or query
     *
     * If there are multiple intervals, they must have been optimized using QueryInterval.optimizeIntervals()
     * before calling this method.
     *
     * @param queryIntervals Intervals to bound the iteration (reads must overlap one of these intervals). If null, iteration is unbounded.
     * @return Iterator over all reads in this data source, limited to overlap with the supplied intervals
     */
    private Iterator<SAMRecord> prepareIteratorsForTraversal( final QueryInterval[] queryIntervals ) {
        // htsjdk requires that only one iterator be open at a time per reader, so close out
        // any previous iterations
        closePreviousIterationsIfNecessary();

        // Set up an iterator for each reader, bounded to overlap with the supplied intervals if there are any
        for ( Map.Entry<SamReader, CloseableIterator<SAMRecord>> readerEntry : readers.entrySet() ) {
            readerEntry.setValue(queryIntervals == null ? readerEntry.getKey().iterator() :
                                                          readerEntry.getKey().queryOverlapping(queryIntervals));
        }

        // Create a merging iterator over all readers if necessary. In the case where there's only a single reader,
        // return its iterator directly to avoid the overhead of the merging iterator.
        Iterator<SAMRecord> startingIterator = null;
        if ( readers.size() == 1 ) {
            startingIterator = readers.entrySet().iterator().next().getValue();
        }
        else {
            startingIterator= new MergingSamRecordIterator(headerMerger, readers, true);
        }

        // Apply any additional transformations on the read stream before returning it
        return applyDecoratingIterators(startingIterator);
    }

    /**
     * Apply arbitrary transformations to the read stream (such as read filtering, downsampling, etc.)
     *
     * @param startingIterator basic iterator over all reads in this data source
     * @return an iterator (or chain of iterators) that wraps the startingIterator to transform the read stream in some way
     */
    private Iterator<SAMRecord> applyDecoratingIterators( final Iterator<SAMRecord> startingIterator ) {
        // For now, just return the iterator we are given. But, in the future, if we want
        // to inject additional iterators into the read stream to modify it in some way,
        // this is the place to do it.
        return startingIterator;
    }

    /**
     * Converts our intervals from GATK format into htsjdk-compatible format suitable for querying
     * overlapping reads
     *
     * @return intervals converted into htsjdk format
     */
    private QueryInterval[] prepareIntervalsForTraversal() {
        QueryInterval[] convertedIntervals = new QueryInterval[intervals.size()];

        // Convert each GenomeLoc to a QueryInterval
        int intervalIndex = 0;
        for ( GenomeLoc interval : intervals ) {
            convertedIntervals[intervalIndex] = new QueryInterval(interval.getContigIndex(), interval.getStart(), interval.getStop());
            ++intervalIndex;
        }

        // Intervals must be optimized (sorted and merged) in order to use the htsjdk query API
        return QueryInterval.optimizeIntervals(convertedIntervals);
    }

    /**
     * Create a header merger from the individual SAM/BAM headers in our readers
     *
     * @return a header merger containing all individual headers in this data source
     */
    private SamFileHeaderMerger createHeaderMerger() {
        List<SAMFileHeader> headers = new ArrayList<>(readers.size());
        for ( Map.Entry<SamReader, CloseableIterator<SAMRecord>> readerEntry : readers.entrySet() ) {
            headers.add(readerEntry.getKey().getFileHeader());
        }

        // TODO: don't require coordinate ordering
        SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(SAMFileHeader.SortOrder.coordinate, headers, true);
        return headerMerger;
    }

    /**
     * Shut down this data source permanently, closing all iterations and readers.
     */
    public void close() {
        closePreviousIterationsIfNecessary();

        try {
            for ( Map.Entry<SamReader, CloseableIterator<SAMRecord>> readerEntry : readers.entrySet() ) {
                readerEntry.getKey().close();
            }
        }
        catch ( IOException e ) {
            throw new GATKException("Error closing SAMReader");
        }
    }

    /**
     * Close any previously-opened iterations over our readers (htsjdk allows only one open iteration per reader).
     */
    private void closePreviousIterationsIfNecessary() {
        if ( mergingIterator != null ) {
            mergingIterator.close();
            mergingIterator = null;
        }

        for ( Map.Entry<SamReader, CloseableIterator<SAMRecord>> readerEntry : readers.entrySet() ) {
            CloseableIterator<SAMRecord> readerIterator = readerEntry.getValue();
            if ( readerIterator != null ) {
                readerIterator.close();
                readerEntry.setValue(null);
            }
        }
    }
}
