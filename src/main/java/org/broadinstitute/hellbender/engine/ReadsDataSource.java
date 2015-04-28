package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Manages traversals and queries over sources of reads (for now, SAM/BAM files only).
 *
 * Two basic operations are available:
 *
 * -Iteration over all reads, optionally restricted to reads that overlap a set of intervals
 * -Targeted queries by one interval at a time
 */
public final class ReadsDataSource implements GATKDataSource<SAMRecord>, AutoCloseable {
    protected static final Logger logger = LogManager.getLogger(ReadsDataSource.class);

    /**
     * Mapping from SamReaders to iterators over the reads from each reader. Only one
     * iterator can be open from a given reader at a time (this is a restriction
     * in htsjdk). Iterator is set to null for a reader if no iteration is currently
     * active on that reader.
     */
    private Map<SamReader, CloseableIterator<SAMRecord>> readers;

    /**
     * Hang onto the input files so that we can print useful errors about them
     */
    final private Map<SamReader, File> backingFiles;

    /**
     * Interval set to bound iteration over this data source. Null if iteration is unbounded.
     * Only reads that overlap these intervals will be returned during a full iteration.
     * Individual queries are unaffected by these intervals, however.
     */
    private List<SimpleInterval> intervals;

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
     * Are indices available for all files?
     */
    private boolean indicesAvailable;



    /**
     * Initialize this data source with a single SAM/BAM file
     *
     * @param samFile SAM/BAM file, not null.
     */
    public ReadsDataSource( final File samFile ) {
        this(samFile != null ? Arrays.asList(samFile) : null);
    }

    /**
     * Initialize this data source with multiple SAM/BAM files
     *
     * @param samFiles SAM/BAM files, not null.
     */
    public ReadsDataSource( final List<File> samFiles ) {
        if ( samFiles == null || samFiles.size() == 0 ) {
            throw new IllegalArgumentException("ReadsDataSource cannot be created from empty file list");
        }

        readers = new LinkedHashMap<>(samFiles.size() * 2);
        backingFiles = new LinkedHashMap<>(samFiles.size() *2);
        indicesAvailable = true;

        for ( File samFile : samFiles ) {
            // Ensure each file can be read
            try {
                IOUtil.assertFileIsReadable(samFile);
            }
            catch ( SAMException|IllegalArgumentException e ) {
                throw new UserException.CouldNotReadInputFile(samFile, e);
            }

            // TODO: allow SamReader settings to be customized by the client
            SamReader reader = SamReaderFactory.makeDefault().validationStringency(getValidationStringency()).open(samFile);

            // Ensure that each file has an index
            if ( ! reader.hasIndex() ) {
                indicesAvailable = false;
            }

            readers.put(reader, null);
            backingFiles.put(reader, samFile);
        }

        // Prepare a header merger only if we have multiple readers
        headerMerger = samFiles.size() > 1 ? createHeaderMerger() : null;
        mergingIterator = null;
    }

    private ValidationStringency getValidationStringency() {
        return ValidationStringency.SILENT;
    }

    /**
     * Restricts a traversal of this data source via {@link #iterator} to only return reads which overlap the given intervals.
     * Calls to {@link #query} are not affected by setting these intervals.
     *
     * @param intervals Our next full traversal will return only reads overlapping these intervals
     */
    public void setIntervalsForTraversal( final List<SimpleInterval> intervals ){
        // Treat null and empty interval lists the same
        this.intervals = (intervals != null && ! intervals.isEmpty()) ? intervals : null;

        if ( this.intervals != null ) {
            if ( ! indicesAvailable ) {
                raiseExceptionForMissingIndex("Traversal by intervals was requested but some input files are not indexed.");
            }

            logger.info("Preparing intervals for traversal");
            preparedIntervals = prepareIntervalsForTraversal();
            logger.info("Done preparing intervals for traversal");
        }
        else {
            preparedIntervals = null;
        }
    }

    private void raiseExceptionForMissingIndex(String reason) {
        String commandsToIndex = backingFiles.entrySet().stream()
                .filter(f -> !f.getKey().hasIndex())
                .map(Map.Entry::getValue)
                .map(File::getAbsolutePath)
                .map(f -> "samtools index " + f)
                .collect(Collectors.joining("\n","\n","\n"));

        throw new UserException(reason + "\nPlease index all input files:\n" + commandsToIndex);
    }

    /**
     * Iterate over all reads in this data source. If intervals were provided via {@link #setIntervalsForTraversal(List)},
     * iteration is limited to reads that overlap that set of intervals.
     *
     * @return An iterator over the reads in this data source, limited to reads that overlap the intervals supplied via {@link #setIntervalsForTraversal(List)} (if intervals were provided)
     */
    @Override
    public Iterator<SAMRecord> iterator() {
        logger.info("Preparing readers for traversal");
        final Iterator<SAMRecord> traversalIter = prepareIteratorsForTraversal(preparedIntervals);
        logger.info("Done preparing readers for traversal");

        return traversalIter;
    }

    /**
     * Query reads over a specific interval. This operation is not affected by the intervals supplied at construction.
     *
     * @param interval The interval over which to query
     * @return Iterator over reads overlapping the query interval
     */
    @Override
    public Iterator<SAMRecord> query( final SimpleInterval interval ) {
        if ( ! indicesAvailable )
            raiseExceptionForMissingIndex("Cannot query reads data source by interval unless all files are indexed");

        final QueryInterval[] queryInterval = { convertIntervalToQueryInterval(interval) };
        return prepareIteratorsForTraversal(queryInterval);
    }

    /**
     * Returns the SAM header for this data source. Will be a merged header if there are multiple readers.
     * If there is only a single reader, returns its header directly.
     *
     * @return SAM header for this data source
     */
    public SAMFileHeader getHeader() {
        return headerMerger != null ? headerMerger.getMergedHeader() : readers.entrySet().iterator().next().getKey().getFileHeader();
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
     * Converts our intervals from GATK format into htsjdk-compatible "QueryInterval" format suitable for
     * querying overlapping reads
     *
     * @return intervals converted into htsjdk QueryInterval format
     */
    private QueryInterval[] prepareIntervalsForTraversal() {
        QueryInterval[] convertedIntervals = new QueryInterval[intervals.size()];

        // Convert each SimpleInterval to a QueryInterval
        int intervalIndex = 0;
        for ( SimpleInterval interval : intervals ) {
            convertedIntervals[intervalIndex] = convertIntervalToQueryInterval(interval);
            ++intervalIndex;
        }

        // Intervals must be optimized (sorted and merged) in order to use the htsjdk query API
        return QueryInterval.optimizeIntervals(convertedIntervals);
    }

    /**
     * Converts an interval in SimpleInterval format into an htsjdk QueryInterval.
     *
     * In doing so, a header lookup is performed to convert from contig name to index
     *
     * @param interval interval to convert
     * @return an equivalent interval in QueryInterval format
     */
    private QueryInterval convertIntervalToQueryInterval( final SimpleInterval interval ) {
        final int contigIndex = getHeader().getSequenceIndex(interval.getContig());
        if ( contigIndex == -1 ) {
            throw new UserException("Contig " + interval.getContig() + " not present in reads sequence dictionary");
        }

        return new QueryInterval(contigIndex, interval.getStart(), interval.getEnd());
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

    /**
     * Get the sequence dictionary for this ReadsDataSource
     *
     * @return SAMSequenceDictionary from the SAMReader backing this if there is only 1 input file, otherwise the merged SAMSequenceDictionary from the merged header
     */
    public SAMSequenceDictionary getSequenceDictionary() {
        return getHeader().getSequenceDictionary();
    }

}
