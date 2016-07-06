package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.iterators.SAMRecordToReadIterator;
import org.broadinstitute.hellbender.utils.iterators.SamReaderQueryingIterator;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadConstants;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Manages traversals and queries over sources of reads (for now, SAM/BAM/CRAM files only).
 *
 * Two basic operations are available:
 *
 * -Iteration over all reads, optionally restricted to reads that overlap a set of intervals
 * -Targeted queries by one interval at a time
 */
public final class ReadsDataSource implements GATKDataSource<GATKRead>, AutoCloseable {
    protected static final Logger logger = LogManager.getLogger(ReadsDataSource.class);

    /**
     * Mapping from SamReaders to iterators over the reads from each reader. Only one
     * iterator can be open from a given reader at a time (this is a restriction
     * in htsjdk). Iterator is set to null for a reader if no iteration is currently
     * active on that reader.
     */
    private final Map<SamReader, CloseableIterator<SAMRecord>> readers;

    /**
     * Hang onto the input files so that we can print useful errors about them
     */
    private final Map<SamReader, File> backingFiles;

    /**
     * Only reads that overlap these intervals (and unmapped reads, if {@link #traverseUnmapped} is set) will be returned
     * during a full iteration. Null if iteration is unbounded.
     *
     * Individual queries are unaffected by these intervals -- only traversals initiated via {@link #iterator} are affected.
     */
    private List<SimpleInterval> intervalsForTraversal;

    /**
     * If true, restrict traversals to unmapped reads (and reads overlapping any {@link #intervalsForTraversal}, if set).
     * False if iteration is unbounded or bounded only by our {@link #intervalsForTraversal}.
     *
     * Note that this setting covers only unmapped reads that have no position -- unmapped reads that are assigned the
     * position of their mates will be returned by queries overlapping that position.
     *
     * Individual queries are unaffected by this setting  -- only traversals initiated via {@link #iterator} are affected.
     */
    private boolean traverseUnmapped;

    /**
     * Used to create a merged Sam header when we're dealing with multiple readers. Null if we only have a single reader.
     */
    private final SamFileHeaderMerger headerMerger;

    /**
     * Are indices available for all files?
     */
    private boolean indicesAvailable;

    /**
     * Initialize this data source with a single SAM/BAM file without a reference and validation stringency SILENT.
     *
     * @param samFile SAM/BAM file, not null.
     */
    public ReadsDataSource( final File samFile ) {
        this(samFile != null ? Arrays.asList(samFile) : null, null);
    }

    /**
     * Initialize this data source with multiple SAM/BAM files without a reference and validation stringency SILENT.
     *
     * @param samFiles SAM/BAM files, not null.
     */
    public ReadsDataSource( final List<File> samFiles ) {
        this(samFiles, null);
    }

    /**
     * Initialize this data source with a single SAM/BAM file and a custom SamReaderFactory
     *
     * @param samFile SAM/BAM file, not null.
     * @param customSamReaderFactory SamReaderFactory to use, if null a default factory with no reference and validation
     *                               stringency SILENT is used.
     */
    public ReadsDataSource( final File samFile, SamReaderFactory customSamReaderFactory) {
        this(samFile != null ? Arrays.asList(samFile) : null, customSamReaderFactory);
    }

    /**
     * Initialize this data source with multiple SAM/BAM files and a custom SamReaderFactory
     *
     * @param samFiles SAM/BAM files, not null.
     * @param customSamReaderFactory SamReaderFactory to use, if null a default factory with no reference and validation
     *                               stringency SILENT is used.
     */
    public ReadsDataSource( final List<File> samFiles, SamReaderFactory customSamReaderFactory ) {
        Utils.nonNull(samFiles);
        Utils.nonEmpty(samFiles, "ReadsDataSource cannot be created from empty file list");

        readers = new LinkedHashMap<>(samFiles.size() * 2);
        backingFiles = new LinkedHashMap<>(samFiles.size() * 2);
        indicesAvailable = true;

        final SamReaderFactory samReaderFactory =
                customSamReaderFactory == null ?
                    SamReaderFactory.makeDefault().validationStringency(ReadConstants.DEFAULT_READ_VALIDATION_STRINGENCY) :
                    customSamReaderFactory;

        for ( final File samFile : samFiles ) {
            // Ensure each file can be read
            try {
                IOUtil.assertFileIsReadable(samFile);
            }
            catch ( SAMException|IllegalArgumentException e ) {
                throw new UserException.CouldNotReadInputFile(samFile, e);
            }

            SamReader reader = samReaderFactory.open(samFile);

            // Ensure that each file has an index
            if ( ! reader.hasIndex() ) {
                indicesAvailable = false;
            }

            readers.put(reader, null);
            backingFiles.put(reader, samFile);
        }

        // Prepare a header merger only if we have multiple readers
        headerMerger = samFiles.size() > 1 ? createHeaderMerger() : null;
    }

    /**
     * Restricts a traversal of this data source via {@link #iterator} to only return reads which overlap the given intervals.
     * Calls to {@link #query} are not affected by setting these intervals.
     *
     * @param intervals Our next full traversal will return only reads overlapping these intervals
     */
    public void setTraversalBounds( final List<SimpleInterval> intervals ) {
        setTraversalBounds(intervals, false);
    }

    /**
     * Restricts a traversal of this data source via {@link #iterator} to only return reads that overlap the given intervals,
     * and to unmapped reads if specified.
     *
     * Calls to {@link #query} are not affected by this method.
     *
     * @param traversalParameters set of traversal parameters to control which reads get returned by the next call
     *                            to {@link #iterator}
     */
    public void setTraversalBounds( final TraversalParameters traversalParameters ) {
        setTraversalBounds(traversalParameters.getIntervalsForTraversal(), traversalParameters.traverseUnmappedReads());
    }

    /**
     * Restricts a traversal of this data source via {@link #iterator} to only return reads that overlap the given intervals,
     * and to unmapped reads if specified.
     *
     * Calls to {@link #query} are not affected by this method.
     *
     * @param intervals Our next full traversal will return reads overlapping these intervals
     * @param traverseUnmapped Our next full traversal will return unmapped reads (this affects only unmapped reads that
     *                         have no position -- unmapped reads that have the position of their mapped mates will be
     *                         included if the interval overlapping that position is included).
     */
    public void setTraversalBounds( final List<SimpleInterval> intervals, final boolean traverseUnmapped ) {
        // Set intervalsForTraversal to null if intervals is either null or empty
        this.intervalsForTraversal = intervals != null && ! intervals.isEmpty() ? intervals : null;
        this.traverseUnmapped = traverseUnmapped;

        if ( traversalIsBounded() && ! indicesAvailable ) {
            raiseExceptionForMissingIndex("Traversal by intervals was requested but some input files are not indexed.");
        }
    }

    /**
     * @return True if traversals initiated via {@link #iterator} will be restricted to reads that overlap intervals
     *         as configured via {@link #setTraversalBounds}, otherwise false
     */
    public boolean traversalIsBounded() {
        return intervalsForTraversal != null || traverseUnmapped;
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
     * Iterate over all reads in this data source. If intervals were provided via {@link #setTraversalBounds},
     * iteration is limited to reads that overlap that set of intervals.
     *
     * @return An iterator over the reads in this data source, limited to reads that overlap the intervals supplied
     *         via {@link #setTraversalBounds} (if intervals were provided)
     */
    @Override
    public Iterator<GATKRead> iterator() {
        logger.debug("Preparing readers for traversal");
        return prepareIteratorsForTraversal(intervalsForTraversal, traverseUnmapped);
    }

    /**
     * Query reads over a specific interval. This operation is not affected by prior calls to
     * {@link #setTraversalBounds}
     *
     * @param interval The interval over which to query
     * @return Iterator over reads overlapping the query interval
     */
    @Override
    public Iterator<GATKRead> query( final SimpleInterval interval ) {
        if ( ! indicesAvailable ) {
            raiseExceptionForMissingIndex("Cannot query reads data source by interval unless all files are indexed");
        }

        return prepareIteratorsForTraversal(Arrays.asList(interval));
    }

    /**
     * @return An iterator over just the unmapped reads with no assigned position. This operation is not affected
     *         by prior calls to {@link #setTraversalBounds}. The underlying file must be indexed.
     */
    public Iterator<GATKRead> queryUnmapped() {
        if ( ! indicesAvailable ) {
            raiseExceptionForMissingIndex("Cannot query reads data source by interval unless all files are indexed");
        }

        return prepareIteratorsForTraversal(null, true);
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
    private Iterator<GATKRead> prepareIteratorsForTraversal( final List<SimpleInterval> queryIntervals ) {
        return prepareIteratorsForTraversal(queryIntervals, false);
    }

    /**
     * Prepare iterators over all readers in response to a request for a complete iteration or query
     *
     * @param queryIntervals Intervals to bound the iteration (reads must overlap one of these intervals). If null, iteration is unbounded.
     * @return Iterator over all reads in this data source, limited to overlap with the supplied intervals
     */
    private Iterator<GATKRead> prepareIteratorsForTraversal( final List<SimpleInterval> queryIntervals, final boolean queryUnmapped ) {
        // htsjdk requires that only one iterator be open at a time per reader, so close out
        // any previous iterations
        closePreviousIterationsIfNecessary();

        final boolean traversalIsBounded = (queryIntervals != null && ! queryIntervals.isEmpty()) || queryUnmapped;

        // Set up an iterator for each reader, bounded to overlap with the supplied intervals if there are any
        for ( Map.Entry<SamReader, CloseableIterator<SAMRecord>> readerEntry : readers.entrySet() ) {
            readerEntry.setValue(traversalIsBounded ? new SamReaderQueryingIterator(readerEntry.getKey(), queryIntervals, queryUnmapped) :
                                                      readerEntry.getKey().iterator());
        }

        // Create a merging iterator over all readers if necessary. In the case where there's only a single reader,
        // return its iterator directly to avoid the overhead of the merging iterator.
        Iterator<SAMRecord> startingIterator = null;
        if ( readers.size() == 1 ) {
            startingIterator = readers.entrySet().iterator().next().getValue();
        }
        else {
            startingIterator = new MergingSamRecordIterator(headerMerger, readers, true);
        }

        return new SAMRecordToReadIterator(startingIterator);
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
