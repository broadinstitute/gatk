package org.broadinstitute.hellbender.engine;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.util.concurrent.ThreadFactoryBuilder;
import htsjdk.samtools.*;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.DelegatingIterator;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.htsgetreader.*;
import org.broadinstitute.hellbender.transformers.SamRecordTransformer;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.iterators.ReadFilteringIterator;
import org.broadinstitute.hellbender.utils.iterators.SAMRecordToReadIterator;
import org.broadinstitute.hellbender.utils.iterators.SamRecordTransformingIterator;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadConstants;

import javax.annotation.Nonnull;
import java.io.*;
import java.nio.file.Path;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

/**
 * Manages traversals and queries over sources of reads which are accessible via {@link Path}s pointing to a file
 * behind an htsget server
 * (for now, SAM/BAM/CRAM files only).
 *
 * Two basic operations are available:
 *
 * -Iteration over all reads, optionally restricted to reads that overlap a set of intervals
 * -Targeted queries by one interval at a time
 */
public final class ReadsHtsgetDataSource implements ReadsDataSource {
    private static final Logger logger = LogManager.getLogger(ReadsHtsgetDataSource.class);

    /**
     * The sources provided to this data source
     */
    private final List<GATKPath> sources;

    /**
     * Only reads that overlap these intervals (and unmapped reads, if {@link #traverseUnmapped} is set) will be returned
     * during a full iteration. Null if iteration is unbounded.
     *
     * Individual queries are unaffected by these intervals -- only traversals initiated via {@link #iterator} are affected.
     */
    private List<SimpleInterval> intervals;

    /**
     * If true, restrict traversals to unmapped reads (and reads overlapping any {@link #intervals}, if set).
     * False if iteration is unbounded or bounded only by our {@link #intervals}.
     *
     * Note that this setting covers only unmapped reads that have no position -- unmapped reads that are assigned the
     * position of their mates will be returned by queries overlapping that position.
     *
     * Individual queries are unaffected by this setting  -- only traversals initiated via {@link #iterator} are affected.
     */
    private boolean traverseUnmapped;

    /**
     * Mapping from the paths provided to this data source and their SAM headers
     */
    private final Map<GATKPath, SAMFileHeader> headers;

    /**
     * Used to create a merged Sam header when we're dealing with multiple readers. Null if we only have a single reader.
     */
    private final SamFileHeaderMerger samHeaderMerger;

    /**
     * The merged header for this data source
     */
    private final SAMFileHeader header;

    /**
     * SamReaderFactory used to create the SamReaders used by this data source
     */
    private final SamReaderFactory readerFactory;

    /**
     * The inner iterator this source returns when queried either by interval, for unmapped reads or a full
     * iteration is requested
     */
    private CloseableIterator<GATKRead> iterator;

    /**
     * Initialize this data source with a single SAM/BAM file and validation stringency SILENT.
     *
     * @param source SAM/BAM file, not null.
     */
    public ReadsHtsgetDataSource(final GATKPath source) {
        this(source, null);
    }

    /**
     * Initialize this data source with multiple SAM/BAM files and validation stringency SILENT.
     *
     * @param sources SAM/BAM files, not null.
     */
    public ReadsHtsgetDataSource(final List<GATKPath> sources) {
        this(sources, null);
    }

    /**
     * Initialize this data source with a single SAM/BAM file and a custom SamReaderFactory
     *
     * @param source path to SAM/BAM file, not null.
     * @param customSamReaderFactory SamReaderFactory to use, if null a default factory with no reference and validation
     *                               stringency SILENT is used.
     */
    public ReadsHtsgetDataSource(final GATKPath source, final SamReaderFactory customSamReaderFactory) {
        this(Collections.singletonList(source), customSamReaderFactory);
    }

    /**
     * Initialize this data source with multiple SAM/BAM files and a custom SamReaderFactory
     *
     * @param sources path to SAM/BAM file, not null.
     * @param customSamReaderFactory SamReaderFactory to use, if null a default factory with no reference and validation
     *                               stringency SILENT is used.
     */
    public ReadsHtsgetDataSource(final List<GATKPath> sources, final SamReaderFactory customSamReaderFactory) {
        Utils.nonNull(sources);
        Utils.nonEmpty(sources, "ReadsHtsgetDataSource cannot be created from empty source list");

        this.readerFactory =
            customSamReaderFactory == null
                ? SamReaderFactory.makeDefault().validationStringency(ReadConstants.DEFAULT_READ_VALIDATION_STRINGENCY)
                : customSamReaderFactory;

        this.sources = sources;
        this.headers = new HashMap<>(sources.size());

        for (final GATKPath source : sources) {
            try {
                // Ensure each source can be obtained from htsget server
                final HtsgetRequest req = new HtsgetRequest(source).withDataClass(HtsgetClass.header);
                // Request only the headers and use them to construct SAMFileHeader for each source
                final InputStream headerStream = getRequestInputStream(req);
                this.headers.put(source, readerFactory.open(SamInputResource.of(headerStream)).getFileHeader());
            } catch (final UserException e) {
                throw new UserException(source.toString(), e);
            }
        }

        this.samHeaderMerger = createHeaderMerger(this.headers.values());
        this.header = sources.size() > 1
            ? this.samHeaderMerger.getMergedHeader()
            : headers.values().iterator().next();
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
    @Override
    public void setTraversalBounds(final List<SimpleInterval> intervals, final boolean traverseUnmapped) {
        this.intervals = intervals != null && ! intervals.isEmpty() ? intervals : null;
        this.traverseUnmapped = traverseUnmapped;
    }

    /**
     * @return True if traversals initiated via {@link #iterator} will be restricted to reads that overlap intervals
     *         as configured via {@link #setTraversalBounds}, otherwise false
     */
    @Override
    public boolean traversalIsBounded() {
        return this.intervals != null || this.traverseUnmapped;
    }

    /**
     * This data source can be queried even without index files
     * @return always true
     */
    @Override
    public boolean isQueryableByInterval() {
        return true;
    }

    /**
     * Iterate over all reads in this data source. If intervals were provided via {@link #setTraversalBounds},
     * iteration is limited to reads that overlap that set of intervals.
     *
     * @return An iterator over the reads in this data source, limited to reads that overlap the intervals supplied
     *         via {@link #setTraversalBounds} (if intervals were provided)
     */
    @Override
    @Nonnull
    public Iterator<GATKRead> iterator() {
        this.iterator = new MergingHtsgetIterator(this.intervals, this.traverseUnmapped);
        return this.iterator;
    }

    /**
     * Query reads over a specific interval. This operation is not affected by prior calls to
     * {@link #setTraversalBounds}
     *
     * @param interval The interval over which to query
     * @return Iterator over reads overlapping the query interval
     */
    @Override
    public Iterator<GATKRead> query(final SimpleInterval interval) {
        this.iterator = new MergingHtsgetIterator(Collections.singletonList(interval), false);
        return this.iterator;
    }

    /**
     * @return An iterator over just the unmapped reads with no assigned position. This operation is not affected
     *         by prior calls to {@link #setTraversalBounds}. The underlying file must be indexed.
     */
    @Override
    public Iterator<GATKRead> queryUnmapped() {
        this.iterator = new MergingHtsgetIterator(null, true);
        return this.iterator;
    }

    @Override
    public boolean supportsSerialIteration() {
        return false;
    }

    /**
     * Shut down this data source permanently, closing all iterations and readers.
     */
    @Override
    public void close() {
        this.iterator.close();
    }

    /**
     * Returns the SAM header for this data source. Will be a merged header if there are multiple readers.
     * If there is only a single reader, returns its header directly.
     *
     * @return SAM header for this data source
     */
    @Override
    public SAMFileHeader getHeader() {
        return this.header;
    }

    private SamFileHeaderMerger createHeaderMerger(final Collection<SAMFileHeader> headers) {
        return new SamFileHeaderMerger(identifySortOrder(headers), headers, true);
    }

    private static SAMFileHeader.SortOrder identifySortOrder(final Collection<SAMFileHeader> headers) {
        final Set<SAMFileHeader.SortOrder> sortOrders = headers.stream()
            .map(SAMFileHeader::getSortOrder)
            .collect(Collectors.toSet());
        final SAMFileHeader.SortOrder order;
        if (sortOrders.size() == 1) {
            order = sortOrders.iterator().next();
        } else {
            order = SAMFileHeader.SortOrder.unsorted;
            logger.warn("Inputs have different sort orders. Assuming {} sorted reads for all of them.", order);
        }
        return order;
    }

    /**
     * Returns an InputStream over the data contained in the response to the given htsget request
     * Returns header requests immediately as these are usually small, and downloads larger requests to disk
     * @param req the htsget request to return an InputStream for
     * @return an InputStream over the contents of the request's data
     */
    @VisibleForTesting
    private static InputStream getRequestInputStream(final HtsgetRequest req) {
        return req.getResponse().getDataStream();
    }

    /**
     * Get the sequence dictionary for this ReadsHtsgetDataSource
     *
     * @return SAMSequenceDictionary from the SAMReader backing this if there is only 1 input file, otherwise the merged SAMSequenceDictionary from the merged header
     */
    @Override
    public SAMSequenceDictionary getSequenceDictionary() {
        return this.header.getSequenceDictionary();
    }

    /**
     * Manages iteration over the provided sources and intervals -- schedules downloads on background threads allowing
     * main thread to begin processing reads while background threads initialize readers
     */
    private class MergingHtsgetIterator implements CloseableIterator<GATKRead> {
        /**
         * Only reads that overlap these intervals (and unmapped reads, if {@link #traverseUnmapped} is set) will be returned
         * during a full iteration. Null if iteration is unbounded.
         */
        private final List<SimpleInterval> intervals;

        /**
         * If true, restrict traversals to unmapped reads (and reads overlapping any {@link #intervals}, if set).
         * False if iteration is unbounded or bounded only by our {@link #intervals}.
         *
         * Note that this setting covers only unmapped reads that have no position -- unmapped reads that are assigned the
         * position of their mates will be returned if an interval overlapping that position is provided.
         */
        private final boolean traverseUnmapped;

        /**
         * Executor used to schedule background download threads
         */
        private final ExecutorService executor;

        /**
         * Queue onto which background threads enqueue non-empty iterators and from which the main thread draws
         * once it has exhausted its current iterator -- uses a fixed capacity blocking queue to automatically apply
         * backpressure if background threads are producing or main thread is consuming too quickly
         */
        private final BlockingQueue<CloseableIterator<GATKRead>> queue;

        /**
         * Number of non-empty iterators the main thread expects to be able to draw from the queue
         * Initially set to total number of iterators by main thread -- this may include some empty iterators
         * Producing threads will not enqueue empty iterators, so they must decrement the counter to account for these
         * missing iterators the main thread is expecting
         * Main thread is responsible for decrementing counter to account for non-empty iterators it draws from the queue
         * When this reaches 0, every iterator from every source has either been taken and used up, or was empty
         */
        private final AtomicInteger sourcesWaiting;

        /**
         * The current iterator to draw from
         */
        private CloseableIterator<GATKRead> currentIterator;

        /**
         * Initialize this iterator to only return reads that overlap the given intervals and to unmapped reads if specified.
         * Instantiating this iterator begins downloading the specified sources
         * @param intervals This iterator will return reads overlapping these intervals
         * @param traverseUnmapped This iterator will return unmapped reads (this affects only unmapped reads that
         *                         have no position -- unmapped reads that have the position of their mapped mates will be
         *                         included if the interval overlapping that position is included).
         */
        public MergingHtsgetIterator(final List<SimpleInterval> intervals, final boolean traverseUnmapped) {
            this.executor = Executors.newFixedThreadPool(4, new ThreadFactoryBuilder()
                .setNameFormat("htsget-datasource-thread-%d")
                .setDaemon(true)
                .build());

            this.intervals = intervals;
            this.traverseUnmapped = traverseUnmapped;

            this.queue = new ArrayBlockingQueue<>(4);
            this.sourcesWaiting = new AtomicInteger(0);

            if (this.intervals == null) {
                this.startNoIntervals();
            } else {
                this.startIntervals();
            }
        }

        /**
         * Signal background threads to stop and close all existing iterators in the queue
         */
        @Override
        public void close() {
            this.executor.shutdownNow();
            this.queue.forEach(CloseableIterator::close);
        }

        /**
         * Return the next non-empty iterator, decrementing the count of sources waiting, or null if there are none
         * @return A non-empty iterator
         */
        private CloseableIterator<GATKRead> nextIterator() {
            try {
                CloseableIterator<GATKRead> iterator = null;
                while (iterator == null && this.sourcesWaiting.get() > 0) {
                    iterator = this.queue.poll(1, TimeUnit.SECONDS);
                }
                if (iterator != null) {
                    this.sourcesWaiting.decrementAndGet();
                }
                return iterator;
            } catch (final InterruptedException e) {
                throw new GATKException("Interrupted while dequeueing", e);
            }
        }

        @Override
        public boolean hasNext() {
            if (this.currentIterator != null && this.currentIterator.hasNext()) {
                return true;
            }
            if (this.currentIterator != null) this.currentIterator.close();
            this.currentIterator = this.nextIterator();
            return this.currentIterator != null;
        }

        @Override
        public GATKRead next() {
            return this.currentIterator.next();
        }

        /**
         * Begin fetching sources from htsget server on background threads -- downloads either the whole file if
         * complete iteration was requested, or only the unmapped reads if unmapped read iteration was requested
         */
        private void startNoIntervals() {
            for (final GATKPath source : sources) {
                this.sourcesWaiting.getAndIncrement();
                // Enqueue an iterator over either the entire source or just its unmapped reads
                this.executor.submit(this.traverseUnmapped
                    ? () -> this.getSourceUnmapped(source)
                    : () -> this.getSource(source));
            }
        }

        /**
         * Reduce the intervals down to only include ones that can actually intersect with this source
         * @param source A source provided to this data source
         * @return A list of intervals which overlap this source, based on its header
         */
        private List<SimpleInterval> intervalsOverlappingSource(final GATKPath source) {
            final SAMSequenceDictionary sequenceDictionary = headers.get(source).getSequenceDictionary();
            return this.intervals.stream()
                .filter(interval -> IntervalUtils.intervalIsOnDictionaryContig(interval, sequenceDictionary))
                .collect(Collectors.toList());
        }

        /**
         * Begin fetching sources from htsget server on background threads -- makes a single request per source
         * per overlapping interval, and possibly a request for the unmapped reads if this was requested
         */
        private void startIntervals() {
            for (final GATKPath source : sources) {
                SimpleInterval prevInterval = null;
                // Enqueue iterators over the intervals of this source
                for (final SimpleInterval currInterval : intervalsOverlappingSource(source)) {
                    final SimpleInterval finalPrevInterval = prevInterval;
                    this.sourcesWaiting.getAndIncrement();
                    this.executor.submit(() -> this.getSourceWithInterval(source, finalPrevInterval, currInterval));
                    prevInterval = currInterval;
                }
                // Enqueue an iterator over the unmapped reads only if traverseUnmapped is set to true
                if (this.traverseUnmapped) {
                    this.sourcesWaiting.getAndIncrement();
                    this.executor.submit(() -> this.getSourceUnmapped(source));
                }
            }
        }

        /**
         * Wrap an iterator to allow us to free the backing SamReader without holding onto an explicit reference to it
         * @param iterator the iterator to wrap
         * @param closeable the SamReader to close once the iterator has been used up
         * @return a wrapped CloseableIterator
         */
        private CloseableIterator<GATKRead> wrapIteratorWithClose(final Iterator<GATKRead> iterator, final Closeable closeable) {
            return new DelegatingIterator<GATKRead>(iterator) {
                @Override
                public void close() {
                    try {
                        closeable.close();
                    } catch (final IOException e) {
                        throw new GATKException("Error closing SAMReader", e);
                    }
                }
            };
        }

        /**
         * Wrap an iterator over SAMRecords to resolve its reference indices against this data source's merged header
         * @param samReader the reader to resolve
         * @return a wrapped iterator over SAMRrecords which resolves its records references
         */
        private Iterator<SAMRecord> wrapSamRecordWithMergedHeader(final SamReader samReader) {
            return new SamRecordTransformingIterator(samReader.iterator(), (SamRecordTransformer) record -> {
                // this will resolve the reference indices against the new, merged header
                record.setHeader(samHeaderMerger.getMergedHeader());

                // Fix the read group if needs be
                if (samHeaderMerger.hasReadGroupCollisions()) {
                    final String oldGroupId = (String) record.getAttribute(ReservedTagConstants.READ_GROUP_ID);
                    if (oldGroupId != null) {
                        final String newGroupId = samHeaderMerger.getReadGroupId(samReader.getFileHeader(), oldGroupId);
                        record.setAttribute(ReservedTagConstants.READ_GROUP_ID, newGroupId);
                    }
                }

                // Fix the program group if needs be
                if (samHeaderMerger.hasProgramGroupCollisions()) {
                    final String oldGroupId = (String) record.getAttribute(ReservedTagConstants.PROGRAM_GROUP_ID);
                    if (oldGroupId != null) {
                        final String newGroupId = samHeaderMerger.getProgramGroupId(samReader.getFileHeader(), oldGroupId);
                        record.setAttribute(ReservedTagConstants.PROGRAM_GROUP_ID, newGroupId);
                    }
                }

                return record;
            });
        }

        /**
         * Download all reads of the source and enqueue an iterator over them if the iterator is non-empty
         * @param source the source to download
         */
        private void getSource(final GATKPath source) {
            final HtsgetRequest req = new HtsgetRequest(source);
            final SamReader samReader = readerFactory.open(SamInputResource.of(getRequestInputStream(req)));
            // Resolve reference indices against merged header, fix read group/program id and return iterator
            final Iterator<SAMRecord> transformedReader = this.wrapSamRecordWithMergedHeader(samReader);
            final Iterator<GATKRead> readIterator = new SAMRecordToReadIterator(transformedReader);
            final CloseableIterator<GATKRead> closeWrappedIterator = this.wrapIteratorWithClose(readIterator, samReader);
            this.enqueueIfNotEmpty(closeWrappedIterator);
        }

        /**
         * Download only the unmapped reads of the source and enqueue an iterator over them if the iterator is non-empty
         * @param source the source to download
         */
        private void getSourceUnmapped(final GATKPath source) {
            final HtsgetRequest req = new HtsgetRequest(source).withInterval(HtsgetRequest.UNMAPPED_UNPLACED_INTERVAL);
            final SamReader samReader = readerFactory.open(SamInputResource.of(getRequestInputStream(req)));
            final Iterator<GATKRead> readIterator = new SAMRecordToReadIterator(samReader.iterator());
            final CloseableIterator<GATKRead> closeWrappedIterator = this.wrapIteratorWithClose(readIterator, samReader);
            this.enqueueIfNotEmpty(closeWrappedIterator);
        }

        /**
         * Get only the current interval from a source, excluding duplicate reads that overlap the previous interval
         * and enqueue an iterator over them if the iterator is non-empty
         * @param source the source to download
         * @param prevInterval the previous interval, if it exists, whose reads should be excluded from this interval to avoid duplication
         * @param currInterval the current interval to request
         */
        private void getSourceWithInterval(final GATKPath source, final SimpleInterval prevInterval, final SimpleInterval currInterval) {
            final HtsgetRequest req = new HtsgetRequest(source).withInterval(currInterval);
            final SamReader samReader = readerFactory.open(SamInputResource.of(getRequestInputStream(req)));
            // Resolve reference indices against merged header, fix read group/program id and return iterator
            final Iterator<SAMRecord> transformedReader = this.wrapSamRecordWithMergedHeader(samReader);
            // Convert SAMRecords to GATKReads
            final Iterator<GATKRead> readIterator = new SAMRecordToReadIterator(transformedReader);

            /*
            To remove reads duplicated across two subsequent intervals, we take any read which is in
            the current interval but NOT in the previous interval, unless the current interval is the first,
            in which case we take any read overlapping it
             */
            final Iterator<GATKRead> filteredIterator = new ReadFilteringIterator(
                readIterator,
                read -> currInterval.overlaps(read) &&
                    (prevInterval == null || !prevInterval.overlaps(read))
            );

            // Wrap in CloseableIterator so we can close backing SamReader without explicitly keeping a reference to it
            final CloseableIterator<GATKRead> closeWrappedIterator = this.wrapIteratorWithClose(filteredIterator, samReader);
            // Only enqueue if the iterator is non-empty
            this.enqueueIfNotEmpty(closeWrappedIterator);
        }

        /**
         * Enqueue an iterator onto our queue of iterators the main thread will draw from only if it is non-empty
         * or close it and decrement the count of sources waiting
         */
        private void enqueueIfNotEmpty(final CloseableIterator<GATKRead> iterator) {
            if (!iterator.hasNext()) {
                iterator.close();
                /*
                We need to decrement our count of sources waiting to be submitted because the main thread initially
                expects all iterators to be non-empty and will wait on an empty queue if it is expecting more iterators
                to be enqueued eventually, here we decrement to signal that it should not wait for this iterator since
                it was empty
                 */
                this.sourcesWaiting.getAndDecrement();
            } else {
                try {
                    this.queue.put(iterator);
                } catch (final InterruptedException e) {
                    throw new GATKException("Interrupted exception while enqueueing", e);
                }
            }
        }
    }
}

