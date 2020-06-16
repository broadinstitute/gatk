package org.broadinstitute.hellbender.engine;

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
import java.net.URI;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

public class ReadsHtsgetDataSource implements ReadsDataSource {
    private static final Logger logger = LogManager.getLogger(ReadsPathDataSource.class);

    private final List<GATKPath> sources;
    private List<SimpleInterval> intervals;
    private boolean traverseUnmapped;

    private final Map<GATKPath, SAMFileHeader> headers;
    private final SamFileHeaderMerger samHeaderMerger;
    private final SAMFileHeader header;

    private final SamReaderFactory readerFactory;

    private CloseableIterator<GATKRead> iterator;

    public ReadsHtsgetDataSource(final GATKPath source) {
        this(source, null);
    }

    public ReadsHtsgetDataSource(final GATKPath source, final SamReaderFactory customSamReaderFactory) {
        this(Collections.singletonList(source), customSamReaderFactory);
    }

    public ReadsHtsgetDataSource(final List<GATKPath> sources) {
        this(sources, null);
    }

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
                final URI sourceURI = source.getURI();
                final HtsgetRequest req = new HtsgetRequest(source)
                    .withDataClass(HtsgetClass.header);

                // Request only the headers and use them to construct SAMFileHeader for each source
                final InputStream headerStream = new SequenceInputStream(Collections.enumeration(
                    req.getResponse().streamData().collect(Collectors.toList())));
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

    public void setTraversalBounds(final List<SimpleInterval> intervals, final boolean traverseUnmapped) {
        this.intervals = intervals != null && ! intervals.isEmpty() ? intervals : null;
        this.traverseUnmapped = traverseUnmapped;
    }

    public boolean traversalIsBounded() {
        return this.intervals != null || this.traverseUnmapped;
    }

    /**
     * This data source can be queried even without index files
     * @return always true
     */
    public boolean isQueryableByInterval() {
        return true;
    }

    @Override
    @Nonnull
    public Iterator<GATKRead> iterator() {
        this.iterator = new MergingHtsgetIterator(this.intervals, this.traverseUnmapped);
        return this.iterator;
    }

    @Override
    public Iterator<GATKRead> query(final SimpleInterval interval) {
        this.iterator = new MergingHtsgetIterator(Collections.singletonList(interval), false);
        return this.iterator;
    }

    public Iterator<GATKRead> queryUnmapped() {
        this.iterator = new MergingHtsgetIterator(null, true);
        return this.iterator;
    }

    public boolean supportsSerialIteration() {
        return false;
    }

    public void close() {
        this.iterator.close();
    }

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
     * Get the sequence dictionary for this ReadsHtsgetDataSource
     *
     * @return SAMSequenceDictionary from the SAMReader backing this if there is only 1 input file, otherwise the merged SAMSequenceDictionary from the merged header
     */
    @Override
    public SAMSequenceDictionary getSequenceDictionary() {
        return this.header.getSequenceDictionary();
    }


    private class MergingHtsgetIterator implements CloseableIterator<GATKRead> {
        private final List<SimpleInterval> intervals;
        private final boolean traverseUnmapped;
        private final ExecutorService executor;
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

        private CloseableIterator<GATKRead> currentIterator;

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

        private void startNoIntervals() {
            for (final GATKPath source : sources) {
                this.sourcesWaiting.getAndIncrement();
                // Enqueue an iterator over either the entire source or just its unmapped reads
                this.executor.submit(this.traverseUnmapped
                    ? () -> this.getSourceUnmapped(source)
                    : () -> this.getSource(source));
            }
        }

        private List<SimpleInterval> intervalsOverlappingSource(final GATKPath source) {
            final SAMSequenceDictionary sequenceDictionary = headers.get(source).getSequenceDictionary();
            return this.intervals.stream()
                .filter(interval -> IntervalUtils.intervalIsOnDictionaryContig(interval, sequenceDictionary))
                .collect(Collectors.toList());
        }

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
            final HtsgetRequest req = new HtsgetRequest(source).withInterval(new SimpleInterval("*"));
            final Path tempPath = IOUtils.createTempPath("htsget-temp", "");
            try (final OutputStream ostream = Files.newOutputStream(tempPath)) {
                req.getResponse().streamData().forEachOrdered(stream -> {
                    try {
                        org.apache.commons.io.IOUtils.copy(stream, ostream);
                    } catch (final IOException e) {
                        throw new GATKException("Error while downloading block", e);
                    }
                });
            } catch (final IOException e) {
                throw new UserException(source.toString(), e);
            }
            final SamReader samReader = readerFactory.open(tempPath);
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
            final HtsgetRequest req = new HtsgetRequest(source).withInterval(new SimpleInterval("*"));
            final Path tempPath = IOUtils.createTempPath("htsget-temp", "");
            try (final OutputStream ostream = Files.newOutputStream(tempPath)) {
                req.getResponse().streamData().forEachOrdered(stream -> {
                    try {
                        org.apache.commons.io.IOUtils.copy(stream, ostream);
                    } catch (final IOException e) {
                        throw new GATKException("Error while downloading block", e);
                    }
                });
            } catch (final IOException e) {
                throw new UserException(source.toString(), e);
            }
            final SamReader samReader = readerFactory.open(tempPath);
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
            final Path tempPath = IOUtils.createTempPath("htsget-temp", "");
            try (final OutputStream ostream = Files.newOutputStream(tempPath)) {
                req.getResponse().streamData().forEachOrdered(stream -> {
                    try {
                        org.apache.commons.io.IOUtils.copy(stream, ostream);
                    } catch (final IOException e) {
                        throw new GATKException("Error while downloading block", e);
                    }
                });
            } catch (final IOException e) {
                throw new UserException(source.toString(), e);
            }

            final SamReader samReader = readerFactory.open(tempPath);
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

